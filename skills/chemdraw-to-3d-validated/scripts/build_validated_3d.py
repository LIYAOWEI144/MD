#!/usr/bin/env python
"""Build a deterministic 3D ligand and verify its chemistry against a 2D source."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--smiles", help="Confirmed isomeric SMILES.")
    source.add_argument("--input", type=Path, help="Approved 2D SDF, MOL, or MOL2 file.")
    parser.add_argument("--name", default="LIG", help="Residue and output molecule name.")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20261009)
    parser.add_argument("--conformers", type=int, default=20)
    return parser.parse_args()


def load_source(args: argparse.Namespace) -> Chem.Mol:
    if args.smiles:
        mol = Chem.MolFromSmiles(args.smiles)
    else:
        suffix = args.input.suffix.lower()
        if suffix in {".sdf", ".sd"}:
            supplier = Chem.SDMolSupplier(str(args.input), removeHs=False)
            mol = next((item for item in supplier if item is not None), None)
        elif suffix == ".mol2":
            mol = Chem.MolFromMol2File(str(args.input), removeHs=False)
        else:
            mol = Chem.MolFromMolFile(str(args.input), removeHs=False)
    if mol is None:
        raise ValueError("Could not read a valid molecular graph from the supplied source.")
    Chem.SanitizeMol(mol)
    return Chem.RemoveHs(mol)


def graph_signature(mol: Chem.Mol) -> dict:
    atoms = [
        {
            "index": atom.GetIdx(),
            "element": atom.GetSymbol(),
            "formal_charge": atom.GetFormalCharge(),
            "is_aromatic": atom.GetIsAromatic(),
            "chiral_tag": str(atom.GetChiralTag()),
        }
        for atom in mol.GetAtoms()
    ]
    bonds = [
        {
            "begin": min(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            "end": max(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()),
            "bond_type": str(bond.GetBondType()),
            "is_aromatic": bond.GetIsAromatic(),
            "stereo": str(bond.GetStereo()),
        }
        for bond in mol.GetBonds()
    ]
    bonds.sort(key=lambda item: (item["begin"], item["end"]))
    return {"atoms": atoms, "bonds": bonds}


def isomeric_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True, canonical=True)


def inchi_key(mol: Chem.Mol) -> str:
    molecule = Chem.RemoveHs(Chem.Mol(mol))
    molecule.RemoveAllConformers()
    for atom in molecule.GetAtoms():
        atom.SetAtomMapNum(0)
    return Chem.MolToInchiKey(molecule)


def build_3d(source: Chem.Mol, seed: int, conformers: int) -> Chem.Mol:
    molecule = Chem.AddHs(Chem.Mol(source), addCoords=True)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.25
    ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=conformers, params=params))
    if not ids:
        raise RuntimeError("RDKit could not embed a 3D conformer for this molecular graph.")
    energies = []
    for conformer_id in ids:
        try:
            AllChem.MMFFOptimizeMolecule(molecule, confId=conformer_id)
            properties = AllChem.MMFFGetMoleculeProperties(molecule)
            forcefield = AllChem.MMFFGetMoleculeForceField(molecule, properties, confId=conformer_id)
        except Exception:
            AllChem.UFFOptimizeMolecule(molecule, confId=conformer_id)
            forcefield = AllChem.UFFGetMoleculeForceField(molecule, confId=conformer_id)
        energies.append((forcefield.CalcEnergy(), conformer_id))
    best = min(energies)[1]
    for conformer in list(molecule.GetConformers()):
        if conformer.GetId() != best:
            molecule.RemoveConformer(conformer.GetId())
    return molecule


def main() -> int:
    args = parse_args()
    source = load_source(args)
    for index, atom in enumerate(source.GetAtoms(), start=1):
        atom.SetAtomMapNum(index)
    product = build_3d(source, args.seed, args.conformers)
    for index, atom in enumerate(source.GetAtoms()):
        product.GetAtomWithIdx(index).SetChiralTag(atom.GetChiralTag())
    product.SetProp("_Name", args.name)
    product.SetProp("residue_name", args.name)
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    source_smiles = isomeric_smiles(source)
    product_heavy = Chem.RemoveHs(product)
    product_smiles = isomeric_smiles(product_heavy)
    source_inchi = inchi_key(source)
    product_inchi = inchi_key(product_heavy)
    report = {
        "name": args.name,
        "source_isomeric_smiles": source_smiles,
        "product_isomeric_smiles": product_smiles,
        "source_inchi_key": source_inchi,
        "product_inchi_key": product_inchi,
        "formula": rdMolDescriptors.CalcMolFormula(source),
        "formal_charge": Chem.GetFormalCharge(source),
        "heavy_atom_count": source.GetNumAtoms(),
        "atom_count_with_hydrogen": product.GetNumAtoms(),
        "molecular_weight": Descriptors.MolWt(product),
        "graph_match": graph_signature(source) == graph_signature(product_heavy),
        "inchi_key_match": source_inchi == product_inchi,
        "stereochemistry_match": source_smiles == product_smiles,
        "seed": args.seed,
    }
    report["validated"] = all(
        report[key] for key in ("graph_match", "inchi_key_match", "stereochemistry_match")
    )

    writer = Chem.SDWriter(str(output / f"{args.name}_3d.sdf"))
    writer.write(product)
    writer.close()
    Chem.MolToPDBFile(product, str(output / f"{args.name}_3d.pdb"))
    with (output / "atom_map.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["source_heavy_atom_index_1based", "element", "product_heavy_atom_index_1based"])
        for atom in source.GetAtoms():
            writer.writerow([atom.GetIdx() + 1, atom.GetSymbol(), atom.GetIdx() + 1])
    (output / "validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["validated"] else 2


if __name__ == "__main__":
    sys.exit(main())
