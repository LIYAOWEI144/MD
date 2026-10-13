#!/usr/bin/env python
"""Build a deterministic 3D ligand and verify its chemistry against a 2D source."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import xml.etree.ElementTree as ElementTree
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--smiles", help="Confirmed isomeric SMILES.")
    source.add_argument("--input", type=Path, help="Approved 2D SDF, MOL, or MOL2 file.")
    parser.add_argument(
        "--expected-smiles",
        help="Independent, confirmed isomeric SMILES used to audit a CDXML/SDF/MOL source before 3D generation.",
    )
    parser.add_argument("--expected-formula", help="Optional independently recorded molecular formula.")
    parser.add_argument("--expected-charge", type=int, help="Optional independently recorded total formal charge.")
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
        if suffix == ".cdxml":
            mol = load_cdxml(args.input)
        elif suffix in {".sdf", ".sd"}:
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


def load_cdxml(path: Path) -> Chem.Mol:
    """Use RDKit's CDXML reader, which handles both wedge ends and nested fragments."""
    molecules = Chem.MolsFromCDXMLFile(str(path), removeHs=False)
    if len(molecules) != 1 or molecules[0] is None:
        raise ValueError("The CDXML source must contain exactly one valid molecule.")
    return molecules[0]


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


def identity_summary(mol: Chem.Mol) -> dict:
    """Return identity fields that catch graph, charge, and stereochemical drift."""
    clean = Chem.RemoveHs(Chem.Mol(mol))
    Chem.AssignStereochemistry(clean, cleanIt=True, force=True)
    return {
        "isomeric_smiles": isomeric_smiles(clean),
        "inchi_key": inchi_key(clean),
        "formula": rdMolDescriptors.CalcMolFormula(clean),
        "formal_charge": Chem.GetFormalCharge(clean),
        "heavy_atom_count": clean.GetNumAtoms(),
        "graph": graph_signature(clean),
    }


def compare_identity(reference: Chem.Mol, candidate: Chem.Mol) -> dict:
    left = identity_summary(reference)
    right = identity_summary(candidate)
    return {
        "graph_match": left["graph"] == right["graph"],
        "inchi_key_match": left["inchi_key"] == right["inchi_key"],
        "stereochemistry_match": left["isomeric_smiles"] == right["isomeric_smiles"],
        "formula_match": left["formula"] == right["formula"],
        "formal_charge_match": left["formal_charge"] == right["formal_charge"],
    }


def build_3d(source: Chem.Mol, seed: int, conformers: int) -> Chem.Mol:
    molecule = Chem.AddHs(Chem.Mol(source), addCoords=True)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.25
    params.maxIterations = 4000
    ids = list(AllChem.EmbedMultipleConfs(molecule, numConfs=conformers, params=params))
    if not ids:
        params.useRandomCoords = True
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
    Chem.AssignStereochemistry(source, cleanIt=True, force=True)
    expected = Chem.MolFromSmiles(args.expected_smiles) if args.expected_smiles else None
    if args.input and args.input.suffix.lower() == ".cdxml" and expected is None:
        raise ValueError("CDXML input requires --expected-smiles to independently confirm its parsed 2D graph.")
    if expected is not None:
        Chem.SanitizeMol(expected)
        Chem.AssignStereochemistry(expected, cleanIt=True, force=True)
        source_reference = compare_identity(expected, source)
    else:
        source_reference = {"status": "approved_source_used_as_reference"}
    if args.expected_formula and identity_summary(source)["formula"] != args.expected_formula:
        raise ValueError("Source formula differs from --expected-formula; stop before 3D generation.")
    if args.expected_charge is not None and Chem.GetFormalCharge(source) != args.expected_charge:
        raise ValueError("Source formal charge differs from --expected-charge; stop before 3D generation.")
    if expected is not None and not all(source_reference.values()):
        raise ValueError("Parsed 2D source differs from --expected-smiles; stop before 3D generation.")
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
    product_check = compare_identity(source, product_heavy)
    source_identity = identity_summary(source)
    product_identity = identity_summary(product_heavy)
    report = {
        "name": args.name,
        "source_isomeric_smiles": source_smiles,
        "product_isomeric_smiles": product_identity["isomeric_smiles"],
        "source_inchi_key": source_identity["inchi_key"],
        "product_inchi_key": product_identity["inchi_key"],
        "formula": source_identity["formula"],
        "formal_charge": source_identity["formal_charge"],
        "heavy_atom_count": source_identity["heavy_atom_count"],
        "atom_count_with_hydrogen": product.GetNumAtoms(),
        "molecular_weight": Descriptors.MolWt(product),
        "source_reference_check": source_reference,
        "product_check": product_check,
        "graph_match": product_check["graph_match"],
        "inchi_key_match": product_check["inchi_key_match"],
        "stereochemistry_match": product_check["stereochemistry_match"],
        "formula_match": product_check["formula_match"],
        "formal_charge_match": product_check["formal_charge_match"],
        "seed": args.seed,
    }
    report["validated"] = all(
        report[key]
        for key in (
            "graph_match", "inchi_key_match", "stereochemistry_match",
            "formula_match", "formal_charge_match",
        )
    )

    source_copy = Chem.Mol(source)
    AllChem.Compute2DCoords(source_copy)
    writer = Chem.SDWriter(str(output / f"{args.name}_approved_2d.sdf"))
    writer.write(source_copy)
    writer.close()
    writer = Chem.SDWriter(str(output / f"{args.name}_3d.sdf"))
    writer.write(product)
    writer.close()
    Chem.MolToPDBFile(product, str(output / f"{args.name}_3d.pdb"))
    with (output / "atom_map.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "source_heavy_atom_index_1based", "element", "formal_charge", "aromatic",
            "product_heavy_atom_index_1based",
        ])
        for atom in source.GetAtoms():
            writer.writerow([
                atom.GetIdx() + 1, atom.GetSymbol(), atom.GetFormalCharge(),
                atom.GetIsAromatic(), atom.GetIdx() + 1,
            ])
    with (output / "bond_audit.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["begin_1based", "end_1based", "bond_type", "aromatic", "stereo"])
        for bond in source.GetBonds():
            writer.writerow([
                bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1,
                str(bond.GetBondType()), bond.GetIsAromatic(), str(bond.GetStereo()),
            ])
    (output / "validation.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["validated"] else 2


if __name__ == "__main__":
    sys.exit(main())
