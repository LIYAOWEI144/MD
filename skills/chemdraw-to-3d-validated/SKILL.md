---
name: chemdraw-to-3d-validated
description: Convert a user-supplied ChemDraw structure image or confirmed 2D chemical record into a validated 3D ligand while preserving a traceable atom order. Use for ligand preparation before docking, CGenFF, CHARMM-GUI, GROMACS, or OpenMM; do not use when the requested molecule identity is ambiguous.
---

# ChemDraw To 3D Validated

Create a reproducible 3D ligand from a ChemDraw-derived structure and prove that the generated molecule has the same molecular graph as the approved 2D source.

## Chemical Source of Truth

An image is a visual reference, not an authoritative machine-readable molecular graph. Inspect the ChemDraw image and extract a proposed SMILES only when every bond order, aromatic ring, formal charge, isotope, stereocenter, and attachment point is unambiguous. Present the proposed molecular identity for confirmation when any of these are unclear.

Prefer a user-provided ChemDraw export (`.cdx`, `.cdxml`), molfile, SDF, MOL2, or confirmed isomeric SMILES. Treat that approved record as the source of truth. Do not claim that a 3D structure is identical to an image until the visual and graph checks below have passed.

## Workflow

1. Create a clean 2D source record from the confirmed isomeric SMILES or exported SDF/MOL. Keep original heavy-atom order unless the user explicitly requests canonical renumbering. For CDXML or image-derived input, obtain an independently confirmed `--expected-smiles`; the parser's own output cannot validate itself.
2. Inspect and record a visual checklist: molecular formula, count of rings, key functional groups, halogens, formal charge, protonation/tautomer state, and all indicated stereochemistry. Compare it directly with the image and the approved 2D record.
3. Use `scripts/build_validated_3d.py` to embed and optimize conformers. It preserves heavy-atom order, appends hydrogens, and writes the approved 2D SDF, 3D SDF/PDB, atom mapping, bond audit, and a JSON validation report.
4. Accept output only when the report says `graph_match: true`, `inchi_key_match: true`, `stereochemistry_match: true`, `formula_match: true`, and `formal_charge_match: true`. For CDXML, the `source_reference_check` must also pass. Inspect the rendered 2D and 3D structures for obvious misparsed rings, bond orders, protonation changes, or inverted stereochemistry.
5. For force-field preparation, use the generated SDF/MOL2 only after the full audit passes. After CGenFF conversion, verify its atom order against `atom_map.csv` and compare all topology bonds to `bond_audit.csv` before merging coordinates with a protein.

## Commands

Use the dedicated local Conda environment:

```powershell
& 'D:\Miniforge\envs\chemdraw3d\python.exe' \
  'C:\Users\92336\.codex\skills\chemdraw-to-3d-validated\scripts\build_validated_3d.py' \
  --smiles 'CONFIRMED_ISOMERIC_SMILES' \
  --name LIG \
  --output-dir 'D:\path\to\ligand_3d'
```

For an approved ChemDraw `.cdxml`, 2D SDF, MOL, or MOL2 file, use `--input path\\to\\ligand.cdxml`. CDXML requires an independent reference identity:

```powershell
--input 'D:\path\to\ligand.cdxml' --expected-smiles 'CONFIRMED_ISOMERIC_SMILES'
```

Read `references/verification.md` before treating an image-derived structure as ready for docking or molecular dynamics.

## Boundaries

- Never invent bond orders, protonation states, stereochemistry, or formal charges from a low-resolution image.
- Do not use a 3D file as CGenFF/GROMACS input merely because embedding succeeded. An independently checked 2D chemical record and the generated 3D record must agree first.
- A valid graph check proves chemical connectivity, not that the chosen protonation state, tautomer, docking pose, or force-field parameters are scientifically optimal.
- CGenFF penalties still need review after topology generation, especially for new torsions or flexible, conjugated, and charged ligands.
