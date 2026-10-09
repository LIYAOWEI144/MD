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

1. Create a clean 2D source record from the confirmed isomeric SMILES or exported SDF/MOL. Keep original heavy-atom order unless the user explicitly requests canonical renumbering.
2. Inspect and record a visual checklist: molecular formula, count of rings, key functional groups, halogens, formal charge, and all indicated stereochemistry. Compare it directly with the image.
3. Use `scripts/build_validated_3d.py` to embed and optimize conformers. It preserves heavy-atom order, appends hydrogens, and writes SDF, PDB, atom mapping, and a JSON validation report.
4. Accept output only when the report says `graph_match: true`, `inchi_key_match: true`, and `stereochemistry_match: true`. Inspect the 3D structure in PyMOL or another molecular viewer for obvious clashes or inverted stereochemistry.
5. For force-field preparation, use the generated SDF/MOL2 as the CGenFF input. After CGenFF conversion, verify its atom order against `atom_map.csv` before merging coordinates with a protein.

## Commands

Use the dedicated local Conda environment:

```powershell
& 'D:\Miniforge\envs\chemdraw3d\python.exe' \
  'C:\Users\92336\.codex\skills\chemdraw-to-3d-validated\scripts\build_validated_3d.py' \
  --smiles 'CONFIRMED_ISOMERIC_SMILES' \
  --name LIG \
  --output-dir 'D:\path\to\ligand_3d'
```

For an approved 2D SDF or MOL file, use `--input path\\to\\ligand.sdf` instead of `--smiles`.

Read `references/verification.md` before treating an image-derived structure as ready for docking or molecular dynamics.

## Boundaries

- Never invent bond orders, protonation states, stereochemistry, or formal charges from a low-resolution image.
- A valid graph check proves chemical connectivity, not that the chosen protonation state, tautomer, docking pose, or force-field parameters are scientifically optimal.
- CGenFF penalties still need review after topology generation, especially for new torsions or flexible, conjugated, and charged ligands.
