# Structure Verification

Before converting a ChemDraw image into a simulation ligand, inspect the drawing and capture these facts in the task output:

1. Proposed isomeric SMILES and molecular formula.
2. Ring count and aromatic systems.
3. Functional groups and bond orders, including carbonyls, nitriles, sulfonamides, phosphates, and conjugated linkers.
4. Halogen identity and attachment positions.
5. Formal charge, intended protonation state, and counterions that are or are not part of the ligand.
6. Wedge/dash stereocenters and E/Z double-bond stereochemistry.

The validation script compares canonical isomeric SMILES, InChIKey, atom/bond graph, formula, formal charge, and heavy-atom stereo labels between its 2D source and 3D product. It also writes `atom_map.csv`; atom `n` in the approved heavy-atom source maps to atom `n` in the output before hydrogens are added.

If the molecule was read from a screenshot, visually compare the source drawing to a rendered 2D depiction of the approved SMILES. Stop for user confirmation if any listed property is ambiguous or inconsistent. A report that passes graph checks does not validate chemistry that was mistranscribed from the screenshot.
