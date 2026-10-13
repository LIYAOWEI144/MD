# Detailed Workflow

## Structure and Chemical State

1. Inspect the input structure before editing. Record source path, residue
   name, atom count, formal charge, disconnected fragments, and any covalent
   receptor attachment.
2. Establish an approved 2D chemical source of truth before trusting a 3D
   structure. For drawings/CDXML, require an independently reviewed isomeric
   SMILES rather than a SMILES derived from the same parser output. Re-export
   a MOL2 or SDF from a chemistry-aware representation, preserve
   stereochemistry, and create an atom map between source atom labels and
   topology indices. Use explicit hydrogens. Do not infer an arbitrary
   protonation state from a docking PDB; identify the intended pH/chemical
   state and record the assumption.
3. Before CGenFF, compare the approved 2D record and 3D record for formula,
   total formal charge, per-atom element/charge/aromaticity, bond order,
   double-bond stereo, tetrahedral stereo, and heavy-atom mapping. Save the
   checked `validation.json`, `atom_map.csv`, and `bond_audit.csv`. Stop if
   any field differs; never repair a discrepancy only by renaming topology
   atoms.
4. For ligands with unusual sulfur, phosphorus, metals, radicals, or covalent
   linkages, stop automated parameterization if the selected tool reports an
   unsupported atom type. Repair the chemical state or choose an appropriate
   specialized method before continuing.

## Initial CHARMM/CGenFF Parameterization

1. Generate the CGenFF/CHARMM-compatible topology from the validated SDF/MOL2.
   Save all returned RTF/PRM/STR or GROMACS conversion files.
2. Convert to GROMACS using a tool/version compatible with the intended
   CHARMM force field. Keep the original generated `LIG.itp` immutable.
3. Extract all penalty scores and classify them by charge, bond, angle, and
   dihedral. Rank terms by penalty, then inspect whether each term appears in
   a flexible linker, aromatic system, fused ring, or unusual functional group.

## QM Refinement

Use a local Gaussian installation only when it is configured by the user and
the calculation is authorized. Create a job folder with `.gjf`, `.chk`, `.log`,
and `.fchk` files; verify normal termination before fitting.

### Charges

- Build chemically meaningful water-probe complexes at donor/acceptor sites.
- Optimize the chosen QM target and fit only the designated atom charges while
  preserving the total ligand charge and charge symmetry where chemically
  justified.
- Record before/after charges and the fitting objective. Validate that topology
  atom order is unchanged.

### Torsions

- Derive scan coordinates from the topology atom map, not atom names alone.
- For a normal rotatable bond, use constrained relaxed scans spanning an
  adequate angular range. For fused/aromatic rings, define a coupled distortion
  and explicitly verify the final constrained dihedrals from the formatted
  checkpoint; Gaussian can otherwise silently retain the planar input geometry.
- Fit a CHARMM-compatible periodic form to completed QM relative energies.
  Replace only the matching local topology entries with explicit parameters;
  never edit generic atom-type terms globally for a ligand-specific correction.
- Exclude failed, unconverged, or constraint-ignored scans. State the excluded
  points and why.

### Bonds and Angles

- A QM Hessian provides a reference geometry and normal-mode target. In coupled
  rings, projected internal-coordinate force constants are diagnostics, not
  automatically transferable independent force constants.
- Use a true bonded-parameter optimizer (for example CHARMM ffTK or a
  compatible Hessian/normal-mode fitting workflow) before replacing bonds or
  angles. If no compatible optimizer is available, retain the initial CGenFF
  terms and explicitly document that only directly supported parameters changed.

## GROMACS Packaging

1. Construct the ligand-only `LIG.top` and run `gmx grompp` with `LIG.gro`.
2. For receptor-ligand systems, preserve the receptor topology and coordinate
   atom order. Create a new `PROA_LIG.top` that includes the refined `LIG.itp`.
   Do not include the original and refined ligand topologies together.
3. Create receptor-only and ligand-only PDB exports from the delivery GRO files
   when requested. Confirm their atom counts match the relevant topology blocks.
4. Run complex `gmx grompp`, energy minimization, and a short NVT/NPT only
   when the user has supplied the necessary solvent/ion system and MDP files.
   Never label an unsolvated preprocessor-only package as a production MD run.

## Required Report

The delivery README must list: source files, force-field source/version, total
charge, atom count, parameter changes, QM jobs used, validation result for
ligand and complex separately, and commands required to reproduce the package.
