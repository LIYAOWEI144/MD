---
name: gromacs-ligand-qm-refinement
description: Prepare a user-supplied 3D small-molecule ligand for receptor-ligand GROMACS analysis, including chemistry validation, CGenFF/CHARMM parameter review, targeted QM refinement, and delivery of reproducible GRO/ITP/TOP files. Use when a ligand needs parameterization or repair before GROMACS MD; do not use for proteins alone.
---

# GROMACS Ligand QM Refinement

## Outcome

Turn a supplied 3D ligand into a reproducible, reviewable GROMACS delivery
package. The final output is not a Gaussian log: it must include a ligand
`LIG.itp`, matching `LIG.gro`, ligand `LIG.top`, and, when receptor coordinates
are supplied, a receptor-ligand `.gro`, receptor `.itp`, and complex `.top`.
Create these in a new delivery folder inside the user's working directory.
Never overwrite the input topology or a global force-field file.

## First Pass

1. Identify the input ligand, receptor, GROMACS version, force-field family,
   and whether the requested system is ligand-only or a complex.
2. Preserve the original structure and generate a chemistry-valid source
   `SDF` or `MOL2` with explicit hydrogens, bond orders, formal charge, and a
   declared protonation/tautomer state. A PDB alone is insufficient evidence
   for a novel ligand's bond orders. If the 3D ligand originated from a 2D
   drawing/CDXML, require an independent approved 2D identity and a passing
   2D-to-3D audit before parameterization.
3. Atom-map the source structure to the eventual topology. Verify element,
   formal charge, aromaticity, bond order, stereo, and total charge against
   the approved 2D record before replacing charges or bonded terms.
4. Use a ligand-specific parameterizer compatible with the selected force
   field, such as CGenFF/CHARMM-GUI for a CHARMM workflow. Do not run
   `pdb2gmx` on a nonstandard ligand.

## Parameter Triage

Read [the detailed workflow](references/workflow.md) when performing the
parameterization. In particular:

- inspect CGenFF penalties and map every high-penalty term back to ligand atom
  numbers and chemistry;
- treat aromatic/fused ring terms as coupled, not independent rotors;
- use QM water-interaction fitting only for clearly scoped charge refinement;
- use constrained QM scans for torsions and a compatible bonded-parameter
  optimizer for true Hessian fitting;
- do not copy a projected QM Hessian element directly into a bond or angle
  parameter for a fused ring.
- stop before CGenFF/QM work if `validation.json`, `atom_map.csv`, or
  `bond_audit.csv` identifies a 2D-to-3D discrepancy.

## Delivery and Validation

Create a folder such as `GROMACS_delivery_<ligand>` containing standard names:

```text
LIG.itp
LIG.gro
LIG.top
forcefield.itp
PROA.itp                 # when a receptor is supplied
PROA_LIG.gro             # when a receptor is supplied
PROA_LIG.top             # when a receptor is supplied
reproducible_scripts/
README.md
```

The `reproducible_scripts` folder must retain the structure conversion,
atom-map, QM-input generation, fitting, and packaging scripts actually used.
Put QM inputs/outputs and fit reports in a separate `QM_refinement/` folder so
the deliverable stays usable while provenance remains inspectable.

Run a ligand-only `gmx grompp` check. For a complex, also run full-complex
`gmx grompp`; if it fails, run the same command against the unmodified control
topology. Report whether the defect predates the refined ligand. Do not claim
the complex is production-ready without a generated complex `.tpr`.

## MD Handoff Contract

Read [GROMACS handoff](references/gromacs-handoff.md) before preparing files
for solvation, equilibration, or production. The delivery must preserve a
complete compatible include chain and must not activate a ligand restraint
block through a generic protein `POSRES` macro by accident. When a user
supplies an MD command list, hand it to `gromacs-analysis-assistant` for a
Chinese command review before declaring the package ready for manual use.

State in `README.md`: input identity and protonation state, force-field and
tool versions, parameter changes, QM level of theory, known limitations, exact
validation commands, and which outputs have actually passed validation.
