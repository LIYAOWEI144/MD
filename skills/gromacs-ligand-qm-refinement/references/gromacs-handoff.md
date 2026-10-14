# GROMACS Handoff Requirements

## Include Closure

The delivered `LIG.top`/`LIG.itp` must identify the exact CHARMM/CGenFF
force-field source used for parameterization. Do not repair a missing include
by substituting a stock GROMACS force field. Before complex assembly, resolve
all quoted includes and confirm that the receptor and ligand share one
compatible nonbonded convention.

## Position Restraints

- Do not append a `[ position_restraints ]` block that is activated by a
  generic `POSRES` macro unless the complete simulation protocol explicitly
  defines numerical force constants for it.
- Prefer `POSRES_LIG` for a ligand-specific restraint. It must be included
  conditionally and should be inactive by default.
- A generated `posre_lig.itp` is not sufficient by itself; record the matching
  `LIG.itp` include and required MDP `define` statement.

## Solvated-System Handoff

When a user requests subsequent MD preparation, validate the assembled system
after each major stage: dry complex, solvated system, ionized system, EM, NVT,
and NPT. Water and ion molecule types must match coordinate residue names and
the selected force field. A successful dry-complex TPR alone does not prove
that `SOL`, `NA`, or `CL` will be recognized later.

## Production Readiness

Use a GROMACS preflight without `-maxwarn` before production. Preserve any
fast-bond stability warning and resolve it with an explicit timestep or
constraint decision. For normal structural analyses, plan compressed
trajectory output and avoid unnecessary full-precision coordinate and velocity
trajectories. The analysis-assistant skill is responsible for producing the
user-facing command review and storage plan.
