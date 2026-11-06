# solve

## Purpose
Iterative solvers and preconditioners. `Multigrid` is a geometric multigrid solver on a hierarchy of refined meshes with the smoother 'Richardson' or 'Jacobi' (character vectors) or a user-defined function handle, a damping factor, M pre- and M post-smoothing steps per level and W coarse grid corrections per level and cycle (W = 1 gives a V-cycle, W = 2 a W-cycle); its method solve stops when the residual norm is at most the absolute tolerance TOL. The `Schwarz` classes are additive and multiplicative Schwarz preconditioners built on the abstract `Schwarz` class. `pdaSet` is a primal-dual active set algorithm for min 1/2 x'Ax - b'x subject to the bound constraints low <= x <= up, as in obstacle problems, and `pinvit` an inexact preconditioned inverse iteration for one eigenpair (usually the smallest eigenvalue) of the symmetric positive definite pair (A, B), with cycles of a `Multigrid` solver as preconditioner.

## Key concepts
- Mesh hierarchy: the levels are added with `addLevel` from the coarsest to the finest mesh, each with its system matrix and the prolongation from the previous level (`fe.prolongation(feOld)`, see `FE` in `elements` and the refinement data of `mesh`).
- Matrix conventions: `Multigrid.addLevel` takes the matrix and the prolongation already restricted to the free degrees of freedom, `A(fe.freeNodes,fe.freeNodes)` and `P(fe.freeNodes,feOld.freeNodes)`, or full ones together with the mask `fixedNodes`. The `addLevel` of the `Schwarz` classes takes the spaces FE and FEOLD, the full prolongation and the full matrix and restricts them to `fe.freeNodes` itself; the preconditioner then acts on the free degrees of freedom of the finest level, e.g. `@(r)applyMult(S,r)` in `pcg`.
- Domain decomposition preconditioners: `SchwarzBPS` (Bramble-Pasciak-Schatz), `SchwarzBDDC` (balancing domain decomposition by constraints), `SchwarzNN` (Neumann-Neumann). They need the side dofs `dof4ed` (in 3d `dof4f`) and reject elements without them, such as the discontinuous ones and `HPk3d`. `SchwarzBDDC` and `SchwarzNN` also need the local element matrices A4E (the second output of `assemble` of `BilinearForm`). All three assume the vertex dofs (for `SchwarzBDDC` all primal dofs) numbered first and the element interior dofs last, as `dofManager` numbers elements with one component (`Pk3d` of degree >= 3 reorients its edge dofs and does not fit `SchwarzBDDC`); `SchwarzBPS` instead uses `dof4n` and the interior dof columns `idof4ed`, `idof4f`, `idof4e` where the element has them (`PkVector`).
- Multilevel preconditioners: `SchwarzBPX` (Bramble-Pasciak-Xu), `SchwarzHB` (hierarchical basis, no coarse solve). `SchwarzHB` needs nested numbering: the free dofs of the previous level keep their numbers, as for `Pk` of degree 1 on meshes refined with `refine`.
- Vertex patch preconditioner: `SchwarzVO` (overlapping patches around the vertices; for elements without `dof4ed` (in 3d `dof4f`), such as the discontinuous ones and `HPk3d`, the patches hold all degrees of freedom of the patch elements).
- Local solves: with the property `pseudoinverse` set before `addLevel`, `SchwarzBDDC`, `SchwarzBPS` and `SchwarzVO` invert the local matrices with `pinv`; `SchwarzNN` ignores it.

## Relation to other parts
The system matrices are assembled with the `forms` classes for element spaces from `elements` on meshes from `mesh`. In the examples, `Multigrid` solves Poisson problems directly, the `Schwarz` classes precondition the CG method (`pcg`), and `pdaSet` solves the discrete obstacle problems of `exObstacleP1Adaptive`.

## Contents

{{MEMBERS}}

## Related examples

{{RELATED_EXAMPLES}}
