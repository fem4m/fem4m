# fem4m

fem4m is a MATLAB toolbox for finite element methods. It provides finite element spaces (`elements`), variational forms written as strings and assembled into matrices and vectors (`forms`), mesh classes with adaptive refinement and marking strategies (`mesh`), ready-made problem definitions for many PDEs (`problems`), iterative solvers and preconditioners (`solve`), domain geometries (`geometries`) and helper functions (`internal`).

## Main parts
- **elements** — finite element spaces in 2d and 3d: Lagrange, DG and hp families, nonconforming, H(div) and H(curl) spaces, biharmonic (plate) and elasticity elements.
- **forms** — bilinear and linear forms, error and estimator forms, parsed from strings and assembled for an element space.
- **mesh** — triangular, tetrahedral, polygonal and surface meshes with uniform and adaptive refinement and marking strategies.
- **problems** — problem templates and concrete problems: scalar elliptic (Poisson), vector elliptic, biharmonic and plate, elasticity, Stokes, Maxwell, Helmholtz, eigenvalue, obstacle, transport and surface problems.
- **solve** — geometric multigrid, Schwarz domain decomposition and multilevel preconditioners, a primal-dual active set solver and a preconditioned inverse iteration for eigenvalue problems.
- **geometries** — 2d, 3d and surface domains, built in the code or loaded from data files, used by problems and meshes.
- **internal** — numerical, plotting and repository helpers used by the above.
- **apps** — the app `fem4m.mlapp` for adaptive computations with selectable problems, elements and solvers, and `plotApp.mlapp`, which plots saved results; these are not documented on these pages.

## Examples
The [Examples](examples/examples_overview.html) cover the Poisson problem (conforming, nonconforming, adaptive, hp-adaptive, mixed, DG, LDG and HDG methods, also in 3d), the Laplace-Beltrami problem on a surface, eigenvalue, elasticity, Stokes, biharmonic, Maxwell, Helmholtz, obstacle and transport problems, multigrid, smoothers and Schwarz preconditioners, and the handling of meshes and basis functions.
