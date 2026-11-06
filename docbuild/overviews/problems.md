# problems

## Purpose
Ready-made problem definitions, grouped by PDE type: scalar elliptic (Poisson), biharmonic and plate, elasticity, Stokes, Maxwell, Helmholtz, eigenvalue, obstacle, transport, vector elliptic and surface problems. A problem object provides the geometry (`geom`), the coefficients and data of the PDE (e.g. `f`, `u4Db`, `g`) and, for the problems with a known solution (most have Exact in the name), the exact solution (e.g. `exactU`, `exactGradU`). Eigenvalue problems provide `geom` and reference eigenvalues `eigenvalue(m)` instead of PDE data, with NaN for a mode without a reference value.

## Key concepts
- Templates and concrete problems: most subfolders contain an abstract template such as `ProblemScalarEllipticTemplate` (and `ProblemScalarEllipticTemplateExact` for problems with exact solution) that defines the interface and the array layouts of the data; the `problem*` classes implement it for a domain, coefficients and boundary conditions. The `transport` and `vectorelliptic` problems derive from the scalar elliptic templates; the `maxwell` problems define their data directly, without a template.
- Combined templates: the obstacle problems combine a scalar elliptic template with `ProblemObstacleTemplate`, which adds the bounds `lower` and `upper` (defaults -inf and inf). In `surface`, `ProblemSurfaceTemplate` describes a surface in 3d by a level-set function; classes such as `problemSurfaceSphere` define only the surface; most surface PDE problems combine a surface class or `ProblemSurfaceTemplate` with a PDE template (`ProblemScalarEllipticTemplateExact`, `ProblemStokesTemplate`, `ProblemStokesTemplateExact`), while the vector Laplace problems define their data directly. The eigenvalue problems derive from `ProblemEigenvalue`.
- The domain is a `Geometry` (from `geometries`), selected by name and options, e.g. `Geometry(struct('name','Lshape'))`.
- Parameter classes: `ElasticityParameters` (material parameters E, nu, lambda, mu and the elasticity tensors) and `StokesParameters` (viscosity nu and compliance matrices). `ProblemElasticityTemplate`, `ProblemStokesTemplate` and the elasticity eigenvalue problems derive from them, so every elasticity problem needs the material parameters in its constructor, e.g. `problemElasticitySquareExact('E',1,'nu',0.3)` or `('lambda',1,'mu',1)`, and a Stokes problem takes the viscosity (default 1).
- Problems with exact solutions are used to measure convergence rates of discretizations.

## Relation to other parts
A script creates a mesh from `p.geom` (from `mesh`; for surface problems a `MeshSurface` that also takes `@p.surfaceFunction` and `@p.gradSurfaceFunction`), chooses an element space (from `elements`) and passes the data of the problem as function handles to the `assemble` methods of the `forms` classes, e.g. `b.assemble(fe,'f',@p.f)`.

## Contents

{{SUBFOLDERS}}

## Related examples

{{RELATED_EXAMPLES}}
