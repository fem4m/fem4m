# Examples

The examples are tutorials. Each page states the problem and its weak formulation, derives the discretization and the linear system, and relates every integral to its fem4m form string. Then it shows the commented script with its output (convergence tables and figures), compares the observed with the expected convergence, and ends with possibilities for extensions. The scripts are plain-text live scripts in `toolbox/examples`: open one in MATLAB to run it and to change its parameters.

### Learning path
Start with the basics, then follow the branches you are interested in. Every example names the examples it builds on and the ones that continue it.

- Basics: `exMeshEnumeration`, `exElementBasisFunctions`, `exPoissonP1`
- Adaptivity: `exPoissonPkAdaptive`, then `exPoissonHPAdaptive`, `exEllipticCRAdaptive`, `exPoisson3dPkAdaptive`, `exLaplaceBeltramiPkAdaptive`, `exEigenvaluePkAdaptive` (and `exEigenvalueDiffusionPkAdaptive`), `exElasticityPkAdaptive`, `exHelmholtzPkAdaptive`, `exObstacleP1Adaptive`
- Mixed methods: `exPoissonMixedP1P0Unstable`, `exPoissonRTAdaptive`, then `exElasticityAW`, `exBiharmonicHHJ`, `exStokesTaylorHood`
- Discontinuous Galerkin methods: `exPoissonDG`, then `exPoissonDGAdaptive`, `exPoissonLDG`, `exPoissonMixedDG` (and `exPoissonHDG`), `exStokesDG`, `exTransportDG`
- Conforming elements for other spaces: `exBiharmonicArgyris` (C1), `exMaxwellNd` (H(rot))
- Solvers: `exSolverSmoothing`, `exSolverTwoGrid`, `exSolverMultigrid`, `exSolverSchwarz`

### Notation
- The triangulation of the domain Ω consists of the triangles T (tetrahedra in 3D) with the edges E (faces in 3D); an interior edge lies inside Ω.
- h_T and h_E are the diameters of T and E, n_E is a unit normal of E, and [v] and {v} are the jump and the average of v across an interior edge.
- P_k are the piecewise polynomials of degree k on the triangulation, V_h is a discrete space and u_h a discrete solution.
- Γ_D and Γ_N are the Dirichlet and the Neumann boundary with the data u_D and g.
- N is the number of degrees of freedom. The convergence rates are rates in N; in 2D, N is about h^-2 on uniform meshes, so a rate 1/2 in N is order 1 in h.
- η is an error estimator and η_T its contribution of the element T.

## Basics

- exMeshEnumeration
- exElementBasisFunctions
- exPoissonP1

## Poisson and scalar elliptic problems: adaptive methods

- exPoissonPkAdaptive
- exPoissonHPAdaptive
- exEllipticCRAdaptive

## Poisson problem: mixed methods

- exPoissonMixedP1P0Unstable
- exPoissonRTAdaptive

## Poisson problem: discontinuous Galerkin methods

- exPoissonDG
- exPoissonDGAdaptive
- exPoissonLDG
- exPoissonMixedDG
- exPoissonHDG

## 3D and surface problems

- exPoisson3dPkAdaptive
- exLaplaceBeltramiPkAdaptive

## Eigenvalue problems

- exEigenvaluePkAdaptive
- exEigenvalueDiffusionPkAdaptive

## Linear elasticity

- exElasticityPkAdaptive
- exElasticityAW

## Stokes problem

- exStokesTaylorHood
- exStokesDG

## Biharmonic problem

- exBiharmonicArgyris
- exBiharmonicHHJ

## Maxwell and Helmholtz problems

- exMaxwellNd
- exHelmholtzPkAdaptive

## Obstacle problem

- exObstacleP1Adaptive

## Transport equation

- exTransportDG

## Solvers and preconditioners

- exSolverSmoothing
- exSolverTwoGrid
- exSolverMultigrid
- exSolverSchwarz
