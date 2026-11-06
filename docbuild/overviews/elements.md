# elements

## Purpose
This folder provides the finite element spaces of the toolbox. `FE` is the abstract base class of the element objects defined on mesh elements. It bundles the mesh and degree-of-freedom information with the evaluation routines for basis functions and finite element functions (`getValue`, `getGradValue`, `getDivValue`, `getD2Value` and the further get...Value routines of the derivatives), `interpolate`, `project`, `prolongation`, `prolongate` and `plot`.

`FEEdge` and `FEFace` are the abstract base classes of spaces that live on the edges or faces of a mesh, e.g. the Lagrange multiplier spaces `DGDkEdge` and `DGDkFace` of hybridized methods; they do not derive from `FE`. `FEVector`, `FESymMatrix`, `FESkewSymMatrix` (2d) and `FEEdgeVector` build vector- and matrix-valued spaces from scalar component spaces, `FEMatrix` from one scalar element per entry or from two vector-valued row elements such as `BDMk` (2d). The subfolders contain the concrete element families.

## Key concepts
- Element object: basis functions on each mesh element, the local-to-global map of the degrees of freedom (`dof4e`) and their number (`nrDoF`).
- Meshes: the 2d elements live on the triangular `Mesh`; `Pk`, `HPk`, `RTk`, `BDMk` (with their subclasses such as `DGPk`, `HRTk`, `HBDMk`) and `DGDk` also work on a `MeshSurface`. The elements ending in `3d` and the face spaces `DGDkFace` and `DGPkFace` live on the tetrahedral `Mesh3d`, and `PDGk` (and `PDGkVector`, `PDGkMatrix`, `PDGkSymMatrix` in 2d) on the polygonal `MeshPoly`.
- Conformity: H1-conforming spaces are globally continuous (`Pk`, `Pk3d`, `HPk`, `HPk3d`). H(div)-conforming spaces have continuous normal components (`RTk`, `BDMk`, `RTk3d`, `BDMk3d`), H(curl)-conforming spaces continuous tangential components (`Nd1k`, `Nd2k`, `Nd1k3d`, `Nd2k3d`). H2-conforming are `Argyris` and `TUBAk` (every degree >= 5). Symmetric matrix-valued stresses: H(div)-conforming for elasticity (`AW0`, `AW1`, `AWk`, `AWH`, `HZ1`) with the piecewise rigid body motions `RBM` as displacements, and with continuous normal-normal component (`HHJk`, `HHJB`).
- Nonconforming spaces (`CRk`, `AWNC0`, `AWNC1`, `GGNC1`) and discontinuous spaces (`DGPk`, `DGDk`, `DGDk3d`, `PDGk`, `PDGk3d`, `Nd1DGk3d`) for nonconforming and discontinuous Galerkin methods; the discontinuous spaces number their degrees of freedom per element only (`dof4e`), also those derived from a continuous sibling (`DGPk`, `DGPk3d`, `DGHPk`, `Nd1DGk3d`). The hybridized `HRTk` and `HBDMk` also have element-local degrees of freedom, coupled through Lagrange multipliers, and are not H(div)-conforming themselves. The scaled monomials of `PDGk` and `PDGk3d` live on the physical elements (`isRefElemBased` false), so these spaces have no `prolongation`.
- Construction: every element is created by `Name(degree,mesh)` with its polynomial degree; only the composite `FEVector`, `FEMatrix`, `FESymMatrix`, `FESkewSymMatrix` and `FEEdgeVector` take component spaces. Elements of one fixed degree check it (`Argyris(5,mesh)`, `AW1(3,mesh)`, `AWNC1(2,mesh)`, `RBM(1,mesh)`; `CRk` takes 1 or 3), and variants are subclasses: `AWH` hybridizes the vertex constraints of `AWk`, `HHJB` adds the TDNNS bubbles to `HHJk`, `HRTk` and `HBDMk` are the hybridized `RTk` and `BDMk`. So for a uniform degree `feval(class(fe),fe.degree,mesh)` builds the same space on another mesh, which the prolongation uses.
- Vector- and matrix-valued spaces of one scalar family (`PkVector`, `DGPkMatrix`, `DGPkSymMatrix`, `PDGkVector`) take `(degree,mesh)` and stack copies of the scalar element.
- Bubble spaces (`Bk`, `B2k`, `CurlBk`, `CurlBRotk`) enrich other spaces through `plus`, e.g. the MINI element `PkVector(1,mesh) + FEVector(Bk(3,mesh),Bk(3,mesh))` or `FEMatrix(ROW1,ROW2) + CurlBRotk(degree,mesh)`; the prolongation of such a sum is a local L2 projection.
- hp spaces with a polynomial degree per element (`HPk`, `DGHPk`); `DGDk` and the bubble spaces built on it also accept one degree per element, `HPk3d` takes one degree for all elements.
- Interpolation and projection: `interpolate` and `project` take a function of the physical points, e.g. `@p.u4Db`, or `@(pts,pts_ref,parts)` e.g. to interpolate or project another finite element function by its `getValue`, the same conventions as the functions in the forms (`FormEval.evalFunction`). `interpolate` evaluates it at the reference dof points `c4refDof` of the Lagrange elements (`Pk`, `DGPk`, `Pk3d`, `DGPk3d` and the vector- and matrix-valued spaces built from them). `project` solves the local L2 projection on every element and averages shared degrees of freedom, which is the L2 projection in the discontinuous spaces.

## Relation to other parts
Element objects are created on a mesh (from `mesh`), e.g. `Pk(1,mesh)`, and passed to the `assemble` methods of the `forms` classes. The problem classes (from `problems`) provide the data; which element is used is chosen in the calling script. The multilevel solvers of `solve` (`Multigrid`, `Schwarz`) take the prolongation matrices of the elements (`prolongation`) between refinement levels.

## Contents

{{SUBFOLDERS}}

Classes directly in the elements folder:

{{MEMBERS}}

## Related examples

{{RELATED_EXAMPLES}}
