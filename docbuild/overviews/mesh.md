# mesh

## Purpose
Triangular, tetrahedral, polygonal and surface mesh classes with uniform and adaptive refinement, plus marking strategies that select the elements or edges to refine in an adaptive loop.

`Mesh` is the triangular 2d mesh. `refine(obj)` refines the marked edges by red-green-blue refinement (uniformly, with a warning, if nothing is marked), `refine(obj,'uniform')` refines all elements, `refine(obj,'NVB')` uses newest vertex bisection, and `refine(obj,'graded',...)`, `refine(obj,'gradedbeta',...)` and `refine(obj,'gradedred',...)` build meshes graded towards given points. All of these keep the mesh conforming (no hanging nodes); `coarse` undoes the last newest vertex bisection at marked elements by merging sibling pairs. The class also provides geometric queries such as `area4e`, `length4ed`, `point4e`, `normal4ed` and `patch4n`.

## Key concepts
- Mesh data: `c4n` coordinates of the nodes, `n4e` node ids of each element, `n4Db`/`n4Nb`/`n4Sb` Dirichlet, Neumann and symmetry boundary sides, which together cover the boundary (in 2d every edge in the node order of its element; `Mesh(c4n,n4e)` without boundary lists makes all boundary edges Dirichlet edges), plus derived edge data (`n4ed`, `ed4e`, `e4ed`, ...).
- Connectivity convention: every `x4y` has one row per `y`, and the nonzero entries of row `y` are the ids of the `x` of `y`, e.g. `n4e(e,:)` the nodes and `ed4e(e,:)` the edges of the element `e`, `e4ed(ed,:)` the elements at the edge `ed` and `e4n(n,:)` the elements at the node `n` (`n4e` and `ed4e` are cell arrays in `MeshPoly`). Code for all meshes reads rows only; loops read the columns of the transpose once, e.g. `e4nT = mesh.e4n'`. Relations with a bounded number of entries are dense lists padded with 0 (`n4e`, `ed4e`, `n4ed`, `f4e`, `n4f`, `ed4f`, `e4f` and the 2d `e4ed`), the others sparse incidence matrices of size nrY x nrX with the id as value (`e4n` in all meshes, `f4n`, the 3d `e4ed` and `f4ed`, computed on access in `Mesh3d`). The one exception is `ed4n`, the symmetric sparse nrNodes x nrNodes lookup `ed4n(i,j) = ed4n(j,i)` of the edge between the nodes `i` and `j`; its rows also list the edges at a node.
- `Mesh3d`: tetrahedral meshes with recursive newest vertex bisection after Maubach and Traxler: `refine(marked4e)` bisects every marked element once, `refine(marked4ed)` halves the marked edges, `refine('uniform',n)` bisects every element three times per round (8 children on the initial or a uniformly refined mesh; after adaptive steps the completion can add bisections), and each call completes the mesh to a conforming one. The initial bisection labels are given (`t4e`, e.g. the `*_t4e.dat` files of the 3d geometries) or come from the longest edges (Arnold, Mukherjee, Pouly: every element is bisected first at its longest edge, elements of the types M and O once by their rules), or with `'init','coloring'` from a greedy vertex coloring (Diening, Gehring, Storn); all elements stay positively oriented. Every boundary face is in exactly one of `n4Db`, `n4Nb` and `n4Sb`, whose fourth column is the fixed displacement component of a gliding (symmetry) boundary.
- `MeshPoly`, `MeshPolyT`, `MeshPolyQ`, `MeshPolyH`: polygon meshes for polygonal methods. `MeshPoly` itself cannot be refined; `MeshPolyT` and `MeshPolyQ` refine marked triangles and quadrilaterals with hanging nodes, stored as extra polygon vertices, and `MeshPolyH` (Geometry type 'H') refines without hanging nodes.
- `MeshSurface`: triangular meshes of surfaces embedded in 3d; the vertices are projected onto the zero level set of a given function, also after every refinement, optionally with isoparametric elements of higher degree.
- Marking strategies: `MarkUniform`, `MarkBulk` (Dörfler marking with parameter theta), `MarkMaximum`, `MarkMean`, `MarkFixedFraction` (a fixed number of the largest indicators, scaled by the polynomial degree), `MarkP` (pure p-refinement) and the hp strategies `MarkHP` and `MarkHP3Edges` (for elements with an elementwise degree `degree4e`, e.g. `HPk`), all derived from the abstract `Mark`, whose constructor sets their parameters. `mark` stores the marks in `marked4ed`, `marked4e` and `markedP4e`, read by `mesh = mesh.refine(mark.marked4ed)` in 2d or `mesh = mesh.refine(mark.marked4e)` for a `Mesh3d`.
- Refinement data: every refine call of `Mesh`, `MeshSurface` and `Mesh3d` stores `parent4e`, the ancestor of every element in the mesh before the call (also over several rounds), and `prolongation`, the P1 prolongation matrix; `coarse` sets both to [] and stores `children4e` instead. The polygon meshes keep no refinement data. `parent4e` and `prolongation` hold no finite element information; `FE.prolongation` computes the prolongation of any element from them, e.g. for multigrid.

## Relation to other parts
Element spaces (from `elements`) are created on a mesh; marking strategies feed `refine` in the adaptive examples; `Multigrid` and the `Schwarz` preconditioners (from `solve`) use the hierarchy of refined meshes. Meshes are usually created from a `Geometry` (from `geometries`).

## Contents

{{MEMBERS}}

## Related examples

{{RELATED_EXAMPLES}}
