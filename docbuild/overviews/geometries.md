# geometries

## Purpose
This folder contains the `Geometry` class and, in one subfolder per domain, the coarse mesh data of the domains. `Geometry()` creates the coarse triangulation of the unit square with two triangles and Dirichlet conditions on all edges. `Geometry(OPTS)` selects the domain by `OPTS.name`: `'Square'`, `'Lshape'`, `'Slit'`, `'Circle'` and `'Cube'` (case-insensitive) are built in the code, every other name loads the data files of that domain. For the built-in domains `OPTS.type` selects triangles, quadrilaterals or H-type polygons (2d domains only, `'Cube'` is tetrahedral) and `OPTS.reftype` the pattern of the coarse triangular or tetrahedral mesh; for every domain `OPTS.length`, `OPTS.rotation` and `OPTS.midpoint` scale, rotate and shift the coarse mesh; the `Geometry` help lists all options and their defaults.

## Key concepts
- `OPTS.name` selects the domain, e.g. `'Square'` (up to 4 different boundary conditions, set with `OPTS.bdry`).
- Data files per domain: `*_c4n.dat` coordinates of the nodes, `*_n4e.dat` node ids of each element, `*_n4Db.dat`/`*_n4Nb.dat`/`*_n4Sb.dat` nodes of the Dirichlet, Neumann and symmetry boundary sides (edges in 2d, faces in 3d). The last column of `*_n4Sb.dat` is the fixed displacement component (1, 2, in 3d also 3, or 0 for a general normal); `Geometry` and the meshes reject files without it. In 3d `*_t4e.dat` holds the bisection labels of the tetrahedra: with the vertex order of `*_n4e.dat` (refinement edge between the first and last vertex) they fix the newest vertex bisection of `Mesh3d`; without the file `Mesh3d` takes the labels from the longest edges.

## Relation to other parts
Geometry objects initialize `Mesh`, `Mesh3d`, `MeshSurface` and the `MeshPoly` classes (from `mesh`); `MeshPolyQ` and `MeshPolyH` rebuild a geometry of another type with `OPTS.type` `'Q'` or `'H'`. Problem classes (from `problems`) store their domain as a `Geometry` in the property `geom`.

## Contents

{{MEMBERS}}

### Geometry data folders

The built-in `'Circle'` shadows the files of the folder circle, and `'Cube'` with `OPTS.type` `'T'` those of the folder cube; these files are read through the file name prefix, e.g. `Mesh3d('Cube')`.

| Folder | Data files |
|---|---|
| backwardfacingstep | BackwardFacingStep_c4n, _n4Db, _n4Nb, _n4e (.dat) |
| circle | Circle_c4n, _n4Db, _n4e (.dat) |
| circularinclusion | CircularInclusion_c4n, _n4Db, _n4Nb, _n4Sb, _n4e (.dat) |
| cooks | Cooks_c4n, _n4Db, _n4Nb, _n4e (.dat) |
| cube | Cube_c4n, _n4Db, _n4e, _t4e (.dat) |
| dblconnected | DblConnected_c4n, _n4Db, _n4e (.dat) |
| dumbbell | Dumbbell_c4n, _n4Db, _n4e (.dat) |
| fichera | Fichera_c4n, _n4Db, _n4e, _t4e (.dat) |
| hexagon | Hexagon_c4n, _n4Db, _n4e (.dat) |
| isospectraldrums | IsospectralDrumsA_c4n, _n4Db, _n4e; IsospectralDrumsB_c4n, _n4Db, _n4e (.dat) |
| lshape3d | Lshape3d_c4n, _n4Db, _n4e, _t4e (.dat) |
| pentagon | Pentagon_c4n, _n4Db, _n4e (.dat) |
| sevenquarter | SevenQuarter_c4n, _n4Db, _n4e (.dat) |
| surface | subfolders ellipsoid, g88, sphere, spherel, squaresurface, torus (each _c4n, _n4e; spherel and squaresurface also _n4Db; the torus files are named TorusFine_*) |
| triangle | Triangle_c4n, _n4Db, _n4e (.dat) |

## Related examples

{{RELATED_EXAMPLES}}
