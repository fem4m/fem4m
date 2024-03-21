# internal

## Purpose
Helper functions and small support classes used by the rest of the toolbox: numerical utilities (Jacobi polynomials, monomials, prolongation matrices, integration, matrix entries at index pairs), evaluation of functions on element edges and faces, plotting of convergence results, results files and batch runs of the experiments, and tools for the repository itself.

## Key concepts
- Polynomials and level transfer: `JacobiPolynomials` serves `HPk`, `AWk` and the Dubiner DG bases (`DGDk`, `DGDk3d`, `DGDkEdge`, `DGDkFace`); `BarycentricJacobi` provides the derivatives of the barycentric Jacobi factors of the Beuchler-Pillwein-Schoeberl-Zaglmayr bases of `HPk3d`, `HHJk` and the H(div) and H(curl) elements (`BDMk`, `BDMk3d`, `RTk3d`, `Nd1k`, `Nd2k`, `Nd1k3d`, `Nd2k3d`); `Monomials` (monomials and the integer matrices of their derivatives) serves `Pk`, `Pk3d`, `PDGk` and `PDGk3d`, and its multi-indices order the derivatives in `FE` and `BarycentricJacobi`; `prolongationMatrix` computes the prolongation matrices of `FE.prolongation`.
- Assembly helpers: `integrate` (Gauss quadrature on mesh parts), `localDoFtoGlobalDoF` (index pairs for sparse assembly), `matMul` (matrix products along the first dimension), `rowaddr` (matrix entries at row and column index pairs) and `intUniqueRows` (stable unique rows as integer indices).
- Function evaluation on sides: `getEdgeValues4ElemFunktion`, `getEdgeValues4RefElemFunktion`, `getFaceValues4ElemFunktion` and `getFaceValues4RefElemFunktion` return the values on both neighboring elements of each edge (2d and surface meshes; `getEdgeValues4ElemFunktion` also polygonal meshes) or face (3d meshes), zero on the boundary; the forms use them for jumps and averages, the elements and `MeshSurface` for values on edges and faces; `projectVertices2Surface` projects vertices onto an implicit surface, `SurfaceTransformation` lifts the right-hand side to a `MeshSurface` and gives the coefficient of the geometric error estimator.
- Plotting: `plotConvergence`, `drawConvergenceSlope`, `drawTriangleSlope`, `rearrangeLegend` for convergence plots; `plotFunc4RefinedMesh`, `plotRefFunc4RefinedMesh`, `spyEps`, `resetFigurePosition` and `saveFigureHighDpiPng` for other figures; `AitkenExtrapolation` to extrapolate sequences of values.
- Results and batch runs: `fem4mResult` (file name and contents of the level-wise results read by `plotApp` and `mat2data`), `getBatchMatrix` (all parameter combinations of a batch run) and `getBatchMatrixIndicesToCompute` (the experiments without a results file); `mat2data` writes .dat files of the results referenced in a LaTeX file.
- Repository tools: `fem4mPath` (root, results and tests folders and the MATLAB path), `codeIndentation` (re-indent all M-files below a folder), `testCoverage` (run the tests with a code coverage report) and `fem4mLogo` (render and save the logo).

## Relation to other parts
Called from the examples and from the code in `mesh`, `elements`, `forms`, `problems` and `solve`; most functions are not meant to be called directly by users.

## Contents

{{MEMBERS}}

## Related examples

{{RELATED_EXAMPLES}}
