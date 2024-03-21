% INTERNAL  Numerical, plotting and repository helper classes and functions
%
%   Classes
%   BarycentricJacobi              -  Derivatives of barycentric Jacobi factors
%   fem4mPath                      -  Locate and prepare the fem4m repository
%   fem4mResult                    -  File name and contents of level-wise results
%   JacobiPolynomials              -  Jacobi polynomials and their integrals
%   Monomials                      -  Monomials and the matrices of their derivatives
%   SurfaceTransformation          -  Lifted surface data and geometric error coefficients
%
%   Functions
%   AitkenExtrapolation            -  Aitken extrapolation of a sequence
%   codeIndentation                -  Re-indent all M-files below a folder
%   drawConvergenceSlope           -  Draw reference slope in convergence plot
%   drawTriangleSlope              -  Draw slope triangle in loglog plot
%   fem4mLogo                      -  Render and save the fem4m logo
%   getBatchMatrix                 -  All combinations of the entries of several vectors
%   getBatchMatrixIndicesToCompute -  Experiments without a results file
%   getEdgeValues4ElemFunktion     -  Element function values on both elements of edges
%   getEdgeValues4RefElemFunktion  -  Reference element function values on all edges
%   getFaceValues4ElemFunktion     -  Element function values on both elements of faces
%   getFaceValues4RefElemFunktion  -  Reference element function values on all faces
%   integrate                      -  Vectorized Gauss integration over mesh parts
%   intUniqueRows                  -  Stable unique rows as integer indices
%   localDoFtoGlobalDoF            -  Global index pairs of local element dofs
%   mat2data                       -  Write .dat files of the results referenced in a LaTeX file
%   matMul                         -  Elementwise matrix product along first dimension
%   plotConvergence                -  Plot error curve in loglog scale
%   plotFunc4RefinedMesh           -  Plot element function on subtriangulation
%   plotRefFunc4RefinedMesh        -  Plot reference function on subtriangulation
%   projectVertices2Surface        -  Project vertices onto implicit surface
%   prolongationMatrix             -  Prolongation between finite element spaces on refined meshes
%   rearrangeLegend                -  Reorder or fade legend of current axes
%   resetFigurePosition            -  Restore default figure position
%   rowaddr                        -  Entries of a matrix at row and column index pairs
%   saveFigureHighDpiPng           -  Save figure as high-resolution image
%   spyEps                         -  Display sparsity pattern ignoring tiny entries
%   testCoverage                   -  Run the tests with a code coverage report