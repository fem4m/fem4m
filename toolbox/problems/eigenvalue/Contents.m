% EIGENVALUE  Biharmonic, elasticity, Laplace and surface eigenvalue problems
%
%   Classes
%   problemBiharmonicEigenvalueBucklingLshapeClamped   -  Buckling on an L-shape
%   problemBiharmonicEigenvalueBucklingSlitClamped     -  Buckling on a slit
%   problemBiharmonicEigenvalueBucklingSquareClamped   -  Buckling on a square
%   problemBiharmonicEigenvalueVibrationLshapeClamped  -  Vibration on an L-shape
%   problemBiharmonicEigenvalueVibrationSlitClamped    -  Vibration on a slit
%   problemBiharmonicEigenvalueVibrationSquareClamped  -  Vibration on a square
%   ProblemEigenvalue                                  -  Abstract base class for eigenvalue problems
%   problemElasticityEigenvalueCircle                  -  Elastic disk, Dirichlet boundary
%   problemElasticityEigenvalueCircleNeumann           -  Elastic disk, Neumann boundary
%   problemElasticityEigenvalueCircleSb                -  Elastic disk, symmetry boundary
%   problemElasticityEigenvalueCooks                   -  Eigenvalues of the Cook membrane
%   problemElasticityEigenvalueLshape                  -  Elastic L-shape, Dirichlet boundary
%   problemElasticityEigenvalueLshapeNeumann           -  Elastic L-shape, Neumann
%   problemElasticityEigenvalueLshapeSb                -  Elastic L-shape, symmetry boundary
%   problemElasticityEigenvalueSquare                  -  Elastic square, Dirichlet boundary
%   problemElasticityEigenvalueSquareNeumann           -  Elastic square, Neumann
%   problemElasticityEigenvalueSquareSb                -  Elastic square, symmetry boundary
%   problemElasticityEigenvalueStripeX                 -  Elastic square, symmetry at x = 0
%   problemElasticityEigenvalueStripeY                 -  Elastic square, symmetry at y = 0
%   problemScalarEllipticEigenvalueDblConnected        -  Doubly connected domain eigenvalues
%   problemScalarEllipticEigenvalueDumbbell            -  Dumbbell domain eigenvalues
%   problemScalarEllipticEigenvalueIsospectralDrumsA   -  Isospectral drum, A
%   problemScalarEllipticEigenvalueIsospectralDrumsB   -  Isospectral drum, B
%   problemScalarEllipticEigenvalueLshape              -  L-shape Dirichlet eigenvalues
%   problemScalarEllipticEigenvalueLshapeNeumann       -  L-shape Neumann eigenvalues
%   problemScalarEllipticEigenvalueSlit                -  Slit domain Dirichlet eigenvalues
%   problemScalarEllipticEigenvalueSquare              -  Unit square Dirichlet eigenvalues
%   problemScalarEllipticEigenvalueSquareNeumann       -  Square Neumann eigenvalues
%   problemScalarEllipticEigenvalueSquareTwist         -  Twisted diffusion on the unit square
%   problemSurfaceEigenvalueEllipsoid                  -  Ellipsoid surface eigenvalues
%   problemSurfaceEigenvalueG88                        -  G88 surface eigenvalues
%   problemSurfaceEigenvalueSphere                     -  Sphere surface eigenvalues
%   problemSurfaceEigenvalueSphereL                    -  SphereL surface eigenvalues
%   problemSurfaceEigenvalueSquare                     -  Flat square surface eigenvalues
%   problemSurfaceEigenvalueTorus                      -  Torus surface eigenvalues