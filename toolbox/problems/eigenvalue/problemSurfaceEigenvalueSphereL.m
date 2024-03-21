classdef problemSurfaceEigenvalueSphereL < problemSurfaceSphereL & ProblemEigenvalue
    %problemSurfaceEigenvalueSphereL - SphereL surface eigenvalues
    %   OBJ = problemSurfaceEigenvalueSphereL() supplies the exact
    %   Laplace-Beltrami eigenvalues nu*(nu+1), nu = 2*k/3+n with k >= 1 and
    %   n >= 0, of the three-quarter unit sphere SphereL with Dirichlet
    %   boundary: the modes 1 to 11; other entries are NaN.
    %
    %   problemSurfaceEigenvalueSphereL functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   See also problemSurfaceSphereL, ProblemEigenvalue,
    %   problemSurfaceEigenvalueSphere

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max([m(:);11]),1);
            eigenvalue(1) = 10/9;
            eigenvalue(2) = 28/9;
            eigenvalue(3) = 40/9;
            eigenvalue(4) = 6;
            eigenvalue(5) = 70/9;
            eigenvalue(6:7) = 88/9;
            eigenvalue(8) = 12;
            eigenvalue(9:10) = 130/9;
            eigenvalue(11) = 154/9;
            val = eigenvalue(m);
        end

    end

end