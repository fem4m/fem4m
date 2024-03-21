classdef problemSurfaceEigenvalueEllipsoid < problemSurfaceEllipsoid & ProblemEigenvalue
    %problemSurfaceEigenvalueEllipsoid - Ellipsoid surface eigenvalues
    %   OBJ = problemSurfaceEigenvalueEllipsoid() supplies reference
    %   Laplace-Beltrami eigenvalues of the ellipsoid surface
    %   x^2+y^2+z^2/4 = 1 with the semi-axes 1, 1 and 2: the mode 1 (zero)
    %   and the mode 2; other entries are NaN.
    %
    %   The reference eigenvalues belong to the default ellipsoid of
    %   problemSurfaceEllipsoid (R = 1, S = 4) only.
    %
    %   problemSurfaceEigenvalueEllipsoid functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   See also problemSurfaceEllipsoid, ProblemEigenvalue

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max([m(:);2]),1);
            eigenvalue(1) = 0;
            eigenvalue(2) = 0.7287949;
            val = eigenvalue(m);
        end

    end

end