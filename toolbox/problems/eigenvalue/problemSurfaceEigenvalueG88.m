classdef problemSurfaceEigenvalueG88 < problemSurfaceG88 & ProblemEigenvalue
    %problemSurfaceEigenvalueG88 - G88 surface eigenvalues
    %   OBJ = problemSurfaceEigenvalueG88() supplies reference
    %   Laplace-Beltrami eigenvalues of the closed G88 surface
    %   (x-z^2)^2+y^2+z^2 = 1: the modes 1 (zero) to 4; other entries are
    %   NaN.
    %
    %   problemSurfaceEigenvalueG88 functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   See also problemSurfaceG88, ProblemEigenvalue

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max([m(:);4]),1);
            eigenvalue(1) = 0;
            eigenvalue(2) = 1.6746269;
            eigenvalue(3) = 1.70108230;
            eigenvalue(4) = 2.08574803;
            val = eigenvalue(m);
        end

    end

end