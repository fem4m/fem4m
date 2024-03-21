classdef problemSurfaceEigenvalueSphere < problemSurfaceSphere & ProblemEigenvalue
    %problemSurfaceEigenvalueSphere - Sphere surface eigenvalues
    %   OBJ = problemSurfaceEigenvalueSphere() supplies the analytic
    %   Laplace-Beltrami eigenvalues l*(l+1) of multiplicity 2*l+1 of the
    %   unit sphere (radius r = 1) for l = 0,...,10, i.e. the modes 1 to 121;
    %   other entries are NaN.
    %
    %   problemSurfaceEigenvalueSphere functions:
    %       eigenvalue  - Analytic eigenvalue of mode M
    %
    %   See also problemSurfaceSphere, ProblemEigenvalue,
    %   problemSurfaceEigenvalueSphereL

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max([m(:);121]),1);
            eigenvalue(1) = 0;
            eigenvalue(2:4) = 2;
            eigenvalue(5:9) = 6;
            eigenvalue(10:16) = 12;
            eigenvalue(17:25) = 20;
            eigenvalue(26:36) = 30;
            eigenvalue(37:49) = 42;
            eigenvalue(50:64) = 56;
            eigenvalue(65:81) = 72;
            eigenvalue(82:100) = 90;
            eigenvalue(101:121) = 110;
            val = eigenvalue(m);
        end

    end

end