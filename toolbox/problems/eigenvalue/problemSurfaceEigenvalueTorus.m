classdef problemSurfaceEigenvalueTorus < problemSurfaceTorus & ProblemEigenvalue
    %problemSurfaceEigenvalueTorus - Torus surface eigenvalues
    %   OBJ = problemSurfaceEigenvalueTorus() supplies reference
    %   Laplace-Beltrami eigenvalues of the torus surface with the default
    %   radii r0 = 1 and r1 = 1/4: the modes 1 (zero) to 9; other entries
    %   are NaN.
    %
    %   problemSurfaceEigenvalueTorus functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   See also problemSurfaceTorus, ProblemEigenvalue

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max([m(:);9]),1);
            eigenvalue(1) = 0;
            eigenvalue(2:3) = 1.02421785477;
            eigenvalue(4:5) = 3.9994799947;
            eigenvalue(6:7) = 8.6784015631;
            eigenvalue(8:9) = 14.790882933;
            val = eigenvalue(m);
        end

    end

end