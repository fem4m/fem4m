classdef problemSurfaceEigenvalueSquare < problemSurfaceSquare & ProblemEigenvalue
    %problemSurfaceEigenvalueSquare - Flat square surface eigenvalues
    %   OBJ = problemSurfaceEigenvalueSquare() supplies the analytic
    %   Laplace-Beltrami eigenvalues pi^2*(j^2+k^2), j,k >= 1, of the flat
    %   unit square surface [0,1]^2 x {0} with Dirichlet boundary for every
    %   mode.
    %
    %   problemSurfaceEigenvalueSquare functions:
    %       eigenvalue  - Analytic eigenvalue of mode M
    %
    %   See also problemSurfaceSquare, ProblemEigenvalue,
    %   problemScalarEllipticEigenvalueSquare

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = eigenvalue(~,m)
            j = 1:max(m);
            k = (1:max(m))';
            eigenvalue = pi^2*(j.^2+k.^2);
            eigenvalue = sort(eigenvalue(:));
            val = eigenvalue(m);
        end

    end

end