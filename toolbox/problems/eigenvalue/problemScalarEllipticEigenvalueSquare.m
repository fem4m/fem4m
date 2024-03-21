classdef problemScalarEllipticEigenvalueSquare < ProblemEigenvalue
    %problemScalarEllipticEigenvalueSquare - Unit square Dirichlet eigenvalues
    %   OBJ = problemScalarEllipticEigenvalueSquare() supplies the analytic
    %   Dirichlet eigenvalues pi^2*(j^2+k^2), j,k >= 1, of the Laplacian on
    %   the unit square (0,1)^2 for every mode.
    %
    %   problemScalarEllipticEigenvalueSquare functions:
    %       eigenvalue  - Analytic eigenvalue of mode M
    %
    %   problemScalarEllipticEigenvalueSquare properties:
    %       geom        - Unit square geometry
    %
    %   See also ProblemEigenvalue,
    %   problemScalarEllipticEigenvalueSquareNeumann,
    %   problemScalarEllipticEigenvalueSquareTwist,
    %   problemSurfaceEigenvalueSquare

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
    end

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