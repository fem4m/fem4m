classdef problemScalarEllipticEigenvalueDblConnected < ProblemEigenvalue
    %problemScalarEllipticEigenvalueDblConnected - Doubly connected domain eigenvalues
    %   OBJ = problemScalarEllipticEigenvalueDblConnected() supplies
    %   reference Dirichlet eigenvalues of the Laplacian on the doubly
    %   connected domain (0,4)^2 \ [1,3]^2: the modes 1 to 3; other entries
    %   are NaN.
    %
    %   problemScalarEllipticEigenvalueDblConnected functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemScalarEllipticEigenvalueDblConnected properties:
    %       geom        - Doubly connected domain geometry
    %
    %   See also ProblemEigenvalue, problemScalarEllipticEigenvalueDumbbell

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','DblConnected'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 8.93469;
            eigenvalue(2) = 9.2030962;
            eigenvalue(3) = 9.2030962;
            val = eigenvalue(m);
        end

    end

end