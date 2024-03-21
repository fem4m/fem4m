classdef problemScalarEllipticEigenvalueDumbbell < ProblemEigenvalue
    %problemScalarEllipticEigenvalueDumbbell - Dumbbell domain eigenvalues
    %   OBJ = problemScalarEllipticEigenvalueDumbbell() supplies reference
    %   Dirichlet eigenvalues of the Laplacian on a dumbbell-shaped domain,
    %   two unit squares joined by a channel of width and length 1/4: the
    %   mode 1; other entries are NaN.
    %
    %   problemScalarEllipticEigenvalueDumbbell functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemScalarEllipticEigenvalueDumbbell properties:
    %       geom        - Dumbbell domain geometry
    %
    %   See also ProblemEigenvalue, problemScalarEllipticEigenvalueDblConnected

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Dumbbell'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 1.9557938*pi^2;
            val = eigenvalue(m);
        end

    end

end