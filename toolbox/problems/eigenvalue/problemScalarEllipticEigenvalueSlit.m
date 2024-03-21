classdef problemScalarEllipticEigenvalueSlit < ProblemEigenvalue
    %problemScalarEllipticEigenvalueSlit - Slit domain Dirichlet eigenvalues
    %   OBJ = problemScalarEllipticEigenvalueSlit() supplies reference
    %   Dirichlet eigenvalues of the Laplacian on the slit domain
    %   (-1,1)^2 \ [0,1]x{0}: the mode 1; other entries are NaN.
    %
    %   problemScalarEllipticEigenvalueSlit functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemScalarEllipticEigenvalueSlit properties:
    %       geom        - Slit domain geometry
    %
    %   See also ProblemEigenvalue, problemScalarEllipticEigenvalueLshape

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Slit'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 8.3713297112;
            val = eigenvalue(m);
        end

    end

end