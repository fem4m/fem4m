classdef problemScalarEllipticEigenvalueSquareNeumann < ProblemEigenvalue
    %problemScalarEllipticEigenvalueSquareNeumann - Square Neumann eigenvalues
    %   OBJ = problemScalarEllipticEigenvalueSquareNeumann() supplies the
    %   analytic Neumann eigenvalues pi^2*(j^2+k^2), j,k >= 0, of the
    %   Laplacian on the unit square (0,1)^2 for every mode; mode 1 is the
    %   zero eigenvalue of the constants.
    %
    %   problemScalarEllipticEigenvalueSquareNeumann functions:
    %       eigenvalue  - Analytic eigenvalue of mode M
    %
    %   problemScalarEllipticEigenvalueSquareNeumann properties:
    %       geom        - Unit square geometry
    %
    %   See also ProblemEigenvalue, problemScalarEllipticEigenvalueSquare,
    %   problemScalarEllipticEigenvalueLshapeNeumann

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry','N'))
    end

    methods

        function val = eigenvalue(~,m)
            j = 0:max(m);
            k = (0:max(m))';
            eigenvalue = pi^2*(j.^2+k.^2);
            eigenvalue = sort(eigenvalue(:));
            val = eigenvalue(m);
        end

    end

end