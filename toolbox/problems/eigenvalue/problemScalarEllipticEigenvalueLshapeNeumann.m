classdef problemScalarEllipticEigenvalueLshapeNeumann < ProblemEigenvalue
    %problemScalarEllipticEigenvalueLshapeNeumann - L-shape Neumann eigenvalues
    %   OBJ = problemScalarEllipticEigenvalueLshapeNeumann() supplies
    %   reference Neumann eigenvalues of the Laplacian on the L-shaped domain
    %   (-1,1)^2 \ [0,1]x[-1,0]: the modes 1 to 6, starting with the zero
    %   eigenvalue of the constants; other entries are NaN.
    %
    %   problemScalarEllipticEigenvalueLshapeNeumann functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemScalarEllipticEigenvalueLshapeNeumann properties:
    %       geom        - L-shaped domain geometry
    %
    %   See also ProblemEigenvalue, problemScalarEllipticEigenvalueLshape,
    %   problemScalarEllipticEigenvalueSquareNeumann

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape','bdry','N'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 0;
            eigenvalue(2) = 0.147562182408E+01;
            eigenvalue(3) = 0.353403136678E+01;
            eigenvalue(4) = 0.986960440109E+01;
            eigenvalue(5) = 0.986960440109E+01;
            eigenvalue(6) = 0.113894793979E+02;
            val = eigenvalue(m);
        end

    end

end