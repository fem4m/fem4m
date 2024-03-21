classdef problemScalarEllipticSquare < ProblemScalarEllipticTemplate
    %problemScalarEllipticSquare - Scalar elliptic problem on a square
    %   problemScalarEllipticSquare is a scalar elliptic problem on the
    %   unit square with the constant source f = 1 and homogeneous
    %   Dirichlet data on the whole boundary; the zero Neumann data g are
    %   not used with this geometry.
    %
    %   problemScalarEllipticSquare functions:
    %       f     - Source term
    %       u4Db  - Dirichlet data
    %       g     - Neumann data
    %
    %   problemScalarEllipticSquare properties:
    %       geom  - Unit square geometry
    %
    %   See also ProblemScalarEllipticTemplate,
    %   problemScalarEllipticSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry','D'))
    end

    methods

        function val = f(~,x)
            val = ones(size(x,1),1);
        end

        function val = u4Db(~,x)
            val = zeros(size(x,1),1);
        end

        function val = g(~,x,~)
            val = zeros(size(x,1),1);
        end

    end

end