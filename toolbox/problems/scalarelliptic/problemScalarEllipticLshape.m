classdef problemScalarEllipticLshape < ProblemScalarEllipticTemplate
    %problemScalarEllipticLshape - Scalar elliptic problem on an L-shape
    %   problemScalarEllipticLshape is a scalar elliptic problem on the
    %   L-shape (-1,1)^2 \ [0,1]x[-1,0] with the constant source f = 1 and
    %   homogeneous Dirichlet data on the whole boundary; the zero Neumann
    %   data g are not used with this geometry.
    %
    %   problemScalarEllipticLshape functions:
    %       f     - Source term
    %       u4Db  - Dirichlet data
    %       g     - Neumann data
    %
    %   problemScalarEllipticLshape properties:
    %       geom  - L-shaped domain geometry
    %
    %   See also ProblemScalarEllipticTemplate,
    %   problemScalarEllipticLshapeExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape','bdry','D'))
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