classdef problemScalarEllipticLbrick < ProblemScalarEllipticTemplate
    %problemScalarEllipticLbrick - Scalar elliptic problem on an L-brick
    %   problemScalarEllipticLbrick is a scalar elliptic problem on the
    %   L-brick ((-1,1)^2 \ [0,1]x[-1,0]) x (0,1) with the constant source
    %   f = 1 and homogeneous Dirichlet data on the whole boundary; the
    %   zero Neumann data g are not used with this geometry.
    %
    %   problemScalarEllipticLbrick functions:
    %       f     - Source term
    %       u4Db  - Dirichlet data
    %       g     - Neumann data
    %
    %   problemScalarEllipticLbrick properties:
    %       geom  - L-shaped 3d geometry
    %
    %   See also ProblemScalarEllipticTemplate, problemScalarEllipticCube,
    %   problemScalarEllipticLbrickExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape3d'))
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