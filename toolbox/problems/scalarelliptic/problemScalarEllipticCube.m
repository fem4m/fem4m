classdef problemScalarEllipticCube < ProblemScalarEllipticTemplate
    %problemScalarEllipticCube - Scalar elliptic problem on a cube
    %   problemScalarEllipticCube is a scalar elliptic problem on the unit
    %   cube with the constant source f = 1 and homogeneous Dirichlet data
    %   on the whole boundary; the zero Neumann data g are not used with
    %   this geometry.
    %
    %   problemScalarEllipticCube functions:
    %       f     - Source term
    %       u4Db  - Dirichlet data
    %       g     - Neumann data
    %
    %   problemScalarEllipticCube properties:
    %       geom  - Unit cube geometry
    %
    %   See also ProblemScalarEllipticTemplate,
    %   problemScalarEllipticCubeExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Cube','midpoint',[0.5,0.5,0.5],'length',1))
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