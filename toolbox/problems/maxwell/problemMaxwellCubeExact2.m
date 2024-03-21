classdef problemMaxwellCubeExact2
    %problemMaxwellCubeExact2 - Maxwell on a cube, trigonometric exact field
    %   problemMaxwellCubeExact2 is a Maxwell source problem
    %   curl(curl(u)) = f on the unit cube with the manufactured
    %   trigonometric exact solution
    %   u = (sin(pi*y)sin(pi*z), sin(pi*x)sin(pi*z), sin(pi*x)sin(pi*y)),
    %   whose tangential trace vanishes on the boundary.
    %
    %   problemMaxwellCubeExact2 functions:
    %       f           - Source current
    %       curlf       - Curl of the source current
    %       u4Db        - Dirichlet field (exact solution)
    %       g           - Zero Neumann-type boundary data
    %       exactU      - Exact solution
    %       exactCurlU  - Curl of the exact solution
    %       mu          - Permeability (constant 1)
    %
    %   problemMaxwellCubeExact2 properties:
    %       geom        - Cube geometry
    %
    %   See also problemMaxwellCube, problemMaxwellCubeExact,
    %   problemMaxwellLshape3dExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Cube','midpoint',[0.5,0.5,0.5],'length',1))
    end


    methods

        function val = f(~,x)
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
            val =  [ 2*pi^2*sin(pi*y).*sin(pi*z),...
                2*pi^2*sin(pi*x).*sin(pi*z),...
                2*pi^2*sin(pi*x).*sin(pi*y)];
        end

        function val = curlf(obj,x)
            val = 2*pi^2*obj.exactCurlU(x);
        end

        function val = u4Db(obj,x)
            val = obj.exactU(x);
        end

        function val = g(~,x)
            val = zeros(size(x,1),3);
        end

        function val = exactU(~,x)
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
            val = [sin(pi*y).*sin(pi*z),sin(pi*x).*sin(pi*z),sin(pi*x).*sin(pi*y)];
        end

        function val = exactCurlU(~,x)
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
            val= [pi*sin(pi*x).*(cos(pi*y) - cos(pi*z)),...
                -pi*sin(pi*y).*(cos(pi*x) - cos(pi*z)),...
                pi*sin(pi*z).*(cos(pi*x) - cos(pi*y))];
        end

        function val = mu(~,x)
            val = ones(size(x,1),1);
        end

    end
end