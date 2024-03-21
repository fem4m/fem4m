classdef problemMaxwellCubeExact
    %problemMaxwellCubeExact - Maxwell on a cube, polynomial exact field
    %   problemMaxwellCubeExact is a Maxwell source problem
    %   curl(curl(u)) = f on the unit cube with the manufactured polynomial
    %   exact solution u = (y(1-y)z(1-z), x(1-x)z(1-z), x(1-x)y(1-y)),
    %   whose tangential trace vanishes on the boundary.
    %
    %   problemMaxwellCubeExact functions:
    %       f           - Source current
    %       curlf       - Curl of the source current
    %       u4Db        - Dirichlet field (exact solution)
    %       g           - Zero Neumann-type boundary data
    %       exactU      - Exact solution
    %       exactCurlU  - Curl of the exact solution
    %       mu          - Permeability (constant 1)
    %
    %   problemMaxwellCubeExact properties:
    %       geom        - Cube geometry
    %
    %   See also problemMaxwellCube, problemMaxwellCubeExact2,
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
            val =  [- 2*y.*(y - 1) - 2*z.*(z - 1),...
                - 2*x.*(x - 1) - 2*z.*(z - 1),...
                - 2*x.*(x - 1) - 2*y.*(y - 1)];
        end

        function val = curlf(~,x)
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
            val = 4*[z - y, x - z, y - x];
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
            val = [y.*(1-y).*z.*(1-z),x.*(1-x).*z.*(1-z),x.*(1-x).*y.*(1-y)];
        end

        function val = exactCurlU(~,x)
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
            val= [   2*x.*(y - z).*(x - 1),...
                -2*y.*(x - z).*(y - 1),...
                2*z.*(x - y).*(z - 1)];
        end

        function val = mu(~,x)
            val = ones(size(x,1),1);
        end

    end
end