classdef problemSurfaceG88 < ProblemSurfaceTemplate
    %problemSurfaceG88 - The G88 benchmark surface
    %   problemSurfaceG88 describes the closed G88 surface
    %   (x - z^2)^2 + y^2 + z^2 = 1 of Dziuk (1988), the unit sphere
    %   deformed by (x,y,z) -> (x + z^2,y,z), via an implicit level-set
    %   function and its derivatives.
    %
    %   problemSurfaceG88 functions:
    %       surfaceFunction     - Level-set function
    %       gradSurfaceFunction - Gradient of the level-set function
    %       d2SurfaceFunction   - Hessian of the level-set function
    %
    %   problemSurfaceG88 properties:
    %       geom  - G88 surface geometry
    %
    %   See also ProblemSurfaceTemplate, problemSurfaceStokesG88,
    %   problemSurfaceEllipsoid

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','G88'))
    end


    methods

        function val = surfaceFunction(~,x)
            val = (x(:,1)-x(:,3).^2).^2+x(:,2).^2+x(:,3).^2-1;
        end

        function val = gradSurfaceFunction(~,x)
            val = [ 2*(x(:,1)-x(:,3).^2), 2*x(:,2), 2*x(:,3)-4.*x(:,3).*(x(:,1)-x(:,3).^2)];
        end

        function val = d2SurfaceFunction(~,x)
            val = zeros(size(x,1),3,3);
            val(:,1,1) = 2;
            val(:,1,3) = -4*x(:,3);
            val(:,3,1) = val(:,1,3);
            val(:,2,2) = 2;
            val(:,3,3) = 12*x(:,3).^2 - 4*x(:,1) + 2;
        end

    end

end