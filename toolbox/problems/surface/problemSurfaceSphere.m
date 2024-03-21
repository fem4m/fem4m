classdef problemSurfaceSphere < ProblemSurfaceTemplate
    %problemSurfaceSphere - Surface of a sphere
    %   problemSurfaceSphere describes the sphere |x| = r (property r,
    %   default 1) via the level-set function |x| - r, its derivatives,
    %   and its curvature.
    %
    %   problemSurfaceSphere functions:
    %       surfaceFunction     - Level-set function
    %       gradSurfaceFunction - Gradient of the level-set function
    %       d2SurfaceFunction   - Hessian of the level-set function
    %       curvature           - Principal curvatures (both 1/r)
    %
    %   problemSurfaceSphere properties:
    %       geom  - Sphere geometry
    %       r     - Radius (default 1)
    %
    %   See also ProblemSurfaceTemplate, problemSurfaceSphereExact,
    %   problemSurfaceSphereL, problemSurfaceEllipsoid

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Sphere'))
        r = 1
    end


    methods

        function val = surfaceFunction(obj,x)
            val = sqrt(sum(x.^2,2))-obj.r;
        end

        function val = gradSurfaceFunction(~,x)
            val = [x(:,1)./sqrt(sum(x.^2,2)),x(:,2)./sqrt(sum(x.^2,2)),x(:,3)./sqrt(sum(x.^2,2))];
        end

        function val = d2SurfaceFunction(~,x)
            R = sqrt(sum(x.^2,2));
            val(:,1,1) = 1./R - x(:,1).^2./R.^3;
            val(:,2,2) = 1./R - x(:,2).^2./R.^3;
            val(:,3,3) = 1./R - x(:,3).^2./R.^3;
            val(:,1,2) = -x(:,1).*x(:,2)./R.^3;
            val(:,1,3) = -x(:,1).*x(:,3)./R.^3;
            val(:,2,3) = -x(:,2).*x(:,3)./R.^3;
            val(:,2,1) = val(:,1,2);
            val(:,3,1) = val(:,1,3);
            val(:,3,2) = val(:,2,3);
        end

        function val = curvature(obj,x)
            val = ones(size(x,1),2)./obj.r;
        end

    end

end