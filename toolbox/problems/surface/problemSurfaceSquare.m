classdef problemSurfaceSquare < ProblemSurfaceTemplate
    %problemSurfaceSquare - A flat square surface
    %   problemSurfaceSquare describes the unit square (0,1)^2 in the
    %   plane z = 0 in 3d via the level-set function z, its derivatives,
    %   and its curvature.
    %
    %   problemSurfaceSquare functions:
    %       surfaceFunction     - Level-set function
    %       gradSurfaceFunction - Gradient of the level-set function
    %       d2SurfaceFunction   - Hessian of the level-set function
    %       curvature           - Principal curvatures (both zero)
    %
    %   problemSurfaceSquare properties:
    %       geom  - Square surface geometry
    %
    %   See also ProblemSurfaceTemplate, problemSurfaceSquareExact,
    %   problemSurfaceSphere

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','SquareSurface'))
    end


    methods

        function val = surfaceFunction(~,x)
            val = x(:,3);
        end

        function val = gradSurfaceFunction(~,x)
            val = zeros(size(x,1),3);
            val(:,3) = 1;
        end

        function val = d2SurfaceFunction(~,x)
            val = zeros(size(x,1),3,3);
        end

        function val = curvature(~,x)
            val = zeros(size(x,1),2);
        end
    end

end