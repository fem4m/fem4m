classdef problemSurfaceEllipsoid < ProblemSurfaceTemplate
    %problemSurfaceEllipsoid - Surface of an ellipsoid
    %   OBJ = problemSurfaceEllipsoid describes the ellipsoid
    %   x^2 + y^2 + z^2/S = R with R = 1 and S = 4 (semi-axes 1, 1 and 2)
    %   via a level-set function and its derivatives.
    %
    %   OBJ = problemSurfaceEllipsoid(R) sets the level value R; the
    %   semi-axes are sqrt(R), sqrt(R) and 2*sqrt(R).
    %
    %   OBJ = problemSurfaceEllipsoid(R,S) also sets the scaling S of z^2;
    %   the semi-axes are sqrt(R), sqrt(R) and sqrt(R*S). An empty R or S
    %   selects its default.
    %
    %   problemSurfaceEllipsoid functions:
    %       surfaceFunction     - Level-set function
    %       gradSurfaceFunction - Gradient of the level-set function
    %       d2SurfaceFunction   - Hessian of the level-set function
    %
    %   problemSurfaceEllipsoid properties:
    %       geom  - Ellipsoid geometry
    %       r     - Level value R of x^2 + y^2 + z^2/S = R (default 1)
    %       s     - Scaling S of z^2 in the level-set function (default 4)
    %
    %   See also ProblemSurfaceTemplate, problemSurfaceEllipsoidExact,
    %   problemSurfaceSphere

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Ellipsoid'))
        r
        s
    end


    methods

        function obj = problemSurfaceEllipsoid(r,s)
            if nargin < 1 || isempty(r)
                r = 1;
            end
            if nargin < 2 || isempty(s)
                s = 4;
            end
            obj.r = r;
            obj.s = s;
            obj.geom.c4n(:,3) = [sqrt(obj.s)./2.*obj.geom.c4n(1:8,3);sqrt(obj.s).*obj.geom.c4n(9:10,3)];
        end

        function val = surfaceFunction(obj,x)
            val = x(:,1).^2+x(:,2).^2+x(:,3).^2/obj.s-obj.r;
        end

        function val = gradSurfaceFunction(obj,x)
            val = [2*x(:,1), 2*x(:,2), 2*x(:,3)/obj.s];
        end

        function val = d2SurfaceFunction(obj,x)
            val = zeros(size(x,1),3,3);
            val(:,1,1) = 2;
            val(:,2,2) = 2;
            val(:,3,3) = 2/obj.s;
        end

    end

end