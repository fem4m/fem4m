classdef problemObstacleSquareCusp < ProblemScalarEllipticTemplate & ProblemObstacleTemplate
    %problemObstacleSquareCusp - Obstacle problem with a cusp
    %   problemObstacleSquareCusp is an obstacle problem on the square
    %   (-1,1)^2 whose obstacle 1-50*max(|x1|,|x2|), cut off below at -2,
    %   has a sharp peak (the cusp) of height 1 at the origin. The
    %   Dirichlet data and the source term are those of
    %   problemObstacleSquareExact; no exact solution is provided.
    %
    %   problemObstacleSquareCusp functions:
    %       lower  - Pyramid-shaped obstacle with a peak at the origin
    %       u4Db   - Dirichlet data (r^2-r0^2)^2
    %       g      - Zero Neumann data
    %       f      - Source term
    %
    %   problemObstacleSquareCusp properties:
    %       geom   - Square geometry
    %       r0     - Parameter of the Dirichlet data and f (0 < r0 < 1)
    %       d      - Dimension in the formula of f (must be 2)
    %
    %   See also ProblemObstacleTemplate, problemObstacleSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','midpoint',[0,0],'length',2))
        r0 = 0.7 % Parameter of the Dirichlet data and f (0 < r0 < 1)
        d = 2 % Dimension in the formula of f (must be 2)
    end


    methods

        function val = lower(~,x)
            val = max([-2*ones(size(x,1),1),1-50*max(abs(x),[],2)],[],2);
        end

        function val = u4Db(obj,x)
            [~,r] = cart2pol(x(:,1),x(:,2));
            val = (r.^2-obj.r0^2).^2;
        end

        function val = g(~,x,~)
            val = zeros(size(x,1),1);
        end

        function val = f(obj,x)
            [~,r] = cart2pol(x(:,1),x(:,2));
            val = -4*(2*r.^2 + obj.d*(r.^2-obj.r0^2));
            val(r<=obj.r0) = -8*obj.r0^2*(1 -( r(r<=obj.r0).^2 - obj.r0^2));
        end

    end

end