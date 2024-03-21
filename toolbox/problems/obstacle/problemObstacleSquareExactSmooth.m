classdef problemObstacleSquareExactSmooth < ProblemScalarEllipticTemplateExact & ProblemObstacleTemplate
    %problemObstacleSquareExactSmooth - Smooth obstacle on a square
    %   problemObstacleSquareExactSmooth is an obstacle problem on the
    %   square (-1,1)^2 whose obstacle -(x1^2-1)*(x2^2-1) equals the exact
    %   solution: the contact set is the whole square and f = -div(grad(u)),
    %   so the multiplier vanishes.
    %
    %   problemObstacleSquareExactSmooth functions:
    %       lower       - Obstacle, equal to the exact solution
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemObstacleSquareExactSmooth properties:
    %       geom        - Square geometry
    %
    %   See also ProblemObstacleTemplate, problemObstacleSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','midpoint',[0,0],'length',2))
    end


    methods

        function val = lower(~,x)
            val = -(x(:,1).^2-1).*(x(:,2).^2-1);
        end

        function val = f(~,x)
            val = 2*( (x(:,1).^2-1) + (x(:,2).^2-1) );
        end

        function val = exactU(~,x)
            val = -(x(:,1).^2-1).*(x(:,2).^2-1);
        end

        function val = exactGradU(~,x)
            val = -2*[x(:,1).*(x(:,2).^2-1),(x(:,1).^2-1).*x(:,2)];
        end

    end

end