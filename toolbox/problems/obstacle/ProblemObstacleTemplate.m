classdef (Abstract) ProblemObstacleTemplate
    %ProblemObstacleTemplate - Abstract obstacle problem template
    %   ProblemObstacleTemplate provides default obstacle bounds; combine it
    %   with a differential-operator template to build an obstacle problem.
    %   The defaults are no upper bound (inf) and no lower bound (-inf);
    %   subclasses override the bound that defines their obstacle.
    %
    %   The bounds are evaluated at the points X (N-by-d, one row per
    %   point) and return an N-by-1 vector; the solution u must satisfy
    %   lower(X) <= u(X) <= upper(X).
    %
    %   ProblemObstacleTemplate functions:
    %       upper  - Upper obstacle bound (default: inf)
    %       lower  - Lower obstacle bound (default: -inf)
    %
    %   See also problemObstacleSquareExact, problemObstacleSquareCusp,
    %   problemObstacleSquareDist, problemObstacleSquareExactSmooth,
    %   problemObstacleLshapeExact, pdaSet

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = upper(~,x)
            val = inf(size(x,1),1);
        end

        function val = lower(~,x)
            val = -inf(size(x,1),1);
        end

    end
end