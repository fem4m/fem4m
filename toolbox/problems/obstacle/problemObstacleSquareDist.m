classdef problemObstacleSquareDist < ProblemScalarEllipticTemplate & ProblemObstacleTemplate
    %problemObstacleSquareDist - Obstacle given by distance to boundary
    %   problemObstacleSquareDist is an obstacle problem on the unit square
    %   (0,1)^2 with the upper obstacle dist(x,boundary), zero Dirichlet
    %   data and the source term f = 7.
    %
    %   problemObstacleSquareDist functions:
    %       u4Db  - Zero Dirichlet data
    %       g     - Zero Neumann data
    %       upper - Distance-to-boundary obstacle
    %       f     - Constant source term 7
    %
    %   problemObstacleSquareDist properties:
    %       geom  - Unit square geometry
    %
    %   See also ProblemObstacleTemplate, problemObstacleSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
    end

    methods

        function val = u4Db(~,x)
            val = zeros(size(x,1),1);
        end

        function val = g(~,x,~)
            val = zeros(size(x,1),1);
        end

        function val = upper(~,x)
            val = min(abs([1-x(:,1),x(:,1),x(:,2),1-x(:,2)]),[],2);
        end

        function val = f(~,x)
            val = 7*ones(size(x,1),1);
        end

    end

end