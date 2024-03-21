classdef problemSurfaceStokesSquareExact < ProblemStokesTemplateExact & ProblemSurfaceTemplate
    %problemSurfaceStokesSquareExact - Stokes flow on a flat square
    %   problemSurfaceStokesSquareExact is a surface Stokes problem on the
    %   flat square (-1,1)^2 in the plane z = 0 with the exact velocity
    %   u = (20xy^3, 5x^4 - 5y^4, 0), the exact pressure p = 60x^2y - 20y^3
    %   and zero body force for the viscosity NU = 1. The whole boundary is
    %   a Dirichlet boundary with the data u4Db = exactU. The flat analogue
    %   is problemStokesExactSquareP3Sigma.
    %
    %   problemSurfaceStokesSquareExact functions:
    %       f                   - Zero body force
    %       exactU              - Exact velocity
    %       exactGradU          - Velocity gradient
    %       exactP              - Exact pressure
    %       surfaceFunction     - Level-set function z of the flat square
    %       gradSurfaceFunction - Gradient of the level-set function
    %       d2SurfaceFunction   - Hessian of the level-set function
    %       curvature           - Principal curvatures (both zero)
    %
    %   problemSurfaceStokesSquareExact properties:
    %       geom                - Square surface geometry
    %
    %   See also ProblemStokesTemplateExact, problemSurfaceSquare,
    %   problemStokesExactSquareP3Sigma, problemSurfaceStokesSphereExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','SquareSurface','midpoint',[-1,-1],'length',2))
    end

    methods
        function val = f(~,x)
            val = zeros(size(x,1),3);
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = zeros(size(x,1),3);
            val(:,1) = 20*x.*y.^3;
            val(:,2) = 5*x.^4 - 5*y.^4;
        end

        function val = exactGradU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = zeros(size(x,1),3,3);
            val(:,1,1) = 20*y.^3;
            val(:,1,2) = 20*x.*3.*y.^2;
            val(:,2,1) = 20*x.^3;
            val(:,2,2) = - 20*y.^3;
        end

        function val = exactP(~,x)
            y = x(:,2);
            x = x(:,1);
            val = (60*x.^2.*y - 20*y.^3);
        end

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