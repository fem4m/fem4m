classdef problemSurfaceSquareExact < problemSurfaceSquare & ProblemScalarEllipticTemplateExact
    %problemSurfaceSquareExact - Laplace-Beltrami on a flat square
    %   problemSurfaceSquareExact is a scalar elliptic (Laplace-Beltrami)
    %   problem on the unit square (0,1)^2 in the plane z = 0 with the
    %   exact solution u = x*(1-x)*y*(1-y), which vanishes on the boundary.
    %
    %   problemSurfaceSquareExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   See also problemSurfaceSquare, ProblemScalarEllipticTemplateExact

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = f(~,x)
            val = 2*x(:,1).*(1-x(:,1)) + 2*x(:,2).*(1-x(:,2));
        end

        function val = exactU(~,x)
            val = x(:,1).*(1-x(:,1)).*x(:,2).*(1-x(:,2));
        end

        function val = exactGradU(~,x)
            val = [ (1 - 2*x(:,1)).*x(:,2).*(1-x(:,2)),...
                x(:,1).*(1-x(:,1)).*(1-2*x(:,2)), ...
                zeros(size(x,1),1)];
        end

    end

end