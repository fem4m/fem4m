classdef problemBiharmonicCubeClampedExact2 < ProblemBiharmonicTemplateExact
    %problemBiharmonicCubeClampedExact2 - Clamped on unit cube, inhomogeneous data
    %   problemBiharmonicCubeClampedExact2 defines a clamped biharmonic
    %   problem on the unit cube with the product-of-quadratic exact
    %   solution u = x(1-x)*y(1-y)*z(1-z); u vanishes on the boundary, but
    %   its normal derivative (boundary data g) does not.
    %
    %   problemBiharmonicCubeClampedExact2 functions:
    %       f           - Source term from the exact solution
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %       exactD2     - Exact second derivatives
    %
    %   problemBiharmonicCubeClampedExact2 properties:
    %       geom        - Unit cube geometry
    %
    %   See also ProblemBiharmonicTemplateExact,
    %   problemBiharmonicCubeClampedExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Cube','midpoint',[0.5,0.5,0.5],'length',1))
    end


    methods

        function val = f(~,pts)
            x = pts(:,1);
            y = pts(:,2);
            z = pts(:,3);
            val = - 8.*x.*(x - 1) - 8.*y.*(y - 1) - 8.*z.*(z - 1);
        end

        function val = exactU(~,pts)
            x = pts(:,1);
            y = pts(:,2);
            z = pts(:,3);
            val = x.*(1-x).*y.*(1-y).*z.*(1-z);
        end

        function val = exactGradU(~,pts)
            nrPts = size(pts,1);
            x = pts(:,1);
            y = pts(:,2);
            z = pts(:,3);
            val = zeros(nrPts,3);
            val(:,1) = - x.*y.*z.*(y - 1).*(z - 1) - y.*z.*(x - 1).*(y - 1).*(z - 1);
            val(:,2) = - x.*y.*z.*(x - 1).*(z - 1) - x.*z.*(x - 1).*(y - 1).*(z - 1);
            val(:,3) = - x.*y.*z.*(x - 1).*(y - 1) - x.*y.*(x - 1).*(y - 1).*(z - 1);
        end

        function val = exactD2(~,pts)
            nrPts = size(pts,1);
            x = pts(:,1);
            y = pts(:,2);
            z = pts(:,3);
            val = zeros(nrPts,6);
            val(:,1) = -2.*y.*z.*(y - 1).*(z - 1);
            val(:,2) = -2.*x.*z.*(x - 1).*(z - 1);
            val(:,3) = -2.*x.*y.*(x - 1).*(y - 1);
            val(:,4) = - x.*z.*(y - 1).*(z - 1) - y.*z.*(x - 1).*(z - 1) - z.*(x - 1).*(y - 1).*(z - 1) - x.*y.*z.*(z - 1);
            val(:,5) = - x.*y.*(y - 1).*(z - 1) - y.*z.*(x - 1).*(y - 1) - y.*(x - 1).*(y - 1).*(z - 1) - x.*y.*z.*(y - 1);
            val(:,6) = - x.*y.*(x - 1).*(z - 1) - x.*z.*(x - 1).*(y - 1) - x.*(x - 1).*(y - 1).*(z - 1) - x.*y.*z.*(x - 1);
        end

    end

end