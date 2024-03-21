classdef problemBiharmonicSquareClampedExact2 < ProblemBiharmonicTemplateExact
    %problemBiharmonicSquareClampedExact2 - Clamped on unit square, inhomogeneous data
    %   problemBiharmonicSquareClampedExact2 defines a clamped biharmonic
    %   problem on the unit square with the product-of-quadratic exact
    %   solution u = x*(x-1)*y*(y-1); u vanishes on the boundary, but its
    %   normal derivative (boundary data g) does not.
    %
    %   problemBiharmonicSquareClampedExact2 functions:
    %       f           - Constant source term 8
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %       exactD2     - Exact second derivatives
    %
    %   problemBiharmonicSquareClampedExact2 properties:
    %       geom        - Unit square geometry
    %
    %   See also ProblemBiharmonicTemplateExact,
    %   problemBiharmonicSquareClampedExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
    end


    methods

        function val = f(~,x,varargin)
            val = 8*ones(size(x,1),1);
        end

        function val = exactU(~,x)
            val = x(:,1).*(x(:,1)-1).*x(:,2).*(x(:,2)-1);
        end

        function val = exactGradU(~,x)
            val = [x(:,1).*x(:,2).*(x(:,2) - 1) + x(:,2).*(x(:,1) - 1).*(x(:,2) - 1),...
                x(:,1).*x(:,2).*(x(:,1) - 1) + x(:,1).*(x(:,1) - 1).*(x(:,2) - 1)];
        end

        function val = exactD2(~,x)
            val = [2.*x(:,2).*(x(:,2) - 1),...
                2.*x(:,1).*(x(:,1) - 1),...
                (x(:,1) - 1).*(x(:,2) - 1) + x(:,1).*x(:,2) + x(:,1).*(x(:,2) - 1) + x(:,2).*(x(:,1) - 1)];
        end

    end

end