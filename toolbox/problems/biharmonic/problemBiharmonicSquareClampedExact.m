classdef problemBiharmonicSquareClampedExact < ProblemBiharmonicTemplateExact
    %problemBiharmonicSquareClampedExact - Clamped on unit square, homogeneous data
    %   problemBiharmonicSquareClampedExact defines a clamped biharmonic
    %   problem on the unit square with the polynomial exact solution
    %   u = x^2*(x-1)^2*y^2*(y-1)^2, whose value and normal derivative
    %   vanish on the boundary.
    %
    %   problemBiharmonicSquareClampedExact functions:
    %       f           - Source term from the exact solution
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %       exactD2     - Exact second derivatives
    %
    %   problemBiharmonicSquareClampedExact properties:
    %       geom        - Unit square geometry
    %
    %   See also ProblemBiharmonicTemplateExact,
    %   problemBiharmonicSquareClampedExact2

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
    end


    methods

        function val = f(~,x,varargin)
            val = 24.*x(:,2).^2.*(x(:,2)-1).^2 + 24.*x(:,1).^2.*(x(:,1)-1).^2 ...
                + 2*4.*(6*x(:,1).^2 - 6*x(:,1) + 1).*(6*x(:,2).^2 - 6*x(:,2) + 1);
        end

        function val = exactU(~,x)
            val = x(:,1).^2.*(x(:,1)-1).^2.*x(:,2).^2.*(x(:,2)-1).^2;
        end

        function val = exactGradU(~,x)
            val = [2*x(:,1).*x(:,2).^2.*(x(:,2) - 1).^2.*(2*x(:,1) - 1).*(x(:,1) - 1),...
                2*x(:,1).^2.*x(:,2).*(2*x(:,2) - 1).*(x(:,2) - 1).*(x(:,1) - 1).^2];
        end

        function val = exactD2(~,x)
            val = [2*x(:,2).^2.*(x(:,2) - 1).^2.*(6*x(:,1).^2 - 6*x(:,1) + 1),...
                2*x(:,1).^2.*(6*x(:,2).^2 - 6*x(:,2) + 1).*(x(:,1) - 1).^2,...
                4*x(:,1).*x(:,2).*(2*x(:,2) - 1).*(x(:,2) - 1).*(2*x(:,1) - 1).*(x(:,1) - 1)];
        end

    end

end