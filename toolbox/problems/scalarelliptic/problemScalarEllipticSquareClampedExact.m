classdef problemScalarEllipticSquareClampedExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticSquareClampedExact - Mixed boundary, clamped-type solution
    %   problemScalarEllipticSquareClampedExact is a Poisson problem on
    %   the unit square with mixed boundary conditions (Dirichlet on the
    %   bottom edge, Neumann on the others) and the polynomial exact
    %   solution u = x^2*y^2*(x-1)^2*(y-1)^2, whose value and gradient
    %   vanish on the whole boundary.
    %
    %   problemScalarEllipticSquareClampedExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticSquareClampedExact properties:
    %       geom        - Unit square geometry
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticSquare

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry',['D','N','N','N']))
    end


    methods

        function val = f(~,x)
            y = x(:,2);
            x = x(:,1);
            val = - 2*x.^2.*y.^2.*(x - 1).^2 - 2*x.^2.*y.^2.*(y - 1).^2 ...
                - 2*x.^2.*(x - 1).^2.*(y - 1).^2 - 2*y.^2.*(x - 1).^2.*(y - 1).^2 ...
                - 4*x.*y.^2.*(2*x - 2).*(y - 1).^2 - 4*x.^2.*y.*(2*y - 2).*(x - 1).^2;
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = x.^2.*y.^2.*(x - 1).^2.*(y - 1).^2;
        end

        function val = exactGradU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [ 2*x.*y.^2.*(y - 1).^2.*(2*x.^2 - 3*x + 1),...
                2*x.^2.*y.*(x - 1).^2.*(2*y.^2 - 3*y + 1)];
        end

    end

end