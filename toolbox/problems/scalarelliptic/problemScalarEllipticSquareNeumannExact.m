classdef problemScalarEllipticSquareNeumannExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticSquareNeumannExact - Pure Neumann problem on a square
    %   problemScalarEllipticSquareNeumannExact is a Poisson problem on
    %   the unit square with Neumann boundary conditions on all edges and
    %   the polynomial exact solution u = x^2*y^2*(x-1)^2*(y-1)^2 - 1/900
    %   (zero mean).
    %
    %   problemScalarEllipticSquareNeumannExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticSquareNeumannExact properties:
    %       geom        - Unit square geometry
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticSquare

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry','N'))
    end


    methods

        function val = f(~,x)
            y = x(:,2);
            x = x(:,1);
            val = x.^2.*y.^2.*(x-1.0).^2.*-2.0-x.^2.*y.^2.*(y-1.0).^2.*2.0-x.^2.*(x-1.0).^2.*(y-1.0).^2.*2.0-y.^2.*(x-1.0).^2.*(y-1.0).^2.*2.0-x.*y.^2.*(x.*2.0-2.0).*(y-1.0).^2.*4.0-x.^2.*y.*(y.*2.0-2.0).*(x-1.0).^2.*4.0;
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = x.^2.*y.^2.*(x-1.0).^2.*(y-1.0).^2 - 1/900;
        end

        function val = exactGradU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [ x.*y.^2.*(x-1.0).^2.*(y-1.0).^2.*2.0+x.^2.*y.^2.*(x.*2.0-2.0).*(y-1.0).^2,...
                x.^2.*y.*(x-1.0).^2.*(y-1.0).^2.*2.0+x.^2.*y.^2.*(y.*2.0-2.0).*(x-1.0).^2];
        end

    end

end