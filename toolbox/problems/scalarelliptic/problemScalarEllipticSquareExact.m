classdef problemScalarEllipticSquareExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticSquareExact - Elliptic problem on the unit square
    %   problemScalarEllipticSquareExact is a Poisson problem on the unit
    %   square with the exact solution u = x(1-x)y(1-y) and homogeneous
    %   Dirichlet data.
    %
    %   problemScalarEllipticSquareExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticSquareExact properties:
    %       geom        - Unit square geometry
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticSquare, problemVectorEllipticSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
    end


    methods

        function val = f(~,x)
            val = 2*x(:,1).*(1-x(:,1)) + 2*x(:,2).*(1-x(:,2));
        end

        function val = exactU(~,x)
            val = x(:,1).*(1-x(:,1)).*x(:,2).*(1-x(:,2));
        end

        function val = exactGradU(~,x)
            val = [ (1 - 2*x(:,1)).*x(:,2).*(1-x(:,2)),...
                x(:,1).*(1-x(:,1)).*(1-2*x(:,2))];
        end

    end

end