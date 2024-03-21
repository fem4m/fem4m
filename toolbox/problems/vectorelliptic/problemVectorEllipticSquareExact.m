classdef problemVectorEllipticSquareExact < ProblemScalarEllipticTemplateExact
    %problemVectorEllipticSquareExact - Vector Laplacian on a square
    %   problemVectorEllipticSquareExact is the vector Laplace problem
    %   -Laplace(u) = f on the unit square with the polynomial exact
    %   solution u = (phi,phi), phi = x(1-x)y(1-y), and homogeneous
    %   Dirichlet data on the whole boundary. The inherited Neumann data g
    %   assume a scalar solution and do not apply; exactSigma returns
    %   grad(u) (N-by-2-by-2).
    %
    %   problemVectorEllipticSquareExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemVectorEllipticSquareExact properties:
    %       geom        - Unit square geometry
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
    end


    methods

        function val = f(~,x)
            val = 2*x(:,1).*(1-x(:,1)) + 2*x(:,2).*(1-x(:,2));
            val = [val,val];
        end

        function val = exactU(~,x)
            val = x(:,1).*(1-x(:,1)).*x(:,2).*(1-x(:,2));
            val = [val,val];
        end

        function val = exactGradU(~,x)
            grad = [ (1 - 2*x(:,1)).*x(:,2).*(1-x(:,2)),...
                x(:,1).*(1-x(:,1)).*(1-2*x(:,2))];
            val(:,1,:) = grad;
            val(:,2,:) = grad;
        end

    end

end