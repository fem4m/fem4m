classdef problemScalarEllipticCubeExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticCubeExact - Elliptic problem on a cube, polynomial
    %   problemScalarEllipticCubeExact is a Poisson problem on the unit
    %   cube with the polynomial exact solution u = x(1-x)y(1-y)z(1-z) and
    %   homogeneous Dirichlet data.
    %
    %   problemScalarEllipticCubeExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticCubeExact properties:
    %       geom        - Unit cube geometry
    %
    %   See also ProblemScalarEllipticTemplateExact, problemScalarEllipticCube,
    %   problemScalarEllipticCubeExact2

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Cube','midpoint',[0.5,0.5,0.5],'length',1))
    end


    methods

        function val = f(~,x)
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
            val = 2*(y.*(1-y).*z.*(1-z)+x.*(1-x).*z.*(1-z)+x.*(1-x).*y.*(1-y));
        end

        function val = exactU(~,x)
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
            val = x.*(1-x).*y.*(1-y).*z.*(1-z);
        end

        function val = exactGradU(~,x)
            nrPts = size(x,1);
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
            val = zeros(nrPts,3);
            val(:,1) = (1-2*x).*y.*(1-y).*z.*(1-z);
            val(:,2) = x.*(1-x).*(1-2*y).*z.*(1-z);
            val(:,3) = x.*(1-x).*y.*(1-y).*(1-2*z);
        end

    end
end