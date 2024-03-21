classdef problemScalarEllipticCubeExact2 < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticCubeExact2 - Elliptic problem on a cube, trigonometric
    %   problemScalarEllipticCubeExact2 is a Poisson problem on the unit
    %   cube with the trigonometric exact solution
    %   u = sin(pi*x)sin(pi*y)sin(pi*z) and homogeneous Dirichlet data.
    %
    %   problemScalarEllipticCubeExact2 functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticCubeExact2 properties:
    %       geom        - Unit cube geometry
    %
    %   See also ProblemScalarEllipticTemplateExact, problemScalarEllipticCube,
    %   problemScalarEllipticCubeExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Cube','midpoint',[0.5,0.5,0.5],'length',1))
    end


    methods

        function val = f(~,x)
            val = 3*pi^2*sin(pi*x(:,1)).*sin(pi*x(:,2)).*sin(pi*x(:,3));
        end

        function val = exactU(~,x)
            val = sin(pi*x(:,1)).*sin(pi*x(:,2)).*sin(pi*x(:,3));
        end

        function val = exactGradU(~,x)
            val = [pi*cos(pi*x(:,1)).*sin(pi*x(:,2)).*sin(pi*x(:,3)),...
                pi*sin(pi*x(:,1)).*cos(pi*x(:,2)).*sin(pi*x(:,3)),...
                pi*sin(pi*x(:,1)).*sin(pi*x(:,2)).*cos(pi*x(:,3))];
        end
    end
end