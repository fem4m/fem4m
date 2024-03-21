classdef problemElasticitySquareExactLocking < ProblemElasticityTemplateExact
    %problemElasticitySquareExactLocking - Divergence-free elasticity solution for locking tests
    %   OBJ = problemElasticitySquareExactLocking(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies an elasticity problem on the unit square (0,1)^2 with
    %   Dirichlet boundary on all edges and the divergence-free exact
    %   displacement
    %   u = pi*(cos(pi*y)*sin(pi*x)^2*sin(pi*y),
    %           -cos(pi*x)*sin(pi*x)*sin(pi*y)^2),
    %   which vanishes on the boundary. The body force and the stress do
    %   not depend on lambda, so the problem tests volumetric locking for
    %   nearly incompressible material (lambda >> mu). The material
    %   parameters NAME1,VALUE1,NAME2,VALUE2 are given as in
    %   ElasticityParameters and are required.
    %
    %   problemElasticitySquareExactLocking functions:
    %       f           - Body force
    %       exactU      - Exact displacement
    %       exactSigma  - Exact stress
    %       exactGradU  - Exact displacement gradient
    %       exactD2U    - Second derivatives of the solution
    %
    %   problemElasticitySquareExactLocking properties:
    %       geom        - Unit square geometry
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemElasticityTemplateExact, problemElasticitySquareExact,
    %   problemElasticitySquareExactZeroRhs

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        scale = 0.1
    end


    methods

        function val = f(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = [  -2*obj.mu*pi^3*cos(pi*y).*sin(pi*y).*(2*cos(2*pi*x) - 1) ,...
                2*obj.mu*pi^3*cos(pi*x).*sin(pi*x).*(2*cos(2*pi*y) - 1)];
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [ pi*cos(pi*y).*sin(pi*x).^2.*sin(pi*y),...
                -pi*cos(pi*x).*sin(pi*x).*sin(pi*y).^2 ];
        end

        function val = exactSigma(obj,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = 2*pi^2*cos(pi*x).*cos(pi*y).*sin(pi*x).*sin(pi*y).*(obj.lambda + 2*obj.mu) ...
                - 2*obj.lambda*pi^2*cos(pi*x).*cos(pi*y).*sin(pi*x).*sin(pi*y);
            val(:,2) = 2*obj.lambda*pi^2.*cos(pi*x).*cos(pi*y).*sin(pi*x).*sin(pi*y) ...
                - 2*pi^2*cos(pi*x).*cos(pi*y).*sin(pi*x).*sin(pi*y).*(obj.lambda + 2*obj.mu);
            val(:,3) = -(obj.mu*pi^2*(cos(2*pi*x) - cos(2*pi*y)))/2;
        end

        function val = exactGradU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = zeros(size(x,1),2,2);
            val(:,1,1) = 2*pi^2*cos(pi*x).*cos(pi*y).*sin(pi*x).*sin(pi*y);
            val(:,1,2) =  pi^2*cos(pi*y).^2.*sin(pi*x).^2 - pi^2*sin(pi*x).^2.*sin(pi*y).^2;
            val(:,2,1) =  pi^2*sin(pi*x).^2.*sin(pi*y).^2 - pi^2*cos(pi*x).^2.*sin(pi*y).^2;
            val(:,2,2) = -2*pi^2*cos(pi*x).*cos(pi*y).*sin(pi*x).*sin(pi*y);
        end

        function val = exactD2U(~,x)
            y = x(:,2);
            x = x(:,1);
            val = zeros(size(x,1),2,3);
            val(:,1,1) =  2*pi^3*cos(pi*x).^2.*cos(pi*y).*sin(pi*y) - 2*pi^3*cos(pi*y).*sin(pi*x).^2.*sin(pi*y);
            val(:,1,2) = -4*pi^3*cos(pi*y).*sin(pi*x).^2.*sin(pi*y);
            val(:,1,3) =  2*pi^3*cos(pi*x).*cos(pi*y).^2.*sin(pi*x) - 2*pi^3*cos(pi*x).*sin(pi*x).*sin(pi*y).^2;
            val(:,2,1) =  4*pi^3*cos(pi*x).*sin(pi*x).*sin(pi*y).^2;
            val(:,2,2) =  2*pi^3*cos(pi*x).*sin(pi*x).*sin(pi*y).^2 - 2*pi^3*cos(pi*x).*cos(pi*y).^2.*sin(pi*x);
            val(:,2,3) =  2*pi^3*cos(pi*y).*sin(pi*x).^2.*sin(pi*y) - 2*pi^3*cos(pi*x).^2.*cos(pi*y).*sin(pi*y);
        end

    end
end