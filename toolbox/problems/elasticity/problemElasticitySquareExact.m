classdef problemElasticitySquareExact < ProblemElasticityTemplateExact
    %problemElasticitySquareExact - Elasticity on the unit square
    %   OBJ = problemElasticitySquareExact(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies a linear elasticity problem on the unit square (0,1)^2
    %   with Dirichlet boundary on all edges and the polynomial exact
    %   displacement u = x*y*(1-x)*(1-y)*(1,1), which vanishes on the
    %   boundary. The material parameters NAME1,VALUE1,NAME2,VALUE2 are
    %   given as in ElasticityParameters and are required.
    %
    %   problemElasticitySquareExact functions:
    %       f                - Body force
    %       exactU           - Exact displacement
    %       exactSigma       - Exact stress
    %       exactGradU       - Exact displacement gradient
    %       exactRotSigma    - Row-wise rot of the stress
    %       exactRotRotSigma - Rot of the rot of the stress
    %       exactRotEpsilon  - Row-wise rot of the strain
    %       exactD2U         - Second derivatives of the solution
    %
    %   problemElasticitySquareExact properties:
    %       geom             - Unit square geometry
    %       scale            - Amplification of the displacement plots
    %
    %   See also ProblemElasticityTemplateExact,
    %   problemElasticitySquareExactLocking,
    %   problemElasticitySquareNeumannExact,
    %   problemElasticitySquareMonomExact, problemElasticityLshapeExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        scale = 2
    end


    methods

        function val = f(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = [- 2.*obj.mu.*(2.*x.*y - y - 2.*x + x.^2 + 1/2) ...
                - 2.*y.*(obj.lambda + 2.*obj.mu).*(y - 1) ...
                - obj.lambda.*x.*(2.*y - 1) - obj.lambda.*(2.*y - 1).*(x - 1) ,...
                - 2.*obj.mu.*(2.*x.*y - 2.*y - x + y.^2 + 1/2) ...
                - 2.*x.*(obj.lambda + 2.*obj.mu).*(x - 1)...
                - obj.lambda.*y.*(2.*x - 1) - obj.lambda.*(2.*x - 1).*(y - 1)];
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [x.*y.*(1-x).*(1-y),x.*y.*(1-x).*(1-y)];
        end

        function val = exactSigma(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = [obj.lambda.*x.*(2.*y - 1).*(x - 1) + y.*(2.*x - 1).*(obj.lambda + 2.*obj.mu).*(y - 1) , obj.lambda.*y.*(2.*x - 1).*(y - 1) + x.*(2.*y - 1).*(obj.lambda + 2.*obj.mu).*(x - 1) , 2.*obj.mu.*(x.^2.*y - x.^2/2 + x.*y.^2 - 2.*x.*y + x/2 - y.^2/2 + y/2)];
        end

        function val = exactGradU(~,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1,1) = x.*y.*(y - 1) + y.*(x - 1).*(y - 1);
            val(:,1,2) = x.*y.*(x - 1) + x.*(x - 1).*(y - 1);
            val(:,2,1) = val(:,1,1);
            val(:,2,2) = val(:,1,2);
        end

        function val = exactRotSigma(obj,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = 2.*obj.mu.*(2.*x.*y - 2.*y - x + y.^2 + 1/2) - (2.*x - 1).*(obj.lambda + 2.*obj.mu).*(y - 1) - 2.*obj.lambda.*x.*(x - 1) - y.*(2.*x - 1).*(obj.lambda + 2.*obj.mu);
            val(:,2) = (2.*y - 1).*(obj.lambda + 2.*obj.mu).*(x - 1) - 2.*obj.mu.*(2.*x.*y - y - 2.*x + x.^2 + 1/2) + 2.*obj.lambda.*y.*(y - 1) + x.*(2.*y - 1).*(obj.lambda + 2.*obj.mu);
        end

        function val = exactRotRotSigma(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = 2*(2*x - 1)*(obj.lambda + 2*obj.mu) + 2*(2*y - 1)*(obj.lambda + 2*obj.mu) - 4*obj.mu*(2*x + 2*y - 2);
        end

        function val = exactRotEpsilon(~,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = 2.*x.*y - 2.*y - y.*(2.*x - 1) - x - (2.*x - 1).*(y - 1) + y.^2 + 1/2;
            val(:,2) = 2.*x + y + x.*(2.*y - 1) - 2.*x.*y + (2.*y - 1).*(x - 1) - x.^2 - 1/2;
        end

        function val = exactD2U(~,x)
            y = x(:,2);
            x = x(:,1);
            val = zeros(size(x,1),2,3);
            val(:,1,1) = -2*y.*(1-y);
            val(:,1,2) = -2*x.*(1-x);
            val(:,1,3) = (1-2*x).*(1-2*y);
            val(:,2,:) = val(:,1,:);
        end
    end
end