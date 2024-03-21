classdef problemElasticitySquareNeumannExact < ProblemElasticityTemplateExact
    %problemElasticitySquareNeumannExact - Elasticity with Neumann boundary
    %   OBJ = problemElasticitySquareNeumannExact(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies an elasticity problem on the unit square (0,1)^2 with
    %   Neumann boundary on the edges y = 0, x = 1 and y = 1, Dirichlet
    %   boundary on the left edge x = 0, and the trigonometric exact
    %   displacement u = -(sin(pi*x)*cos(pi*y)^2, cos(pi*x)*sin(pi*y)^2).
    %   The material parameters NAME1,VALUE1,NAME2,VALUE2 are given as in
    %   ElasticityParameters and are required.
    %
    %   problemElasticitySquareNeumannExact functions:
    %       f                - Body force
    %       exactU           - Exact displacement
    %       exactSigma       - Exact stress
    %       exactGradU       - Exact displacement gradient
    %       exactRotEpsilon  - Row-wise rot of the strain
    %       exactD2U         - Second derivatives of the solution
    %
    %   problemElasticitySquareNeumannExact properties:
    %       geom             - Unit square geometry
    %       scale            - Amplification of the displacement plots
    %
    %   See also ProblemElasticityTemplateExact, problemElasticitySquareExact,
    %   problemElasticitySquareExactZeroRhs

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry',['N','N','N','D']))
        scale = 0.1
    end


    methods

        function val = f(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = (pi.^2/2).*[obj.lambda.*cos(pi.*x + 2.*pi.*y) - obj.lambda.*cos(pi.*x - 2.*pi.*y) - 2.*obj.mu.*sin(pi.*x).*sin(2.*pi.*y) - 2.*obj.lambda.*cos(pi.*y).^2.*sin(pi.*x) - 4.*obj.mu.*cos(pi.*y).^2.*sin(pi.*x) - 4.*obj.mu.*cos(2.*pi.*y).*sin(pi.*x),...
                2.*obj.lambda.*cos(pi.*x - 2.*pi.*y) - obj.mu.*cos(pi.*x) + 2.*obj.lambda.*cos(pi.*x + 2.*pi.*y) + 4.*obj.mu.*cos(pi.*x - 2.*pi.*y) + 4.*obj.mu.*cos(pi.*x + 2.*pi.*y) + obj.mu.*cos(pi.*x).*cos(2.*pi.*y) - 2.*obj.mu.*cos(pi.*x).*sin(2.*pi.*y) - 4.*obj.lambda.*cos(pi.*x).*cos(pi.*y).*sin(pi.*y)];
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [cos(pi.*(y + 1)).^2.*sin(pi.*(x + 1)),...
                cos(pi.*(x + 1)).*sin(pi.*(y + 1)).^2];
        end

        function val = exactSigma(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = (pi/2).*[obj.lambda.*sin(pi.*x - 2.*pi.*y) - obj.lambda.*sin(pi.*x + 2.*pi.*y) - 2.*obj.lambda.*cos(pi.*x).*cos(pi.*y).^2 - 4.*obj.mu.*cos(pi.*x).*cos(pi.*y).^2 , ...
                obj.lambda.*sin(pi.*x - 2.*pi.*y) - 2.*obj.lambda.*cos(pi.*x).*cos(pi.*y).^2 - obj.lambda.*sin(pi.*x + 2.*pi.*y) + 2.*obj.mu.*sin(pi.*x - 2.*pi.*y) - 2.*obj.mu.*sin(pi.*x + 2.*pi.*y),...
                obj.mu.*sin(pi.*x) + 2.*obj.mu.*sin(pi.*x).*sin(2.*pi.*y) - obj.mu.*cos(2.*pi.*y).*sin(pi.*x)];
        end

        function val = exactGradU(~,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1,1) = pi.*cos(pi.*(x + 1)).*cos(pi.*(y + 1)).^2;
            val(:,1,2) = (-2).*pi.*cos(pi.*(y + 1)).*sin(pi.*(x + 1)).*sin(pi.*(y + 1));
            val(:,2,1) = -pi.*sin(pi.*(x + 1)).*sin(pi.*(y + 1)).^2;
            val(:,2,2) = 2.*pi.*cos(pi.*(x + 1)).*cos(pi.*(y + 1)).*sin(pi.*(y + 1));
        end


        function val = exactRotEpsilon(~,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = -pi.^2.*(cos(pi.*(x - 2.*y))./8 - cos(pi.*x)./4 + cos(pi.*(x + 2.*y))./8 - sin(pi.*(x - 2.*y))./4 + sin(pi.*(x + 2.*y))./4);
            val(:,2) = -pi.^2.*(cos(pi.*(x + 2.*y))./4 - cos(pi.*(x - 2.*y))./4 + sin(pi.*(x - 2.*y))./2 + sin(pi.*(x + 2.*y))./2);
        end

        function val = exactD2U(~,x)
            y = x(:,2);
            x = x(:,1);
            val = zeros(size(x,1),2,3);
            val(:,1,1) =  pi.^2.*cos(y.*pi).^2.*sin(x.*pi);
            val(:,1,2) =  pi.^2.*cos(y.*pi.*2.0).*sin(x.*pi).*2.0;
            val(:,1,3) =  pi.^2.*cos(x.*pi).*cos(y.*pi).*sin(y.*pi).*2.0;
            val(:,2,1) =  pi.^2.*cos(x.*pi).*sin(y.*pi).^2;
            val(:,2,2) =  pi.^2.*cos(x.*pi).*cos(y.*pi.*2.0).*-2.0;
            val(:,2,3) =  pi.^2.*cos(y.*pi).*sin(x.*pi).*sin(y.*pi).*2.0;
        end

    end
end