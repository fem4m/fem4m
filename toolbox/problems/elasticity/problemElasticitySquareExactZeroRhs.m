classdef problemElasticitySquareExactZeroRhs < ProblemElasticityTemplateExact
    %problemElasticitySquareExactZeroRhs - Elasticity with zero load for H = K = 0
    %   OBJ = problemElasticitySquareExactZeroRhs(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies an elasticity problem on the unit square (0,1)^2
    %   (Dirichlet on the left edge x = 0, Neumann elsewhere) with a
    %   quintic polynomial exact displacement plus the terms
    %   (-H*x^2/2,-H*y^2/2) + K*(x^3+y^3)*(1,1) with the defaults H = 1
    %   and K = 0. The material parameters NAME1,VALUE1,NAME2,VALUE2 are
    %   given as in ElasticityParameters and are required. The quintic
    %   part produces no load, so only the H and K terms enter the body
    %   force: it vanishes for H = K = 0, and the defaults give the
    %   constant body force (lambda+2*mu)*(1,1).
    %
    %   OBJ = problemElasticitySquareExactZeroRhs(...,H) also sets the
    %   coefficient H of the quadratic terms.
    %
    %   OBJ = problemElasticitySquareExactZeroRhs(...,H,K) also sets the
    %   coefficient K of the cubic terms. An empty H or K selects its
    %   default.
    %
    %   problemElasticitySquareExactZeroRhs functions:
    %       f           - Body force (zero for H = K = 0)
    %       exactU      - Exact displacement
    %       exactSigma  - Exact stress
    %       exactGradU  - Exact displacement gradient
    %       exactD2U    - Second derivatives of the solution
    %
    %   problemElasticitySquareExactZeroRhs properties:
    %       geom        - Unit square geometry
    %       scale       - Amplification of the displacement plots
    %       h           - Coefficient H of the quadratic terms (default 1)
    %       k           - Coefficient K of the cubic terms (default 0)
    %
    %   See also ProblemElasticityTemplateExact, problemElasticitySquareExact,
    %   problemElasticitySquareNeumannExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry',['N','N','N','D']))
        scale = 0.1
        h
        k
    end


    methods

        function obj = problemElasticitySquareExactZeroRhs(name1,value1,name2,value2,h,k)
            if nargin < 4
                [name1,value1,name2,value2] = deal([]);
            end
            obj = obj@ProblemElasticityTemplateExact(name1,value1,name2,value2);
            if nargin < 5 || isempty(h)
                obj.h = 1;
            else
                obj.h = h;
            end
            if nargin < 6 || isempty(k)
                obj.k = 0;
            else
                obj.k = k;
            end
        end

        function val = f(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = zeros(size(x,1),2);
            val(:,1) = obj.lambda*obj.h + 2*obj.mu*obj.h - 6*obj.lambda*obj.k*x - 12*obj.mu*obj.k*x - 6*obj.mu*obj.k*y;
            val(:,2) = obj.lambda*obj.h + 2*obj.mu*obj.h - 6*obj.mu*obj.k*x - 6*obj.lambda*obj.k*y - 12*obj.mu*obj.k*y;
        end

        function val = exactU(obj,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = - x.^5/20 + (x.^4.*y)/4 + (x.^3.*y.^2)/2 + obj.k.*x.^3 - (x.^2.*y.^3)/2 - (obj.h.*x.^2)/2 - (x.*y.^4)/4 + y.^5/20 + obj.k.*y.^3;
            val(:,2) = x.^5/20 + (x.^4.*y)/4 - (x.^3.*y.^2)/2 + obj.k.*x.^3 - (x.^2.*y.^3)/2 + (x.*y.^4)/4 + y.^5/20 + obj.k.*y.^3 - (obj.h.*y.^2)/2;
        end

        function val = exactSigma(obj,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = obj.lambda.*(x.^4/4 - x.^3.*y - (3.*x.^2.*y.^2)/2 + x.*y.^3 + y.^4/4 + 3.*obj.k.*y.^2 - obj.h.*y) - (obj.lambda + 2.*obj.mu).*(x.^4/4 - x.^3.*y - (3.*x.^2.*y.^2)/2 - 3.*obj.k.*x.^2 + x.*y.^3 + obj.h.*x + y.^4/4);
            val(:,2) = (obj.lambda + 2.*obj.mu).*(x.^4/4 - x.^3.*y - (3.*x.^2.*y.^2)/2 + x.*y.^3 + y.^4/4 + 3.*obj.k.*y.^2 - obj.h.*y) - obj.lambda.*(x.^4/4 - x.^3.*y - (3.*x.^2.*y.^2)/2 - 3.*obj.k.*x.^2 + x.*y.^3 + obj.h.*x + y.^4/4);
            val(:,3) = 2.*obj.mu.*(x.^4/4 + x.^3.*y - (3.*x.^2.*y.^2)/2 + (3.*obj.k.*x.^2)/2 - x.*y.^3 + y.^4/4 + (3.*obj.k.*y.^2)/2);
        end

        function val = exactGradU(obj,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1,1) = - x.^4/4 + x.^3.*y + (3.*x.^2.*y.^2)/2 + 3.*obj.k.*x.^2 - x.*y.^3 - obj.h.*x - y.^4/4;
            val(:,1,2) = x.^4/4 + x.^3.*y - (3.*x.^2.*y.^2)/2 - x.*y.^3 + y.^4/4 + 3.*obj.k.*y.^2;
            val(:,2,1) = x.^4/4 + x.^3.*y - (3.*x.^2.*y.^2)/2 + 3.*obj.k.*x.^2 - x.*y.^3 + y.^4/4;
            val(:,2,2) = x.^4/4 - x.^3.*y - (3.*x.^2.*y.^2)/2 + x.*y.^3 + y.^4/4 + 3.*obj.k.*y.^2 - obj.h.*y;
        end

        function val = exactD2U(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = zeros(size(x,1),2,3);
            val(:,1,1) = - x.^3 + 3*x.^2.*y + 3*x.*y.^2 + 6*obj.k*x - y.^3 - obj.h;
            val(:,1,2) = x.^3 - 3*x.^2.*y - 3*x.*y.^2 + y.^3 + 6*obj.k*y;
            val(:,1,3) = x.^3 + 3*x.^2.*y - 3*x.*y.^2 - y.^3;
            val(:,2,1) = x.^3 + 3*x.^2.*y - 3*x.*y.^2 + 6*obj.k*x - y.^3;
            val(:,2,2) = - x.^3 - 3*x.^2.*y + 3*x.*y.^2 + y.^3 + 6*obj.k*y - obj.h;
            val(:,2,3) = x.^3 - 3*x.^2.*y - 3*x.*y.^2 + y.^3;
        end

    end
end