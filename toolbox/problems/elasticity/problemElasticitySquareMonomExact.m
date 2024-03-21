classdef problemElasticitySquareMonomExact < ProblemElasticityTemplateExact
    %problemElasticitySquareMonomExact - Elasticity with a monomial solution
    %   problemElasticitySquareMonomExact is a linear elasticity problem on
    %   the unit square (0,1)^2 with the polynomial exact displacement
    %   u = (x^I*y^J, x^J*y^I) and inhomogeneous boundary data (Dirichlet
    %   on all edges by default). It is meant for patch tests; set geom to
    %   change the boundary conditions.
    %
    %   OBJ = problemElasticitySquareMonomExact(NAME1,VALUE1,NAME2,VALUE2,I,J)
    %   sets the material parameters as in ElasticityParameters and the
    %   monomial degrees I and J (nonnegative integers). All arguments are
    %   required.
    %
    %   problemElasticitySquareMonomExact functions:
    %       f           - Body force -div(sigma)
    %       exactU      - Exact displacement
    %       exactSigma  - Exact stress [s11 s22 s12]
    %       exactGradU  - Exact displacement gradient
    %       exactD2U    - Second derivatives of the displacement
    %
    %   problemElasticitySquareMonomExact properties:
    %       geom        - Unit square geometry
    %       scale       - Amplification of the displacement plots
    %       i           - Exponent of x in u1 and of y in u2
    %       j           - Exponent of y in u1 and of x in u2
    %
    %   See also ProblemElasticityTemplateExact, problemElasticitySquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        scale = 1
        i
        j
    end

    methods

        function obj = problemElasticitySquareMonomExact(name1,value1,name2,value2,i,j)
            obj = obj@ProblemElasticityTemplateExact(name1,value1,name2,value2);
            obj.i = i;
            obj.j = j;
        end

        function val = f(obj,x)
            lam = obj.lambda;
            mu = obj.mu;
            val = -[(lam+2*mu)*obj.u1(x,2,0) + mu*obj.u1(x,0,2) + (lam+mu)*obj.u2(x,1,1),...
                mu*obj.u2(x,2,0) + (lam+2*mu)*obj.u2(x,0,2) + (lam+mu)*obj.u1(x,1,1)];
        end

        function val = exactU(obj,x)
            val = [obj.u1(x,0,0),obj.u2(x,0,0)];
        end

        function val = exactSigma(obj,x)
            eps11 = obj.u1(x,1,0);
            eps22 = obj.u2(x,0,1);
            eps12 = (obj.u1(x,0,1) + obj.u2(x,1,0))/2;
            val = [(obj.lambda+2*obj.mu)*eps11 + obj.lambda*eps22,...
                obj.lambda*eps11 + (obj.lambda+2*obj.mu)*eps22,...
                2*obj.mu*eps12];
        end

        function val = exactGradU(obj,x)
            val = zeros(size(x,1),2,2);
            val(:,1,1) = obj.u1(x,1,0);
            val(:,1,2) = obj.u1(x,0,1);
            val(:,2,1) = obj.u2(x,1,0);
            val(:,2,2) = obj.u2(x,0,1);
        end

        function val = exactD2U(obj,x)
            val = zeros(size(x,1),2,3);
            val(:,1,:) = [obj.u1(x,2,0),obj.u1(x,0,2),obj.u1(x,1,1)];
            val(:,2,:) = [obj.u2(x,2,0),obj.u2(x,0,2),obj.u2(x,1,1)];
        end

    end

    methods (Access = private)

        function val = u1(obj,x,dx,dy)
            val = monomial(x,obj.i,obj.j,dx,dy);
        end

        function val = u2(obj,x,dx,dy)
            val = monomial(x,obj.j,obj.i,dx,dy);
        end

    end

end

function val = monomial(x,i,j,dx,dy)
%MONOMIAL - Derivative d^(dx,dy) of x^i y^j
val = derivative1d(x(:,1),i,dx).*derivative1d(x(:,2),j,dy);
end

function val = derivative1d(x,k,d)
%derivative1d - d-th derivative of x^k
if d > k
    val = zeros(size(x));
else
    val = prod(k-d+1:k)*x.^(k-d);
end
end