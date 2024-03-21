classdef problemHelmholtzSquareMonomExact < ProblemHelmholtzTemplateExact
    %problemHelmholtzSquareMonomExact - Helmholtz problem with a monomial
    %   problemHelmholtzSquareMonomExact is a Helmholtz problem
    %   -Laplace(u) - OMEGA^2*u = f on the unit square (0,1)^2 with the
    %   monomial exact solution u = x1^I*x2^J. The bottom and top edges
    %   carry Dirichlet data, the left and right edges impedance (Robin)
    %   data. It is meant for patch tests.
    %
    %   OBJ = problemHelmholtzSquareMonomExact(OMEGA,I,J) sets the angular
    %   frequency OMEGA and the monomial degrees I and J. All three
    %   arguments are required.
    %
    %   problemHelmholtzSquareMonomExact functions:
    %       f           - Source term -Laplace(u) - OMEGA^2*u
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemHelmholtzSquareMonomExact properties:
    %       geom        - Square geometry with Dirichlet and Robin edges
    %       i           - Monomial degree in x1
    %       j           - Monomial degree in x2
    %
    %   See also ProblemHelmholtzTemplateExact,
    %   problemScalarEllipticSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry',['D','N','D','N']))
        i
        j
    end

    methods
        function obj = problemHelmholtzSquareMonomExact(omega,i,j)
            obj = obj@ProblemHelmholtzTemplateExact(omega);
            obj.i = i;
            obj.j = j;
        end

        function val = f(obj,x)
            val = -(monomial(x,obj.i,obj.j,2,0) + monomial(x,obj.i,obj.j,0,2)) ...
                - obj.omega^2*monomial(x,obj.i,obj.j,0,0);
        end

        function val = exactU(obj,x)
            val = monomial(x,obj.i,obj.j,0,0);
        end

        function val = exactGradU(obj,x)
            val = [monomial(x,obj.i,obj.j,1,0),monomial(x,obj.i,obj.j,0,1)];
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