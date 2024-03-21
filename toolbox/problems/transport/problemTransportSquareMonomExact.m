classdef problemTransportSquareMonomExact < ProblemScalarEllipticTemplateExact
    %problemTransportSquareMonomExact - Transport with a monomial solution
    %   problemTransportSquareMonomExact is the transport problem
    %   beta*grad(u) = f on the unit square (0,1)^2 with the constant
    %   velocity beta = (1,1/2) and the monomial exact solution
    %   u = x^I*y^J. The inflow sides are x = 0 and y = 0; the inflow data
    %   u is nonzero there only for I = 0 or J = 0. It is meant for patch
    %   tests.
    %
    %   OBJ = problemTransportSquareMonomExact(I,J) sets the monomial
    %   degrees I and J (nonnegative integers); both are required.
    %
    %   problemTransportSquareMonomExact functions:
    %       f           - Source term beta*grad(u)
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %       beta        - Constant convection velocity (1,1/2)
    %       alpha       - Zero diffusion (pure transport)
    %
    %   problemTransportSquareMonomExact properties:
    %       geom        - Unit square geometry
    %       i           - Monomial degree in x
    %       j           - Monomial degree in y
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemTransportSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        i
        j
    end

    methods

        function obj = problemTransportSquareMonomExact(i,j)
            obj.i = i;
            obj.j = j;
        end

        function val = f(obj,x)
            val = sum(obj.beta(x).*obj.exactGradU(x),2);
        end

        function val = exactU(obj,x)
            val = monomial(x,obj.i,obj.j,0,0);
        end

        function val = exactGradU(obj,x)
            val = [monomial(x,obj.i,obj.j,1,0),monomial(x,obj.i,obj.j,0,1)];
        end

        function val = alpha(~,x)
            d = size(x,2);
            val = zeros(1,d,d);
        end

        function val = beta(~,x)
            val = repmat([1,1/2],size(x,1),1);
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