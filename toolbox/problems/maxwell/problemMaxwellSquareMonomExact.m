classdef problemMaxwellSquareMonomExact
    %problemMaxwellSquareMonomExact - Maxwell on a square, monomial field
    %   problemMaxwellSquareMonomExact is a Maxwell source problem
    %   curl(rot(u)) + ALPHA*u = f on the unit square with the polynomial
    %   exact field u = (x^I y^J, x^J y^I) and inhomogeneous tangential
    %   Dirichlet data. It is meant for patch tests.
    %
    %   OBJ = problemMaxwellSquareMonomExact(ALPHA,I,J) sets the coefficient
    %   ALPHA of the zero-order term and the exponents I >= 0 and J >= 0
    %   (integers) of the exact field.
    %
    %   problemMaxwellSquareMonomExact functions:
    %       f           - Source current curl(rot(u)) + ALPHA*u
    %       rotf        - Rot of the source current
    %       u4Db        - Dirichlet field (exact solution)
    %       g           - Zero Neumann-type boundary data
    %       exactU      - Exact solution
    %       exactRotU   - Rot of the exact solution
    %
    %   problemMaxwellSquareMonomExact properties:
    %       geom        - Square geometry
    %       alpha       - Coefficient of the zero-order term
    %       i           - Exponent I in u = (x^I y^J, x^J y^I)
    %       j           - Exponent J in u = (x^I y^J, x^J y^I)
    %
    %   See also problemMaxwellExactSquare, Nd1k, Nd2k

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        alpha
        i
        j
    end

    methods

        function obj = problemMaxwellSquareMonomExact(alpha,i,j)
            obj.alpha = alpha;
            obj.i = i;
            obj.j = j;
        end

        function val = f(obj,x)
            %F - Source current curl(rot(u)) + ALPHA*u
            %   curl(rot(u)) = (d_xy u2 - d_yy u1, d_xy u1 - d_xx u2)
            val = [obj.u2(x,1,1) - obj.u1(x,0,2),...
                obj.u1(x,1,1) - obj.u2(x,2,0)] + obj.alpha*obj.exactU(x);
        end

        function val = rotf(obj,x)
            val = obj.u1(x,2,1) - obj.u2(x,3,0) - obj.u2(x,1,2) + obj.u1(x,0,3) ...
                + obj.alpha*obj.exactRotU(x);
        end

        function val = u4Db(obj,x)
            val = obj.exactU(x);
        end

        function val = g(~,x)
            val = zeros(size(x,1),2);
        end

        function val = exactU(obj,x)
            val = [obj.u1(x,0,0),obj.u2(x,0,0)];
        end

        function val = exactRotU(obj,x)
            val = obj.u2(x,1,0) - obj.u1(x,0,1);
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