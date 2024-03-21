classdef problemScalarEllipticCubeMonomExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticCubeMonomExact - Elliptic problem with a monomial on a cube
    %   problemScalarEllipticCubeMonomExact is a Poisson problem on the unit
    %   cube with the monomial exact solution u = x^I y^J z^L and
    %   inhomogeneous Dirichlet data. It is meant for patch tests.
    %
    %   OBJ = problemScalarEllipticCubeMonomExact(I,J,L) sets the monomial
    %   degrees I, J and L.
    %
    %   problemScalarEllipticCubeMonomExact functions:
    %       f           - Source term of the monomial solution
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticCubeMonomExact properties:
    %       geom        - Unit cube geometry
    %       i           - Monomial degree in x
    %       j           - Monomial degree in y
    %       l           - Monomial degree in z
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Cube','midpoint',[0.5,0.5,0.5],'length',1))
        i
        j
        l
    end

    methods
        function obj = problemScalarEllipticCubeMonomExact(i,j,l)
            obj.i = i;
            obj.j = j;
            obj.l = l;
        end

        function val = f(obj,x)
            val = -(monomial(x,obj.i,obj.j,obj.l,2,0,0) ...
                + monomial(x,obj.i,obj.j,obj.l,0,2,0) ...
                + monomial(x,obj.i,obj.j,obj.l,0,0,2));
        end

        function val = exactU(obj,x)
            val = monomial(x,obj.i,obj.j,obj.l,0,0,0);
        end

        function val = exactGradU(obj,x)
            val = [monomial(x,obj.i,obj.j,obj.l,1,0,0),...
                monomial(x,obj.i,obj.j,obj.l,0,1,0),...
                monomial(x,obj.i,obj.j,obj.l,0,0,1)];
        end

    end

end

function val = monomial(x,i,j,l,dx,dy,dz)
%MONOMIAL - Derivative d^(dx,dy,dz) of x^i y^j z^l
val = derivative1d(x(:,1),i,dx).*derivative1d(x(:,2),j,dy).*derivative1d(x(:,3),l,dz);
end

function val = derivative1d(x,k,d)
%derivative1d - d-th derivative of x^k
if d > k
    val = zeros(size(x));
else
    val = prod(k-d+1:k)*x.^(k-d);
end
end