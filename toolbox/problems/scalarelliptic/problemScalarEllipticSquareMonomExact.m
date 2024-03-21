classdef problemScalarEllipticSquareMonomExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticSquareMonomExact - Elliptic problem with a monomial
    %   problemScalarEllipticSquareMonomExact is a Poisson problem on the
    %   unit square with the monomial exact solution u = x^I y^J and
    %   inhomogeneous Dirichlet data. It is meant for patch tests.
    %
    %   OBJ = problemScalarEllipticSquareMonomExact(I,J) sets the monomial
    %   degrees I and J.
    %
    %   problemScalarEllipticSquareMonomExact functions:
    %       f           - Source term of the monomial solution
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticSquareMonomExact properties:
    %       geom        - Unit square geometry
    %       i           - Monomial degree in x
    %       j           - Monomial degree in y
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticCubeMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        i
        j
    end

    methods
        function obj = problemScalarEllipticSquareMonomExact(i,j)
            obj.i = i;
            obj.j = j;
        end

        function val = f(obj,x)
            val = -(obj.i.*(obj.i-1).*x(:,1).^(max(0,obj.i-2)).*x(:,2).^(obj.j) ...
                + obj.j.*(obj.j-1).*x(:,1).^(obj.i).*x(:,2).^(max(0,obj.j-2)));
        end

        function val = exactU(obj,x)
            val = x(:,1).^(obj.i).*x(:,2).^(obj.j);
        end

        function val = exactGradU(obj,x)
            val = [ obj.i.*x(:,1).^(max(0,obj.i-1)).*x(:,2).^(obj.j),...
                obj.j.*x(:,1).^(obj.i).*x(:,2).^(max(0,obj.j-1))];
        end

    end

end