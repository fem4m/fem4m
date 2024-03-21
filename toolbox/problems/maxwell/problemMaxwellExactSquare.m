classdef problemMaxwellExactSquare
    %problemMaxwellExactSquare - Maxwell on a square, trigonometric exact field
    %   problemMaxwellExactSquare is a Maxwell source problem on the unit
    %   square with the trigonometric exact field
    %   u = (sin(pi*y), sin(pi*x)), whose tangential trace vanishes on the
    %   boundary, used for convergence studies.
    %
    %   OBJ = problemMaxwellExactSquare() uses ALPHA = 1.
    %
    %   OBJ = problemMaxwellExactSquare(ALPHA) sets the coefficient ALPHA
    %   of the zero-order term in curl(rot(u)) + ALPHA*u = f.
    %
    %   problemMaxwellExactSquare functions:
    %       f           - Source current
    %       rotf        - Rot of the source current
    %       u4Db        - Dirichlet field (exact solution)
    %       g           - Zero Neumann-type boundary data
    %       exactU      - Exact solution
    %       exactRotU   - Rot of the exact solution
    %
    %   problemMaxwellExactSquare properties:
    %       geom        - Square geometry
    %       m           - Number of holes of the domain (0)
    %       alpha       - Coefficient of the zero-order term
    %
    %   See also problemMaxwellExactLshape, problemMaxwellExactDblConnected,
    %   problemMaxwellSquareMonomExact, problemMaxwellCubeExact2

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        m = 0
        alpha
    end


    methods

        function obj = problemMaxwellExactSquare(alpha)
            if nargin < 1 || isempty(alpha)
                alpha = 1;
            end
            obj.alpha = alpha;
        end


        function val = f(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = (pi^2 + obj.alpha)*[sin(pi*y),sin(pi*x)];
        end

        function val = rotf(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = (pi^2 + obj.alpha)*(-pi*cos(pi*y)+pi*cos(pi*x));
        end

        function val = u4Db(obj,x)
            val = obj.exactU(x);
        end

        function val = g(~,x)
            val = zeros(size(x,1),2);
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [sin(pi*y),sin(pi*x)];
        end

        function val = exactRotU(~,x)
            y = x(:,2);
            x = x(:,1);
            val= pi*cos(pi*x)-pi*cos(pi*y);
        end
    end
end