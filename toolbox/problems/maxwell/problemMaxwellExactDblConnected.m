classdef problemMaxwellExactDblConnected
    %problemMaxwellExactDblConnected - Maxwell on a doubly connected domain
    %   problemMaxwellExactDblConnected is a Maxwell source problem on the
    %   doubly connected square ring (0,4)^2 \ [1,3]^2 with the
    %   manufactured polynomial exact field
    %   u = (y(1-y)(3-y)(4-y), x(1-x)(3-x)(4-x)), whose tangential trace
    %   vanishes on the outer and on the inner boundary.
    %
    %   OBJ = problemMaxwellExactDblConnected() uses ALPHA = 1.
    %
    %   OBJ = problemMaxwellExactDblConnected(ALPHA) sets the coefficient
    %   ALPHA of the zero-order term in curl(rot(u)) + ALPHA*u = f.
    %   ALPHA must be nonzero, since the domain is not simply connected.
    %
    %   problemMaxwellExactDblConnected functions:
    %       f           - Source current
    %       rotf        - Rot of the source current
    %       u4Db        - Dirichlet field (exact solution)
    %       u4Db2       - Indicator of the inner boundary (scalar)
    %       g           - Zero Neumann-type boundary data
    %       exactU      - Exact solution
    %       exactRotU   - Rot of the exact solution
    %
    %   problemMaxwellExactDblConnected properties:
    %       geom        - Square ring geometry
    %       m           - Number of holes of the domain (1)
    %       alpha       - Coefficient of the zero-order term
    %
    %   See also problemMaxwellExactSquare, problemMaxwellExactLshape,
    %   problemMaxwellSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','DblConnected'))
        m = 1
        alpha
    end


    methods

        function obj = problemMaxwellExactDblConnected(alpha)
            if nargin < 1 || isempty(alpha)
                alpha = 1;
            end
            obj.alpha = alpha;
        end

        function val = f(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = [12*y.^2 - 48*y + 38, -(48*x - 12*x.^2 - 38)];
            val = val + obj.alpha*[y.*(1-y).*(3-y).*(4-y),x.*(1-x).*(3-x).*(4-x)];
        end

        function val = rotf(obj,x)
            a = obj.alpha;
            y = x(:,2);
            x = x(:,1);
            val = (-2).*(x - y).*(2.*a.*x.^2 + 2.*a.*x.*y - 12.*a.*x + 2.*a.*y.^2 - 12.*a.*y + 19.*a - 12);
        end

        function val = u4Db(obj,x)
            val = obj.exactU(x);
        end

        function val = u4Db2(~,x)
            val = zeros(size(x,1),1);
            val(x(:,1)==1 & x(:,2)<=3 & x(:,2)>=1) = 1;
            val(x(:,1)==3 & x(:,2)<=3 & x(:,2)>=1) = 1;
            val(x(:,2)==1 & x(:,1)<=3 & x(:,1)>=1) = 1;
            val(x(:,2)==3 & x(:,1)<=3 & x(:,1)>=1) = 1;
        end

        function val = g(~,x)
            val = zeros(size(x,1),2);
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [y.*(1-y).*(3-y).*(4-y),x.*(1-x).*(3-x).*(4-x)];
        end

        function val = exactRotU(~,x)
            y = x(:,2);
            x = x(:,1);
            val= 24.*x.^2 - 4.*x.^3 - 38.*x + 4.*y.^3 - 24.*y.^2 + 38.*y;
        end
    end
end