classdef problemConvectionDominatedSquareExact < ProblemScalarEllipticTemplateExact
    %problemConvectionDominatedSquareExact - Convection-diffusion problem with boundary layers
    %   problemConvectionDominatedSquareExact is the convection-diffusion
    %   problem -A*div(grad(u)) + beta*grad(u) = f on the unit square
    %   (0,1)^2 with the constant velocity beta = (1,1) and the diffusivity
    %   A. The exact solution u = v(x)*v(y) with
    %   v(t) = (exp((t-1)/A)-1)/(exp(-1/A)-1) + t - 1 vanishes on the
    %   boundary and has boundary layers at x = 1 and y = 1; the problem
    %   is convection-dominated for small A.
    %
    %   OBJ = problemConvectionDominatedSquareExact() uses A = 1.
    %
    %   OBJ = problemConvectionDominatedSquareExact(A) sets the diffusivity
    %   A > 0; an empty A keeps the default 1.
    %
    %   problemConvectionDominatedSquareExact functions:
    %       f           - Source term -A*div(grad(u)) + beta*grad(u)
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %       alpha       - Diffusion coefficient A*I
    %       beta        - Constant convection velocity (1,1)
    %
    %   problemConvectionDominatedSquareExact properties:
    %       geom        - Unit square geometry
    %       a           - Diffusivity A (default 1)
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemTransportSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        a = 1
    end

    methods

        function obj = problemConvectionDominatedSquareExact(a)
            if nargin > 0 && ~isempty(a)
                obj.a = a;
            end
        end

        function val = f(obj,x)
            alpha = obj.a;
            y = x(:,2);
            x = x(:,1);

            val = (exp((y - 1)/alpha)./(alpha*(exp(-1/alpha) - 1)) + 1).*(x + (exp((x - 1)/alpha) - 1)./(exp(-1/alpha) - 1) - 1) - alpha*((exp((y - 1)/alpha).*(x + (exp((x - 1)/alpha) - 1)./(exp(-1/alpha) - 1) - 1))./(alpha^2*(exp(-1/alpha) - 1)) + (exp((x - 1)/alpha).*(y + (exp((y - 1)/alpha) - 1)/(exp(-1/alpha) - 1) - 1))./(alpha^2*(exp(-1/alpha) - 1))) + (exp((x - 1)/alpha)./(alpha*(exp(-1/alpha) - 1)) + 1).*(y + (exp((y - 1)/alpha) - 1)./(exp(-1/alpha) - 1) - 1);

        end

        function val = exactU(obj,x)
            alpha = obj.a;
            y = x(:,2);
            x = x(:,1);

            val1 = ( exp((x-1)/alpha) -1 )./ ( exp(-1/alpha) -1 )  + x - 1;
            val2 = ( exp((y-1)/alpha) -1 )./ ( exp(-1/alpha) -1 )  + y - 1;

            val = val1 .* val2;
        end

        function val = exactGradU(obj,x)
            alpha = obj.a;
            y = x(:,2);
            x = x(:,1);

            valdx = (exp((x - 1)/alpha)./(alpha*(exp(-1/alpha) - 1)) + 1).*(y + (exp((y - 1)/alpha) - 1)./(exp(-1/alpha) - 1) - 1);
            valdy = (exp((y - 1)/alpha)./(alpha*(exp(-1/alpha) - 1)) + 1).*(x + (exp((x - 1)/alpha) - 1)./(exp(-1/alpha) - 1) - 1);

            val = [valdx,valdy];
        end

        function val = alpha(obj,x)
            val = zeros(size(x,1),2,2);
            val(:,1,1) = obj.a;
            val(:,2,2) = obj.a;
        end

        function val = beta(~,x)
            val = ones(size(x,1),2);
        end


    end

end