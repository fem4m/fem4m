classdef problemTransportSquareExact < ProblemScalarEllipticTemplateExact
    %problemTransportSquareExact - Transport with a manufactured solution
    %   problemTransportSquareExact is the transport problem
    %   beta*grad(u) = f on the unit square (0,1)^2 with the constant
    %   velocity beta = (1,0) and the exact solution
    %   u = x*(1-x)*y*(1-y), which vanishes on the boundary (homogeneous
    %   inflow data).
    %
    %   problemTransportSquareExact functions:
    %       f           - Source term beta*grad(u)
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %       beta        - Constant convection velocity (1,0)
    %       alpha       - Zero diffusion (pure transport)
    %
    %   problemTransportSquareExact properties:
    %       geom        - Unit square geometry
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemTransportSquareMonomExact, problemTransportSquare

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
    end

    methods

        function val = f(obj,x)
            val = matMul(permute(obj.beta(x),[1 3 2]),obj.exactGradU(x));
        end

        function val = exactU(~,x)
            val = x(:,1).*(1-x(:,1)).*x(:,2).*(1-x(:,2));
        end

        function val = exactGradU(~,x)
            val = [ (1 - 2*x(:,1)).*x(:,2).*(1-x(:,2)),...
                x(:,1).*(1-x(:,1)).*(1-2*x(:,2))];
        end

        function val = alpha(~,x)
            d = size(x,2);
            val = zeros(1,d,d);
        end

        function val = beta(~,x)
            val = zeros(size(x,1),2);
            val(:,1) = 1;
            val(:,2) = 0;
        end


    end

end