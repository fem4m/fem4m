classdef problemTransportSquare < ProblemScalarEllipticTemplate
    %problemTransportSquare - Transport on (-1,1)^2 with oscillating inflow data
    %   problemTransportSquare is the transport problem beta*grad(u) = f on
    %   the square (-1,1)^2 with the velocity
    %   beta = (2, 1+4/5*sin(8*pi*x)), whose inflow sides are x = -1 and
    %   y = -1. With r the distance to the point (-3/4,-3/4), the source f
    %   is 16 for r < 1/4 and 0 elsewhere, and the oscillating Dirichlet
    %   data is exp(5*(1-r^2))*sin(16*pi*r^2).
    %
    %   problemTransportSquare functions:
    %       f     - Discontinuous source term
    %       u4Db  - Oscillating Dirichlet inflow data
    %       g     - Zero Neumann data
    %       beta  - Convection velocity (2, 1+4/5*sin(8*pi*x))
    %
    %   problemTransportSquare properties:
    %       geom  - Square geometry (-1,1)^2
    %
    %   See also ProblemScalarEllipticTemplate, problemTransportSquare2,
    %   problemTransportSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','midpoint',[0,0],'length',2))
    end

    methods

        function val = f(~,x)
            y = x(:,2);
            x = x(:,1);
            r = sqrt((x+3/4).^2+(y+3/4).^2);
            val = zeros(size(x,1),1);
            val(r<1/4) = 1/(1/4^2);
        end

        function val = u4Db(~,x)
            y = x(:,2);
            x = x(:,1);
            r = sqrt((x+3/4).^2+(y+3/4).^2);
            val = exp(5*(1-r.^2)).*sin(16*pi*r.^2);
        end

        function val = g(~,x,~)
            val = zeros(size(x,1),1);
        end

        function val = beta(~,x)
            val = [2*ones(size(x,1),1),1+4/5*sin(8*pi*x(:,1))];
        end

    end

end