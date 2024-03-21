classdef problemTransportSquare2 < ProblemScalarEllipticTemplate
    %problemTransportSquare2 - Transport on the unit square with discontinuous inflow data
    %   problemTransportSquare2 is the transport problem beta*grad(u) = 0
    %   on the unit square (0,1)^2 with the velocity
    %   beta = (10*y^2-12*x+1, 1+y), whose inflow sides are x = 0, y = 0
    %   and x = 1. The Dirichlet data is 1 on {x = 0, y <= 1/2} and on
    %   {y = 0, x <= 1/2}, sin(pi*y)^2 on x = 1 and 0 elsewhere.
    %
    %   problemTransportSquare2 functions:
    %       f     - Zero source term
    %       u4Db  - Discontinuous Dirichlet inflow data
    %       g     - Zero Neumann data
    %       beta  - Convection velocity (10*y^2-12*x+1, 1+y)
    %
    %   problemTransportSquare2 properties:
    %       geom  - Unit square geometry
    %
    %   See also ProblemScalarEllipticTemplate, problemTransportSquare

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
    end

    methods

        function val = f(~,x)
            val = zeros(size(x,1),1);
        end

        function val = u4Db(~,x)
            y = x(:,2);
            x = x(:,1);
            val = zeros(size(x,1),1);
            val(x==0 & y <= 0.5) = 1;
            val(x<=0.5 & y == 0) = 1;
            val(x==1) = sin(pi*y(x==1)).^2;
        end

        function val = g(~,x,~)
            val = zeros(size(x,1),1);
        end

        function val = beta(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [10*y.^2-12*x+1,1+y];
        end

    end

end