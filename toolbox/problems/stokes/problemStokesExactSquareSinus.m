classdef problemStokesExactSquareSinus < ProblemStokesTemplateExact
    %problemStokesExactSquareSinus - Stokes problem with a trigonometric flow
    %   problemStokesExactSquareSinus is a Stokes problem on the unit square
    %   (0,1)^2 with the trigonometric exact velocity u = (s,s),
    %   s = sin(pi*(x-y)), and the zero-mean pressure
    %   p = -2*pi*cos(pi*(x-y)) + 8/pi, which does not depend on NU; the
    %   body force is f = -NU*Delta u + grad p for the viscosity NU.
    %
    %   OBJ = problemStokesExactSquareSinus() creates the problem with
    %   the viscosity NU = 1.
    %
    %   OBJ = problemStokesExactSquareSinus(NU) sets the viscosity NU.
    %
    %   problemStokesExactSquareSinus functions:
    %       f                - Body force
    %       exactU           - Exact velocity
    %       exactGradU       - Velocity gradient
    %       exactP           - Exact pressure
    %       exactSigma       - Exact stress 2*NU*eps(u) - p*I as
    %                          [s11 s22 s12]
    %       exactRotEpsilon  - Row-wise rot of the strain
    %
    %   problemStokesExactSquareSinus properties:
    %       geom             - Square geometry
    %
    %   See also ProblemStokesTemplateExact, problemStokesSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
    end


    methods

        function val = f(obj,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = 2*(obj.nu+1)*pi^2*sin(pi*(x - y));
            val(:,2) = 2*(obj.nu-1)*pi^2*sin(pi*(x - y));
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = sin(pi*x).*cos(pi*y) - cos(pi*x).*sin(pi*y);
            val(:,2) = val(:,1);
        end

        function val = exactGradU(~,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1,1) = pi*cos(pi*x).*cos(pi*y) + pi*sin(pi*x).*sin(pi*y);
            val(:,1,2) = -pi*sin(pi*x).*sin(pi*y) - pi*cos(pi*x).*cos(pi*y);
            val(:,2,1) = val(:,1,1);
            val(:,2,2) = val(:,1,2);
        end

        function val = exactP(~,x)
            y = x(:,2);
            x = x(:,1);
            val = (-2)*pi*(cos(pi*x).*cos(pi*y) + sin(pi*x).*sin(pi*y))+8/pi;
        end

        function val = exactSigma(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = cos(pi*x).*cos(pi*y) + sin(pi*x).*sin(pi*y);
            p = obj.exactP([x,y]);
            val = [2*pi*obj.nu*val - p, -2*pi*obj.nu*val - p, zeros(size(x,1),1)];
        end

        function val = exactRotEpsilon(~,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = pi^2*cos(pi*x).*sin(pi*y) - pi^2*cos(pi*y).*sin(pi*x);
            val(:,2) = pi^2*cos(pi*y).*sin(pi*x) - pi^2*cos(pi*x).*sin(pi*y);
        end
    end
end