classdef problemStokesExactCollidingFlow < ProblemStokesTemplateExact
    %problemStokesExactCollidingFlow - Colliding-flow benchmark
    %   problemStokesExactCollidingFlow is a Stokes benchmark on the square
    %   (-1,1)^2 with the known velocity u = (20xy^4 - 4x^5, 20x^4y - 4y^5)
    %   and the zero-mean pressure p = NU*(120x^2y^2 - 20x^4 - 20y^4 - 16/3)
    %   modelling two colliding streams. Pressure and stress scale with the
    %   viscosity NU, so the body force vanishes for every NU.
    %
    %   OBJ = problemStokesExactCollidingFlow() creates the problem with
    %   the viscosity NU = 1.
    %
    %   OBJ = problemStokesExactCollidingFlow(NU) sets the viscosity NU.
    %
    %   problemStokesExactCollidingFlow functions:
    %       f                - Zero body force
    %       exactU           - Exact velocity
    %       exactGradU       - Velocity gradient
    %       exactP           - Exact pressure
    %       exactSigma       - Exact stress 2*NU*eps(u) - p*I as
    %                          [s11 s22 s12]
    %       exactRotEpsilon  - Row-wise rot of the strain
    %
    %   problemStokesExactCollidingFlow properties:
    %       geom             - Square geometry
    %
    %   See also ProblemStokesTemplateExact, problemStokesExactSquareP3Sigma

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','midpoint',[0,0],'length',2))
    end


    methods

        function val = f(~,x)
            val = zeros(size(x));
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [20.*x.*y.^4 - 4.*x.^5,20.*x.^4.*y - 4.*y.^5];
        end

        function val = exactGradU(~,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1,1) = 20.*y.^4 - 20.*x.^4;
            val(:,1,2) = 80.*x.*y.^3;
            val(:,2,1) = 80.*x.^3.*y;
            val(:,2,2) = 20.*x.^4 - 20.*y.^4;
        end

        function val = exactP(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = obj.nu*(120.*x.^2.*y.^2 - 20.*x.^4 - 20.*y.^4 - 32/6);
        end

        function val = exactSigma(obj,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = 60.*y.^4 - 120.*x.^2.*y.^2 - 20.*x.^4 + 32/6;
            val(:,2) = 60.*x.^4 - 120.*x.^2.*y.^2 - 20.*y.^4 + 32/6;
            val(:,3) = 80.0.*x.^3.*y + 80.0.*x.*y.^3;
            val = obj.nu*val;
        end

        function val = exactRotEpsilon(~,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = 120.*x.^2.*y - 40.*y.^3;
            val(:,2) = 40*x.^3 - 120*x.*y.^2;
        end
    end
end