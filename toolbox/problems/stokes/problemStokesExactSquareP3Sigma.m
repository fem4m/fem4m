classdef problemStokesExactSquareP3Sigma < ProblemStokesTemplateExact
    %problemStokesExactSquareP3Sigma - Stokes problem with a P3 stress
    %   problemStokesExactSquareP3Sigma is a Stokes problem on the square
    %   (-1,1)^2 with a cubic (P3) exact stress: a quartic exact velocity
    %   u = (20xy^3, 5x^4 - 5y^4), the cubic exact pressure
    %   p = NU*(60x^2y - 20y^3) for the viscosity NU, and zero body force.
    %
    %   OBJ = problemStokesExactSquareP3Sigma() creates the problem with
    %   the viscosity NU = 1.
    %
    %   OBJ = problemStokesExactSquareP3Sigma(NU) sets the viscosity NU.
    %
    %   problemStokesExactSquareP3Sigma functions:
    %       f                - Zero body force
    %       exactU           - Exact velocity
    %       exactGradU       - Velocity gradient
    %       exactP           - Exact pressure
    %       exactSigma       - Exact stress 2*NU*eps(u) - p*I as
    %                          [s11 s22 s12]
    %       exactRotEpsilon  - Row-wise rot of the strain
    %
    %   problemStokesExactSquareP3Sigma properties:
    %       geom             - Square geometry
    %
    %   See also ProblemStokesTemplateExact, problemStokesExactCollidingFlow,
    %   problemStokesSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','midpoint',[0,0],'length',2))
    end


    methods


        function val = f(~,x)
            %F - Zero body force
            y = x(:,2);
            x = x(:,1);
            %     val(:,1) = 60*(-2*x.*y);
            %     val(:,2) = 60*(y.^2-x.^2);
            val(:,1) = 0*x;
            val(:,2) = 0*y;
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [20*x.*y.^3,5*x.^4 - 5*y.^4];
        end

        function val = exactGradU(~,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1,1) = 20*y.^3;
            val(:,1,2) = 60*x.*y.^2;
            val(:,2,1) = 20*x.^3;
            val(:,2,2) = - 20*y.^3;
        end

        function val = exactP(obj,x)
            y = x(:,2);
            x = x(:,1);
            %     val = 0*(60*x.^2.*y - 20*y.^3);
            val = obj.nu*(60*x.^2.*y - 20*y.^3);
        end

        function val = exactSigma(obj,x)
            y = x(:,2);
            x = x(:,1);
            %     val = [40.*y.^3, -40.*y.^3,  60.*x.*y.^2+20.*x.^3];
            val = obj.nu*[60.*y.^3-60.*x.^2.*y, -20.*y.^3-60.*x.^2.*y,  60.*x.*y.^2+20.*x.^3];
        end

        function val = exactRotEpsilon(~,x)
            y = x(:,2);
            x = x(:,1);
            val(:,1) = 30.*x.^2 - 30.*y.^2;
            val(:,2) = (-60).*x.*y;
        end

    end

end