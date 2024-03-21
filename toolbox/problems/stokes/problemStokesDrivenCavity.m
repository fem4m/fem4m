classdef problemStokesDrivenCavity < ProblemStokesTemplate
    %problemStokesDrivenCavity - Lid-driven cavity flow
    %   problemStokesDrivenCavity is a Stokes flow problem in the cavity
    %   (-1,1)^2 with the lid velocity u = (1,0) on y = 1, no-slip walls
    %   elsewhere and zero body force.
    %
    %   OBJ = problemStokesDrivenCavity() creates the problem with the
    %   viscosity NU = 1.
    %
    %   OBJ = problemStokesDrivenCavity(NU) sets the viscosity NU.
    %
    %   problemStokesDrivenCavity functions:
    %       f     - Zero body force
    %       u4Db  - Lid velocity (1,0) on y = 1, zero elsewhere
    %       g     - Zero Neumann data
    %
    %   problemStokesDrivenCavity properties:
    %       geom  - Square cavity geometry
    %
    %   See also ProblemStokesTemplate, problemStokesBackwardFacingStep

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','midpoint',[0,0],'length',2))
    end


    methods

        function val = f(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [0*x,0*y];
        end

        function val = u4Db(~,x)
            val = zeros(size(x));
            val(x(:,2)==1,1) = 1;
        end

        function val = g(~,x)
            val = zeros(size(x,1),2);
        end

    end
end