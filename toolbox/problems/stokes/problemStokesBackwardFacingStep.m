classdef problemStokesBackwardFacingStep < ProblemStokesTemplate
    %problemStokesBackwardFacingStep - Flow past a backward-facing step
    %   problemStokesBackwardFacingStep is a Stokes flow problem with zero
    %   body force past the backward-facing step (-2,8)x(-1,1) minus
    %   [-2,0]x[-1,0]. The walls are no-slip Dirichlet boundaries, the
    %   inflow side x = -2 and the outflow side x = 8 Neumann boundaries.
    %   u4Db holds the parabolic profiles u1 = y*(1-y)/10 at x = -2 and
    %   u1 = (1-y^2)/80 at x = 8, both with flux 1/60, for solvers that
    %   impose Dirichlet data on all sides; g returns eps(u)*n of these
    %   profiles, without the factor 2*NU and without the pressure.
    %
    %   OBJ = problemStokesBackwardFacingStep() creates the problem with the
    %   viscosity NU = 1.
    %
    %   OBJ = problemStokesBackwardFacingStep(NU) sets the viscosity NU.
    %
    %   problemStokesBackwardFacingStep functions:
    %       f     - Zero body force
    %       u4Db  - Inflow and outflow profiles, zero on the walls
    %       g     - eps(u)*n of the profiles at x = -2 and x = 8
    %
    %   problemStokesBackwardFacingStep properties:
    %       geom  - Backward-facing step geometry
    %
    %   See also ProblemStokesTemplate, problemStokesDrivenCavity

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','BackwardFacingStep'))
    end


    methods

        function val = f(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [0*x,0*y];
        end

        function val = u4Db(~,x)
            val = zeros(size(x));
            y = x(:,2);
            x = x(:,1);
            val(x==-2,1) = -y(x==-2).*(y(x==-2)-1)/10;
            val(x==8,1) = -(y(x==8)+1).*(y(x==8)-1)/80;
        end

        function val = g(~,x,~)
            y = x(:,2);
            x = x(:,1);
            val = zeros(length(x),2);
            val(x==-2,2) = 1/20*(2*y(x==-2)-1);
            val(x== 8,2) = -1/160*2*y(x==8);
        end

    end
end