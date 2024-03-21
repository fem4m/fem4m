classdef problemSurfaceStokesG88 < ProblemStokesTemplate & problemSurfaceG88
    %problemSurfaceStokesG88 - Stokes flow on the G88 surface
    %   problemSurfaceStokesG88 is a surface Stokes problem on the closed
    %   G88 surface (x - z^2)^2 + y^2 + z^2 = 1 of problemSurfaceG88 with the
    %   body force f = (sin(pi*x), sin(pi*y), sin(pi*z)), which is not
    %   projected onto the tangent plane, and zero boundary data, which play
    %   no role on the closed surface. No exact solution is known.
    %
    %   problemSurfaceStokesG88 functions:
    %       f     - Body force sin(pi*x), componentwise
    %       u4Db  - Zero velocity data
    %       g     - Zero traction data
    %
    %   See also ProblemStokesTemplate, problemSurfaceG88,
    %   problemSurfaceStokesSphereExact

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = f(~,x)
            val = sin(pi.*x);
        end

        function val = u4Db(~,x)
            val = zeros(size(x,1),3);
        end

        function val = g(~,x,~)
            val = zeros(size(x,1),3);
        end

    end
end