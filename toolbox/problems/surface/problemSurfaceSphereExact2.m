classdef problemSurfaceSphereExact2 < problemSurfaceSphere & ProblemScalarEllipticTemplateExact
    %problemSurfaceSphereExact2 - Laplace-Beltrami on a sphere, oscillatory
    %   problemSurfaceSphereExact2 is a scalar elliptic (Laplace-Beltrami)
    %   problem on the unit sphere with the oscillatory exact solution
    %   u = cos(5*theta + 3*phi)*sin(theta)^2 in spherical coordinates
    %   (theta polar angle, phi azimuth).
    %   The data require the radius r = 1; other radii raise the error
    %   fem4m:problemSurfaceSphereExact2:radius.
    %
    %   problemSurfaceSphereExact2 functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   See also problemSurfaceSphere, problemSurfaceSphereExact,
    %   problemSurfaceSphereExact3

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = f(obj,x)
            %F - Source term
            % delta(u) = 1/(sin(theta))*d/dtheta(sin(theta) du(x)/dtheta)
            %          + 1/(sin(theta)^2)*d^2u(x)/dphi^2
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereExact2:radius','The data of problemSurfaceSphereExact2 hold only for the radius r = 1.');
            end
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            val =  -((3*cos(3*phi + 3*theta))/2 ...
                - (41*cos(3*phi + 5*theta))/2 + 14*cos(3*phi + 7*theta)) ;
        end

        function val = exactU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereExact2:radius','The data of problemSurfaceSphereExact2 hold only for the radius r = 1.');
            end
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            val = cos(5*theta+3*phi).*sin(theta).^2;
        end

        function val = exactGradU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereExact2:radius','The data of problemSurfaceSphereExact2 hold only for the radius r = 1.');
            end
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            dtheta = -(sin(theta).*(3*cos(3*phi + 4*theta) - 7*cos(3*phi + 6*theta)))/2;
            dphi = -3*sin(3*phi + 5*theta).*sin(theta).^2;
            val = dtheta.*[cos(theta).*cos(phi),cos(theta).*sin(phi),-sin(theta)] ...
                + dphi./sin(theta).*[-sin(phi),cos(phi),zeros(numel(phi),1)];
        end

    end

end