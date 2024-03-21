classdef problemSurfaceSphereExact < problemSurfaceSphere & ProblemScalarEllipticTemplateExact
    %problemSurfaceSphereExact - Laplace-Beltrami on a sphere, zonal solution
    %   problemSurfaceSphereExact is a scalar elliptic (Laplace-Beltrami)
    %   problem on the unit sphere with the axisymmetric (zonal) exact
    %   solution u = (3*sin(theta)^2 - 2)/(8*pi), theta the polar angle.
    %   The data require the radius r = 1; other radii raise the error
    %   fem4m:problemSurfaceSphereExact:radius.
    %
    %   problemSurfaceSphereExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   See also problemSurfaceSphere, ProblemScalarEllipticTemplateExact,
    %   problemSurfaceSphereExact2, problemSurfaceSphereExact3

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = f(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereExact:radius','The data of problemSurfaceSphereExact hold only for the radius r = 1.');
            end
            [~,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            val =  (9*sin(theta).^2 - 6)/(4*pi) ;
        end

        function val = exactU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereExact:radius','The data of problemSurfaceSphereExact hold only for the radius r = 1.');
            end
            [~,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            val = (3*sin(theta).^2-2)/(8*pi);
        end

        function val = exactGradU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereExact:radius','The data of problemSurfaceSphereExact hold only for the radius r = 1.');
            end
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            dtheta = (6*sin(theta).*cos(theta))/(8*pi);
            val = dtheta.*[cos(theta).*cos(phi),cos(theta).*sin(phi),-sin(theta)];
        end

    end

end