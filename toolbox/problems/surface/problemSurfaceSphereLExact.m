classdef problemSurfaceSphereLExact < problemSurfaceSphereL & ProblemScalarEllipticTemplateExact
    %problemSurfaceSphereLExact - Laplace-Beltrami on a three-quarter sphere
    %   problemSurfaceSphereLExact is a scalar elliptic (Laplace-Beltrami)
    %   problem on the three-quarter unit sphere of problemSurfaceSphereL
    %   with the singular exact solution u = sin(theta)^alpha*sin(alpha*phi)
    %   (theta polar angle, phi in [0,3*pi/2] the azimuth), alpha = 2/3.
    %   It vanishes on the boundary, the source term is
    %   f = alpha*(alpha+1)*u, and the gradient is singular at the poles,
    %   the re-entrant corners.
    %   The data require the radius r = 1; other radii raise the error
    %   fem4m:problemSurfaceSphereLExact:radius.
    %
    %   problemSurfaceSphereLExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemSurfaceSphereLExact properties:
    %       alpha       - Singularity exponent (default 2/3)
    %
    %   See also problemSurfaceSphereL, ProblemScalarEllipticTemplateExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        alpha = 2/3
    end

    methods

        function val = f(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereLExact:radius','The data of problemSurfaceSphereLExact hold only for the radius r = 1.');
            end
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            val = obj.alpha.*(obj.alpha+1).*sin(theta).^obj.alpha.*sin(obj.alpha*phi);
        end

        function val = exactU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereLExact:radius','The data of problemSurfaceSphereLExact hold only for the radius r = 1.');
            end
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            val = sin(theta).^obj.alpha.*sin(obj.alpha*phi);
        end

        function val = exactGradU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereLExact:radius','The data of problemSurfaceSphereLExact hold only for the radius r = 1.');
            end
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            dphi = obj.alpha.*cos(obj.alpha*phi).*sin(theta).^obj.alpha;
            dtheta = obj.alpha.*sin(obj.alpha*phi).*cos(theta).*sin(theta).^(obj.alpha - 1);
            val = dtheta.*[cos(theta).*cos(phi),cos(theta).*sin(phi),-sin(theta)] ...
                + dphi./sin(theta).*[-sin(phi),cos(phi),zeros(numel(phi),1)];
        end

    end

end