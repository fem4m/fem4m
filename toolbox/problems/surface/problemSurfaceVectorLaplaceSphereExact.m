classdef problemSurfaceVectorLaplaceSphereExact < ProblemSurfaceTemplate
    %problemSurfaceVectorLaplaceSphereExact - Vector Laplace on a sphere, axisymmetric solution
    %   problemSurfaceVectorLaplaceSphereExact is the surface vector Laplace
    %   problem -P*div_G(grad_G(u)) = f on the unit sphere, with the
    %   tangential projection P = I - n*n', the surface gradient
    %   grad_G(u) = P*grad(u)*P and div_G taken row-wise. The exact solution
    %   is the tangential field u = sin(theta)^2*e_theta for the polar angle
    %   theta and the unit vector e_theta of increasing theta; exactGradU
    %   returns grad_G(u). S returns the basis [e_r,e_theta,e_phi] at the
    %   azimuth phi and the polar angle theta.
    %   The data require the radius r = 1; other radii raise the error
    %   fem4m:problemSurfaceVectorLaplaceSphereExact:radius.
    %
    %   f, exactU and exactGradU depend on the angles only, so the data hold
    %   only for the radius r = 1 (the default). The sphere is closed, so
    %   the problem has no boundary data.
    %
    %   problemSurfaceVectorLaplaceSphereExact functions:
    %       f                   - Body force
    %       exactU              - Exact solution
    %       exactGradU          - Gradient of the solution
    %       surfaceFunction     - Level-set function |x| - r of the sphere
    %       gradSurfaceFunction - Gradient of the level-set function
    %       d2SurfaceFunction   - Hessian of the level-set function
    %       curvature           - Principal curvatures (both 1/r)
    %       S                   - Basis [e_r,e_theta,e_phi] at (phi,theta)
    %
    %   problemSurfaceVectorLaplaceSphereExact properties:
    %       geom                - Sphere geometry
    %       r                   - Radius
    %
    %   See also ProblemSurfaceTemplate, problemSurfaceSphere,
    %   problemSurfaceVectorLaplaceSphereExact2

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Sphere'))
        r = 1
    end


    methods

        function val = f(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceVectorLaplaceSphereExact:radius','The data of problemSurfaceVectorLaplaceSphereExact hold only for the radius r = 1.');
            end
            P = obj.normalProjector(x);
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            val = [cos(phi).*cos(theta).*(cos(theta).^2.*8.0-5.0),cos(theta).*sin(phi).*(cos(theta).^2.*8.0-5.0),sin(theta.*3.0).*-2.0];
            val = -matMul(P,val);
        end

        function val = exactU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceVectorLaplaceSphereExact:radius','The data of problemSurfaceVectorLaplaceSphereExact hold only for the radius r = 1.');
            end
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            null = zeros(size(x,1),1);
            val = [null,sin(theta).^2,null];
            val = matMul(obj.S(phi,theta),val);
        end

        function val = exactGradU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceVectorLaplaceSphereExact:radius','The data of problemSurfaceVectorLaplaceSphereExact hold only for the radius r = 1.');
            end
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            null = zeros(size(x,1),1);
            val = reshape([null,null,null,cos(phi).*sin(theta).*(cos(theta).^2.*3.0-1.0),-sin(phi).*sin(theta).*(sin(theta).^2.*3.0-2.0),cos(theta).*(cos(theta).^2-1.0).*3.0,-cos(theta).*sin(phi).*sin(theta),cos(phi).*cos(theta).*sin(theta),null],[],3,3);
            val = matMul(val,permute(obj.S(phi,theta),[1 3 2]));
            P = obj.normalProjector(x);
            val = matMul(P,val);
        end

        function val = surfaceFunction(obj,x)
            val = sqrt(sum(x.^2,2))-obj.r;
        end

        function val = gradSurfaceFunction(~,x)
            val = [x(:,1)./sqrt(sum(x.^2,2)),x(:,2)./sqrt(sum(x.^2,2)),x(:,3)./sqrt(sum(x.^2,2))];
        end

        function val = d2SurfaceFunction(~,x)
            R = sqrt(sum(x.^2,2));
            val(:,1,1) = 1./R - x(:,1).^2./R.^3;
            val(:,2,2) = 1./R - x(:,2).^2./R.^3;
            val(:,3,3) = 1./R - x(:,3).^2./R.^3;
            val(:,1,2) = -x(:,1).*x(:,2)./R.^3;
            val(:,1,3) = -x(:,1).*x(:,3)./R.^3;
            val(:,2,3) = -x(:,2).*x(:,3)./R.^3;
            val(:,2,1) = val(:,1,2);
            val(:,3,1) = val(:,1,3);
            val(:,3,2) = val(:,2,3);
        end

        function val = curvature(obj,x)
            val = ones(size(x,1),2)./obj.r;
        end

        function S = S(~,phi,theta)
            S = reshape([cos(phi).*sin(theta),sin(phi).*sin(theta),cos(theta),cos(phi).*cos(theta),cos(theta).*sin(phi),-sin(theta),-sin(phi),cos(phi),zeros(size(phi,1),1)],[],3,3);
        end

    end

end