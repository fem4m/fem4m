classdef problemSurfaceStokesSphereExact < ProblemStokesTemplateExact & ProblemSurfaceTemplate
    %problemSurfaceStokesSphereExact - Stokes flow on a sphere, solenoidal velocity
    %   problemSurfaceStokesSphereExact is the surface Stokes problem
    %   -2*P*div_G(E(u)) + P*grad(p) = f, div_G(u) = 0 on the unit sphere,
    %   with the tangential projection P = I - n*n', the surface gradient
    %   grad_G(u) = P*grad(u)*P, E(u) = (grad_G(u) + grad_G(u)')/2 and div_G
    %   taken row-wise. The exact solution is the tangential, solenoidal
    %   velocity u = a*e_theta + b*e_phi with a = sin(phi)*sin(theta)^2 and
    %   b = 3*cos(phi)*cos(theta)*sin(theta)^2 for the azimuth phi, the
    %   polar angle theta and the unit vectors e_theta, e_phi of increasing
    %   angles, and the pressure p = x*y^3 + z; exactGradU returns
    %   grad_G(u). S returns the basis [e_r,e_theta,e_phi] at (phi,theta).
    %   The data require the radius r = 1; other radii raise the error
    %   fem4m:problemSurfaceStokesSphereExact:radius.
    %
    %   exactU, exactGradU and the viscous part of f depend on the angles
    %   only, so the data hold only for the radius r = 1 (the default); f
    %   assumes the viscosity NU = 1. The sphere is closed, so the boundary
    %   data are not used.
    %
    %   problemSurfaceStokesSphereExact functions:
    %       f                   - Body force
    %       exactU              - Exact velocity
    %       exactGradU          - Velocity gradient
    %       exactP              - Exact pressure
    %       surfaceFunction     - Level-set function |x| - r of the sphere
    %       gradSurfaceFunction - Gradient of the level-set function
    %       d2SurfaceFunction   - Hessian of the level-set function
    %       curvature           - Principal curvatures (both 1/r)
    %       S                   - Basis [e_r,e_theta,e_phi] at (phi,theta)
    %
    %   problemSurfaceStokesSphereExact properties:
    %       geom                - Sphere geometry
    %       r                   - Radius
    %
    %   See also ProblemStokesTemplateExact, problemSurfaceSphere,
    %   problemSurfaceStokesSphereExact2, problemSurfaceStokesSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Sphere'));
        r = 1;
    end


    methods

        function val = f(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceStokesSphereExact:radius','The data of problemSurfaceStokesSphereExact hold only for the radius r = 1.');
            end
            P = obj.normalProjector(x);
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            L = [cos(phi).*cos(theta).*sin(phi).*(cos(theta).^2-1.0).*-1.0e+1,-cos(theta).*(cos(phi).^2.*1.0e+1-cos(theta).^2.*5.0-cos(phi).^2.*cos(theta).^2.*1.0e+1+1.0),sin(phi).*sin(theta).*(sin(theta).^2.*5.0-4.0)];
            L = matMul(P,L);
            y = x(:,2);x = x(:,1);
            p = [y.^3,3*x.*y.^2,ones(size(x,1),1)];
            p = matMul(P,p);
            val = -2*L + p;
        end

        function val = exactU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceStokesSphereExact:radius','The data of problemSurfaceStokesSphereExact hold only for the radius r = 1.');
            end
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            null = zeros(size(x,1),1);
            val = [null,sin(phi).*sin(theta).^2,cos(phi).*cos(theta).*sin(theta).^2.*3.0];
            val = matMul(obj.S(phi,theta),val);
        end

        function val = exactGradU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceStokesSphereExact:radius','The data of problemSurfaceStokesSphereExact hold only for the radius r = 1.');
            end
            [phi,theta] = cart2sph(x(:,1),x(:,2),x(:,3));
            theta = pi/2 - theta;
            null = zeros(size(x,1),1);
            val = reshape([null,null,null,cos(phi).*sin(phi).*sin(theta).*(cos(theta).^2.*3.0-1.0).*-2.0,sin(theta).*(sin(phi).^2.*2.0-3.0).*(sin(theta).^2.*3.0-2.0),cos(theta).*sin(phi).*(cos(theta).^2-1.0).*3.0,-cos(phi.*2.0).*sin(theta.*2.0),cos(phi).*cos(theta).*sin(phi).*sin(theta).*-4.0,-cos(phi).*sin(theta).^2],[],3,3);
            val = matMul(val,permute(obj.S(phi,theta),[1 3 2]));
            P = obj.normalProjector(x);
            val = matMul(P,val);
        end

        function val = exactP(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceStokesSphereExact:radius','The data of problemSurfaceStokesSphereExact hold only for the radius r = 1.');
            end
            val = x(:,1).*x(:,2).^3+x(:,3);
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