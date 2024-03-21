classdef problemSurfaceStokesSphereExact2 < ProblemStokesTemplateExact & ProblemSurfaceTemplate
    %problemSurfaceStokesSphereExact2 - Stokes flow on a sphere, non-solenoidal velocity
    %   problemSurfaceStokesSphereExact2 is a surface Stokes problem on the
    %   sphere of radius r with the exact velocity u = P*(-z^2, x, y) and
    %   the exact pressure p = x*y^3 + z, where P = I - n*n' is the
    %   tangential projection. The body force is
    %   f = -2*P*div_G(E(u)) + P*grad(p) with the surface gradient
    %   grad_G(u) = P*grad(u)*P, E(u) = (grad_G(u) + grad_G(u)')/2 and div_G
    %   taken row-wise, but the velocity is not solenoidal: the mass balance
    %   reads div_G(u) = exactDivU instead of div_G(u) = 0. exactGradU
    %   returns grad_G(u). f, exactU, exactGradU and exactP hold for every
    %   radius r; f assumes the viscosity NU = 1. The sphere is closed, so
    %   the boundary data are not used.
    %
    %   problemSurfaceStokesSphereExact2 functions:
    %       f                   - Body force
    %       exactU              - Exact velocity
    %       exactGradU          - Velocity gradient
    %       exactP              - Exact pressure
    %       surfaceFunction     - Level-set function |x| - r of the sphere
    %       gradSurfaceFunction - Gradient of the level-set function
    %       d2SurfaceFunction   - Hessian of the level-set function
    %       curvature           - Principal curvatures (both 1/r)
    %
    %   problemSurfaceStokesSphereExact2 properties:
    %       geom                - Sphere geometry
    %       r                   - Radius
    %
    %   See also ProblemStokesTemplateExact, problemSurfaceStokesSphereExact,
    %   problemSurfaceSphere

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Sphere'))
        r = 1
    end


    methods

        function val = f(obj,x)
            P = obj.normalProjector(x);
            z = x(:,3);y=x(:,2);x=x(:,1);
            L = [1.0./(x.^2+y.^2+z.^2).^2.*(x.^2.*y.^2.*2.0+x.^2.*z.^2.*2.2e+1-y.^2.*z.^2.*6.0-x.^2.*y.*1.1e+1+y.*z.^2.*5.0+y.^3.*5.0+y.^4.*2.0-z.^4.*8.0-x.*y.*z.*1.6e+1).*(-1.0./2.0),1.0./(x.^2+y.^2+z.^2).^2.*(x.*y.^2.*-1.1e+1-x.*y.^3.*2.0-x.^3.*y.*2.0+x.*z.^2.*5.0+x.^2.*z.*5.0-y.^2.*z.*1.1e+1+x.^3.*5.0+z.^3.*5.0+x.*y.*z.^2.*2.8e+1).*(-1.0./2.0),(1.0./(x.^2+y.^2+z.^2).^2.*(x.^2.*y.*-5.0-x.*z.^3.*1.6e+1+x.^3.*z.*1.4e+1+y.*z.^2.*1.1e+1-y.^3.*5.0+x.*y.^2.*z.*1.4e+1+x.*y.*z.*1.6e+1))./2.0];
            L = matMul(P,L);
            p = [y.^3,3*x.*y.^2,ones(size(x,1),1)];
            p = matMul(P,p);
            val = -2*L + p;
        end

        function val = exactU(obj,x)
            val = [-x(:,3).^2,x(:,1),x(:,2)];
            val = matMul(obj.normalProjector(x),val);
        end

        function val = exactGradU(obj,x)
            P = obj.normalProjector(x);
            z = x(:,3);y=x(:,2);x=x(:,1);
            val = reshape([-(x.*y.*2.0+y.*z)./(x.^2+y.^2+z.^2)+x.*1.0./(x.^2+y.^2+z.^2).^2.*(y.^2.*z.^2+x.^2.*y+z.^4+x.*y.*z).*2.0,(y.*z.^2+x.^2.*3.0+z.^2)./(x.^2+y.^2+z.^2)-x.*1.0./(x.^2+y.^2+z.^2).^2.*(x.*z.^2-y.^2.*z+x.^3+x.*y.*z.^2).*2.0,z.*1.0./(x.^2+y.^2+z.^2).^2.*(-x.^2.*z.^2+y.^2.*z.^2+x.^2.*y-y.*z.^2-y.^3+z.^4+x.*y.*z.*2.0),-x.*1.0./(x.^2+y.^2+z.^2).^2.*(-x.*y.^2+x.*z.^2+x.^2.*z-y.^2.*z+x.^3+z.^3+x.*y.*z.^2.*2.0),-(y.*z.*2.0-x.*z.^2)./(x.^2+y.^2+z.^2)-y.*1.0./(x.^2+y.^2+z.^2).^2.*(x.*z.^2-y.^2.*z+x.^3+x.*y.*z.^2).*2.0,(-x.*z+x.^2+y.^2.*3.0)./(x.^2+y.^2+z.^2)-y.*1.0./(x.^2+y.^2+z.^2).^2.*(x.^2.*y+x.*z.^3+y.^3-x.*y.*z).*2.0,-(x.*y+y.^2.*z.*2.0+z.^3.*4.0)./(x.^2+y.^2+z.^2)+z.*1.0./(x.^2+y.^2+z.^2).^2.*(y.^2.*z.^2+x.^2.*y+z.^4+x.*y.*z).*2.0,y.*1.0./(x.^2+y.^2+z.^2).^2.*(-x.^2.*y+x.^3.*z.*2.0+y.*z.^2-y.^3+x.*y.^2.*z.*2.0+x.*y.*z.*2.0),-(x.*y-x.*z.^2.*3.0)./(x.^2+y.^2+z.^2)-z.*1.0./(x.^2+y.^2+z.^2).^2.*(x.^2.*y+x.*z.^3+y.^3-x.*y.*z).*2.0],[],3,3);
            val = matMul(P,val);
            val = matMul(val,P);
        end

        function val = exactP(~,x)
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

    end

end