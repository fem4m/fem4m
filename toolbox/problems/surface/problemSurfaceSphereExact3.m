classdef problemSurfaceSphereExact3 < problemSurfaceSphere & ProblemScalarEllipticTemplateExact
    %problemSurfaceSphereExact3 - Laplace-Beltrami on a sphere, polynomial
    %   problemSurfaceSphereExact3 is a scalar elliptic (Laplace-Beltrami)
    %   problem on the unit sphere with the polynomial exact solution
    %   u = x*y^3 + z.
    %   The data require the radius r = 1; other radii raise the error
    %   fem4m:problemSurfaceSphereExact3:radius.
    %
    %   problemSurfaceSphereExact3 functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   See also problemSurfaceSphere, problemSurfaceSphereExact,
    %   problemSurfaceSphereExact2

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = f(obj,x)
            %F - Source term
            % delta(u) = tr(P*H)-kappa*grad(u)*n
            %          = tr(P*H)-2*grad(u)*x/||x||
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereExact3:radius','The data of problemSurfaceSphereExact3 hold only for the radius r = 1.');
            end
            P = obj.normalProjector(x);
            n = x./sqrt(sum(x.^2,2));
            H = zeros(size(x,1),3,3);
            H(:,1,2) = 3*x(:,2).^2;
            H(:,2,1) = 3*x(:,2).^2;
            H(:,2,2) = 6*x(:,1).*x(:,2);
            grad = [x(:,2).^3,3*x(:,1).*x(:,2).^2,ones(size(x,1),1)];
            L = matMul(P,H);
            L = L(:,1,1) + L(:,2,2) + L(:,3,3) - 2*sum(grad.*n,2);
            val = -L;
        end

        function val = exactU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereExact3:radius','The data of problemSurfaceSphereExact3 hold only for the radius r = 1.');
            end
            val = x(:,1).*x(:,2).^3+x(:,3);
        end

        function val = exactGradU(obj,x)
            if obj.r ~= 1
                error('fem4m:problemSurfaceSphereExact3:radius','The data of problemSurfaceSphereExact3 hold only for the radius r = 1.');
            end
            val = [x(:,2).^3,3*x(:,1).*x(:,2).^2,ones(size(x,1),1)];
            val = matMul(obj.normalProjector(x),val);
        end

    end

end