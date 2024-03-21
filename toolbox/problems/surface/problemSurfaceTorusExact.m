classdef problemSurfaceTorusExact < problemSurfaceTorus & ProblemScalarEllipticTemplateExact
    %problemSurfaceTorusExact - Laplace-Beltrami on a torus
    %   problemSurfaceTorusExact is a scalar elliptic (Laplace-Beltrami)
    %   problem on the torus of problemSurfaceTorus with the exact
    %   solution u = exp(1/(1.85 - x^2))*sin(y).
    %
    %   problemSurfaceTorusExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   See also problemSurfaceTorus, ProblemScalarEllipticTemplateExact
    %
    %   Note: Based on Example 2 of Demlow and Dziuk, An adaptive finite
    %   element method for the Laplace-Beltrami operator on implicitly
    %   defined surfaces, 2007.

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = f(obj,x)
            Uxx = (800*exp(-1./(x(:,1).^2 - 37/20)).*sin(x(:,2)).*(- 1200*x(:,1).^4 + 2280*x(:,1).^2 + 1369))./(20*x(:,1).^2 - 37).^4;
            Uyy = -exp(-1./(x(:,1).^2 - 37/20)).*sin(x(:,2));
            Uxy = (2*x(:,1).*exp(-1./(x(:,1).^2 - 37/20)).*cos(x(:,2)))./(x(:,1).^2 - 37/20).^2;
            % delta(u) = tr(P*H)-kappa*grad(u)*n
            kappa = sum(obj.curvature(x),2);
            n = obj.normal(x);
            P = obj.normalProjector(x);
            H = zeros(size(x,1),3,3);
            H(:,1,2) = Uxy;
            H(:,2,1) = Uxy;
            H(:,1,1) = Uxx;
            H(:,2,2) = Uyy;
            grad = [(2.*x(:,1).*exp(-1./(x(:,1).^2 - 37/20)).*sin(x(:,2)))./(x(:,1).^2 - 37/20).^2, exp(-1./(x(:,1).^2 - 37/20)).*cos(x(:,2)), zeros(numel(x(:,1)),1)];
            L = matMul(P,H);
            L = L(:,1,1) + L(:,2,2) + L(:,3,3) - kappa.*sum(grad.*n,2);
            val = -L;

        end

        function val = exactU(~,x)
            val = exp(1./(1.85-x(:,1).^2)).*sin(x(:,2));
        end

        function val = exactGradU(obj,x)
            val = [(2.*x(:,1).*exp(-1./(x(:,1).^2 - 37/20)).*sin(x(:,2)))./(x(:,1).^2 - 37/20).^2, exp(-1./(x(:,1).^2 - 37/20)).*cos(x(:,2)), zeros(numel(x(:,1)),1)];
            val = matMul(obj.normalProjector(x),val);
        end

    end

end