classdef problemSurfaceEllipsoidExact < problemSurfaceEllipsoid & ProblemScalarEllipticTemplateExact
    %problemSurfaceEllipsoidExact - Laplace-Beltrami on an ellipsoid
    %   OBJ = problemSurfaceEllipsoidExact is a scalar elliptic
    %   (Laplace-Beltrami) problem on the ellipsoid x^2 + y^2 + z^2/S = R
    %   of problemSurfaceEllipsoid with R = 1 and S = 400 (semi-axes 1, 1
    %   and 20) and the exact solution u = sin(y).
    %
    %   OBJ = problemSurfaceEllipsoidExact(R) sets the level value R.
    %
    %   OBJ = problemSurfaceEllipsoidExact(R,S) also sets the scaling S of
    %   z^2; the semi-axes are sqrt(R), sqrt(R) and sqrt(R*S). An empty R
    %   or S selects its default.
    %
    %   problemSurfaceEllipsoidExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   See also problemSurfaceEllipsoid, ProblemScalarEllipticTemplateExact
    %
    %   Note: Based on Example 3 of Demlow and Dziuk, An adaptive finite
    %   element method for the Laplace-Beltrami operator on implicitly
    %   defined surfaces, 2007.

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = problemSurfaceEllipsoidExact(r,s)
            if nargin < 1 || isempty(r)
                r = 1;
            end
            if nargin < 2 || isempty(s)
                s = 400;
            end
            obj = obj@problemSurfaceEllipsoid(r,s);
        end

        function val = f(obj,x)
            %F - Source term
            % delta(u) = tr(P*H)-kappa*grad(u)*n
            kappa = sum(obj.curvature(x),2);
            n = obj.normal(x);
            P = obj.normalProjector(x);
            H = zeros(size(x,1),3,3);
            H(:,2,2) = -sin(x(:,2));
            grad = [zeros(numel(x(:,1)),1), cos(x(:,2)), zeros(numel(x(:,1)),1)];
            L = matMul(P,H);
            L = L(:,1,1) + L(:,2,2) + L(:,3,3) - kappa.*sum(grad.*n,2);
            val = -L;
        end

        function val = exactU(~,x)
            val = sin(x(:,2));
        end

        function val = exactGradU(obj,x)
            val = [zeros(numel(x(:,1)),1), cos(x(:,2)), zeros(numel(x(:,1)),1)];
            val = matMul(obj.normalProjector(x),val);
        end

    end

end