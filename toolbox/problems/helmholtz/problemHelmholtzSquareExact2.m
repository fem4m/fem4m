classdef problemHelmholtzSquareExact2 < ProblemHelmholtzTemplateExact
    %problemHelmholtzSquareExact2 - Helmholtz with a Bessel solution
    %   problemHelmholtzSquareExact2 is a Helmholtz problem on the shifted
    %   square (0,1)x(-1/2,1/2) with the exact solution
    %   u = besselj(ALPHA,OMEGA*r)*cos(ALPHA*phi) in polar coordinates about
    %   the origin, the midpoint of the left edge, phi in [-pi/2,pi/2]. For
    %   non-integer ALPHA (e.g. ALPHA = 2/3) the solution is singular at
    %   this boundary point.
    %
    %   OBJ = problemHelmholtzSquareExact2(OMEGA,ALPHA) sets the angular
    %   frequency OMEGA and the Bessel order ALPHA (default 1). OMEGA is
    %   required.
    %
    %   problemHelmholtzSquareExact2 functions:
    %       exactU      - Exact Bessel solution
    %       exactGradU  - Gradient of the exact solution
    %
    %   problemHelmholtzSquareExact2 properties:
    %       geom        - Shifted square geometry
    %       alpha       - Bessel order (angular exponent)
    %
    %   See also ProblemHelmholtzTemplateExact, problemHelmholtzLshapeExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','midpoint',[0.5,0]))
        alpha
    end

    methods

        function obj = problemHelmholtzSquareExact2(omega,alpha)
            if nargin < 1
                omega = [];
            end
            obj = obj@ProblemHelmholtzTemplateExact(omega);
            if nargin < 2 || isempty(alpha)
                alpha = 1;
            end
            obj.alpha = alpha;
        end

        function val = exactU(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.omega;
            [th,r]=cart2pol(x,y);
            xi = obj.alpha;
            val = besselj(xi,k*r).*cos(xi*th);
        end

        function val = exactGradU(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.omega;
            [th,r]=cart2pol(x,y);
            xi = obj.alpha;
            val = [k*0.5*(besselj(xi-1,k*r)-besselj(xi+1,k*r)).*cos(xi*th).*(x./r)+xi*besselj(xi,k*r).*sin(xi*th).*(y./r.^2),...
                k*0.5*(besselj(xi-1,k*r)-besselj(xi+1,k*r)).*cos(xi*th).*(y./r)-xi*besselj(xi,k*r).*sin(xi*th).*(x./r.^2)];
        end

    end
end