classdef problemHelmholtzLshapeExact < ProblemHelmholtzTemplateExact
    %problemHelmholtzLshapeExact - Helmholtz on an L-shape
    %   problemHelmholtzLshapeExact is a Helmholtz problem on the L-shape
    %   (-1,1)^2 \ [0,1]x[-1,0] with the corner-expansion exact solution
    %   u = besselj(THETA,OMEGA*r)*sin(THETA*phi) in polar coordinates about
    %   the re-entrant corner (the origin), phi in [0,3*pi/2], which is
    %   singular at the re-entrant corner for non-integer THETA (e.g.
    %   THETA = 2/3).
    %
    %   OBJ = problemHelmholtzLshapeExact(OMEGA,THETA) sets the angular
    %   frequency OMEGA and the corner exponent THETA (default THETA = 1,
    %   which gives a smooth solution). OMEGA is required.
    %
    %   problemHelmholtzLshapeExact functions:
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemHelmholtzLshapeExact properties:
    %       geom        - L-shaped domain geometry
    %       theta       - Corner exponent (Bessel order)
    %
    %   See also ProblemHelmholtzTemplateExact,
    %   problemHelmholtzLshapeSmoothExact, problemHelmholtzSquareExact2

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape'))
        theta
    end

    methods

        function obj = problemHelmholtzLshapeExact(omega,theta)
            if nargin < 1
                omega = [];
            end
            obj = obj@ProblemHelmholtzTemplateExact(omega);
            if nargin < 2 || isempty(theta)
                theta = 1;
            end
            obj.theta = theta;
        end

        function val = exactU(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.omega;
            [th,r]=cart2pol(x,y);
            th( th<0 ) = th( th<0 ) + 2*pi;
            xi = obj.theta;
            val = besselj(xi,k*r).*sin(xi*th);
        end

        function val = exactGradU(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.omega;
            [th,r]=cart2pol(x,y);
            th( th<0 ) = th( th<0 ) + 2*pi;
            xi = obj.theta;
            val = [k*0.5*(besselj(xi-1,k*r)-besselj(xi+1,k*r)).*sin(xi*th).*(x./r)-xi*besselj(xi,k*r).*cos(xi*th).*(y./r.^2),...
                k*0.5*(besselj(xi-1,k*r)-besselj(xi+1,k*r)).*sin(xi*th).*(y./r)+xi*besselj(xi,k*r).*cos(xi*th).*(x./r.^2)];
        end

    end
end