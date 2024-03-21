classdef problemHelmholtzGaussianBeam < ProblemHelmholtzTemplateExact
    %problemHelmholtzGaussianBeam - Helmholtz with a Gaussian beam
    %   problemHelmholtzGaussianBeam is a Helmholtz problem whose reference
    %   solution is the paraxial Gaussian beam of waist W0 through the
    %   origin in the direction (cos(THETA),sin(THETA)). The beam is only
    %   an approximate solution of the Helmholtz equation, not an exact
    %   one: with f = 0, the residual -Laplace(u) - omega^2*u is small but
    %   does not vanish, so the errors of discrete solutions measured
    %   against exactU do not converge to zero.
    %
    %   OBJ = problemHelmholtzGaussianBeam(OMEGA,THETA,W0) sets the angular
    %   frequency OMEGA, the beam angle THETA and the waist W0 on the
    %   square (0,4)^2. OMEGA is required; THETA is in radians. The
    %   defaults are THETA = deg2rad(40) and W0 = 8*pi/OMEGA.
    %
    %   problemHelmholtzGaussianBeam functions:
    %       exactU      - Gaussian beam (approximate solution)
    %       exactGradU  - Gradient of the Gaussian beam
    %
    %   problemHelmholtzGaussianBeam properties:
    %       geom        - Square geometry
    %       theta       - Beam propagation angle (radians)
    %       w0          - Beam waist
    %
    %   See also ProblemHelmholtzTemplateExact, problemHelmholtzPlaneWave

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','midpoint',[2,2],'length',4))
        theta
        w0
    end

    methods

        function obj = problemHelmholtzGaussianBeam(omega,theta,w0)
            if nargin < 1
                omega = [];
            end
            obj = obj@ProblemHelmholtzTemplateExact(omega);
            if nargin < 2 || isempty(theta)
                theta = deg2rad(40);
            end
            obj.theta = theta;
            if nargin < 3 || isempty(w0)
                w0 = 4*2*pi/omega;
            end
            obj.w0 = w0;
        end

        function val = exactU(obj,x)
            val = obj.beam(x);
        end

        function val = exactGradU(obj,x)
            [~,val] = obj.beam(x);
        end

    end

    methods (Access = private)

        function [val,grad] = beam(obj,x)
            %BEAM - Gaussian beam and its gradient
            %   u = w0/w(z)*exp(-rho^2/w(z)^2 - i*k*z - i*k*rho^2/(2*R(z))
            %   + i*atan(z/z0)) with the coordinate z along and rho across
            %   the beam axis, the Rayleigh range z0 = k*w0^2/2, the beam
            %   width w(z) = w0*sqrt(1 + (z/z0)^2) and the curvature
            %   1/R(z) = z/(z^2 + z0^2) of the wave fronts.
            k = obj.omega;
            t = [cos(obj.theta),sin(obj.theta)];
            z = x*t';
            rho = x*[-t(2);t(1)];
            z0 = k*obj.w0^2/2;
            d = z.^2 + z0^2;
            q = z./d;              % 1/R = w'/w
            w2 = obj.w0^2*d/z0^2;  % w^2
            val = z0./sqrt(d).*exp(-rho.^2./w2 - 1i*k*z - 1i*k*rho.^2.*q/2 ...
                + 1i*atan(z/z0));
            if nargout > 1
                % derivatives of log(u) along and across the beam axis
                dz = -q + 2*rho.^2.*q./w2 - 1i*k ...
                    - 1i*k*rho.^2.*(z0^2 - z.^2)./(2*d.^2) + 1i*z0./d;
                drho = -2*rho./w2 - 1i*k*rho.*q;
                grad = val.*(dz*t + drho*[-t(2),t(1)]);
            end
        end

    end

end