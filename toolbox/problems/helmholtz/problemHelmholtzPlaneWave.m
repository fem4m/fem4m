classdef problemHelmholtzPlaneWave < ProblemHelmholtzTemplateExact
    %problemHelmholtzPlaneWave - Helmholtz with a plane-wave solution
    %   problemHelmholtzPlaneWave is a Helmholtz problem on the unit square
    %   (0,1)^2 with the plane-wave exact solution
    %   u = exp(1i*OMEGA*(cos(ALPHA)*x1 + sin(ALPHA)*x2)) in the direction
    %   (cos(ALPHA),sin(ALPHA)).
    %
    %   OBJ = problemHelmholtzPlaneWave(OMEGA,ALPHA) sets the angular
    %   frequency OMEGA and the propagation angle ALPHA in radians (default
    %   pi/4). OMEGA is required.
    %
    %   problemHelmholtzPlaneWave functions:
    %       exactU      - Exact plane-wave solution
    %       exactGradU  - Gradient of the exact solution
    %
    %   problemHelmholtzPlaneWave properties:
    %       geom        - Square geometry
    %       alpha       - Propagation angle (radians)
    %
    %   See also ProblemHelmholtzTemplateExact, problemHelmholtzSquareExact,
    %   problemHelmholtzSquareReflection

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        alpha
    end

    methods

        function obj = problemHelmholtzPlaneWave(omega,alpha)
            if nargin < 1
                omega = [];
            end
            obj = obj@ProblemHelmholtzTemplateExact(omega);
            if nargin < 2 || isempty(alpha)
                alpha = pi/4;
            end
            obj.alpha = alpha;
        end

        function val = exactU(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.omega;
            val = exp(1i*k*(cos(obj.alpha)*x+sin(obj.alpha)*y));
        end

        function val = exactGradU(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.omega;
            val = [1i*k*cos(obj.alpha)*exp(1i*k*(cos(obj.alpha)*x+sin(obj.alpha)*y)),...
                1i*k*sin(obj.alpha)*exp(1i*k*(cos(obj.alpha)*x+sin(obj.alpha)*y))];
        end

    end
end