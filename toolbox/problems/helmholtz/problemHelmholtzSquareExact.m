classdef problemHelmholtzSquareExact < ProblemHelmholtzTemplateExact
    %problemHelmholtzSquareExact - Helmholtz with a source at Z0
    %   problemHelmholtzSquareExact is a Helmholtz problem on the unit
    %   square (0,1)^2 with the exact solution besselh(0,OMEGA*|x - Z0|), a
    %   Hankel-function point source at Z0 = [x0,y0], which should lie
    %   outside the closed square.
    %
    %   OBJ = problemHelmholtzSquareExact(OMEGA,Z0) sets the angular
    %   frequency OMEGA and the source location Z0, a vector [x0,y0]
    %   (default [-1/4,0]). OMEGA is required.
    %
    %   problemHelmholtzSquareExact functions:
    %       exactU      - Exact Hankel solution
    %       exactGradU  - Gradient of the exact solution
    %
    %   problemHelmholtzSquareExact properties:
    %       geom        - Square geometry
    %       x0          - x1 coordinate of the source location
    %       y0          - x2 coordinate of the source location
    %
    %   See also ProblemHelmholtzTemplateExact,
    %   problemHelmholtzLshapeSmoothExact, problemHelmholtzPlaneWave

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        x0
        y0
    end

    methods

        function obj = problemHelmholtzSquareExact(omega,z0)
            if nargin < 1
                omega = [];
            end
            obj = obj@ProblemHelmholtzTemplateExact(omega);
            if nargin < 2 || isempty(z0)
                obj.x0 = -1/4;
                obj.y0 = 0;
            else
                obj.x0 = z0(1);
                obj.y0 = z0(2);
            end
        end

        function val = exactU(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = besselh(0,obj.omega*sqrt((x-obj.x0).^2+(y-obj.y0).^2));
        end

        function val = exactGradU(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.omega;
            n = sqrt(sum([x-obj.x0,y-obj.y0].^2,2));
            val = [-k*(x-obj.x0).*besselh(1,k*sqrt((x-obj.x0).^2+(y-obj.y0).^2))./n,...
                -k*(y-obj.y0).*besselh(1,k*sqrt((x-obj.x0).^2+(y-obj.y0).^2))./n];
        end

    end
end