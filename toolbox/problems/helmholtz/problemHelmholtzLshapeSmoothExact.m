classdef problemHelmholtzLshapeSmoothExact < ProblemHelmholtzTemplateExact
    %problemHelmholtzLshapeSmoothExact - Helmholtz with a Hankel solution
    %   problemHelmholtzLshapeSmoothExact is a Helmholtz problem on the
    %   L-shape (-1,1)^2 \ [0,1]x[-1,0] with the exact solution
    %   besselh(0,OMEGA*|x - x0|), a Hankel-function point source at
    %   x0 = (1/4,-1/4) outside the domain, so the solution is smooth in the
    %   domain.
    %
    %   OBJ = problemHelmholtzLshapeSmoothExact(OMEGA) sets the angular
    %   frequency OMEGA (passed to the inherited constructor). OMEGA is
    %   required.
    %
    %   problemHelmholtzLshapeSmoothExact functions:
    %       exactU      - Exact Hankel solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemHelmholtzLshapeSmoothExact properties:
    %       geom        - L-shaped domain geometry
    %
    %   See also ProblemHelmholtzTemplateExact, problemHelmholtzLshapeExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape'))
    end

    methods

        function val = exactU(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = besselh(0,obj.omega*sqrt((x-1/4).^2+(y+1/4).^2));
        end

        function val = exactGradU(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.omega;
            n = sqrt(sum([x-0.25,y+0.25].^2,2));
            val = [-k*(x-0.25).*besselh(1,k*sqrt((x-0.25).^2+(y+0.25).^2))./n,...
                -k*(y+0.25).*besselh(1,k*sqrt((x-0.25).^2+(y+0.25).^2))./n];
        end

    end
end