classdef problemScalarEllipticSquareWaterfallExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticSquareWaterfallExact - Internal layer on a square
    %   problemScalarEllipticSquareWaterfallExact is a Poisson problem on
    %   the unit square with the "waterfall" exact solution
    %   u = x*y*(1-x)*(1-y)*atan(KAPPA*(r-1)), r = |(x,y) - (5/4,-1/4)|,
    %   which has a steep internal layer along a circular arc; KAPPA
    %   controls the steepness of the layer.
    %
    %   OBJ = problemScalarEllipticSquareWaterfallExact() uses KAPPA = 60.
    %
    %   OBJ = problemScalarEllipticSquareWaterfallExact(KAPPA) sets the
    %   steepness KAPPA.
    %
    %   problemScalarEllipticSquareWaterfallExact functions:
    %       f           - Source term of the waterfall solution
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticSquareWaterfallExact properties:
    %       geom        - Unit square geometry
    %       kappa       - Steepness of the internal layer
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticSquareMountainValleyExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        kappa
    end


    methods

        function obj = problemScalarEllipticSquareWaterfallExact(kappa)
            if nargin < 1
                kappa = 60;
            end
            obj.kappa = kappa;
        end

        function val = f(obj,x)
            k = obj.kappa;
            y = x(:,2);
            x = x(:,1);
            val = 2.*y.*(1-y).*atan(k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-1./4.*y.*(1-x).*(1-y).*k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./4.*x.*y.*(1-y).*k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./16.*x.*y.*(1-x).*(1-y).*k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(3./2).*(32.*x-40).^2./(1+k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)-8.*x.*y.*(1-x).*(1-y).*k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)./(1+k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./32.*x.*y.*(1-x).*(1-y).*k.^3./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).*(32.*x-40).^2./(1+k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2).^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1)+2.*x.*(1-x).*atan(k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-1./4.*x.*(1-x).*(1-y).*k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./4.*x.*y.*(1-x).*k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./16.*x.*y.*(1-x).*(1-y).*k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(3./2).*(32.*y+8).^2./(1+k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+1./32.*x.*y.*(1-x).*(1-y).*k.^3./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).*(32.*y+8).^2./(1+k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2).^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1);
        end

        function val = exactU(obj,x)
            k = obj.kappa;
            y = x(:,2);
            x = x(:,1);
            val = x.*y.*(1-x).*(1-y).*atan(k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1));
        end

        function val = exactGradU(obj,x)
            k = obj.kappa;
            y = x(:,2);
            x = x(:,1);
            val = [y.*(1-x).*(1-y).*atan(k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-y.*(1-y).*atan(k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1)).*x+1./8.*x.*y.*(1-x).*(1-y).*k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*x-40)./(1+k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+x-x+y-y,...
                x.*(1-x).*(1-y).*atan(k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))-x.*y.*(1-x).*atan(k.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1))+1./8.*x.*y.*(1-x).*(1-y).*k./(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2).*(32.*y+8)./(1+k.^2.*(1./4.*(16.*x.^2-40.*x+26+16.*y.^2+8.*y).^(1./2)-1).^2)+x-x+y-y];

        end

    end

end