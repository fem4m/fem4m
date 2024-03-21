classdef problemScalarEllipticSquareMountainValleyExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticSquareMountainValleyExact - Mountain-valley solution on a square
    %   problemScalarEllipticSquareMountainValleyExact is a Poisson
    %   problem on the square (-3,3)^2 with a smooth mountain-valley exact
    %   solution: a Gaussian peak, a Gaussian dip, and the oscillating
    %   ridge term 0.5*tanh(10*sin(x^2-y^2)*cos(x^2+y^2)). The Dirichlet
    %   data are taken from u on the whole boundary.
    %
    %   problemScalarEllipticSquareMountainValleyExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticSquareMountainValleyExact properties:
    %       geom        - Square geometry (-3,3)^2
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticSquareWaterfallExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','midpoint',[0,0],'length',6))
    end

    methods

        function val = f(~,x)
            y = x(:,2);
            x = x(:,1);
            val = 12.*exp(- 3.*(x + 3./10).^2 - 3.*(y + 1./5).^2) - 12.*exp(- 3.*(x - 3./10).^2 - 3.*(y - 2./5).^2) - (20.*cos(2.*x.^2) - 80.*x.^2.*sin(2.*x.^2))./(cosh(10.*sin(2.*x.^2) - 10.*sin(2.*y.^2)) + 1) + (20.*cos(2.*y.^2) - 80.*y.^2.*sin(2.*y.^2))./(cosh(10.*sin(2.*x.^2) - 10.*sin(2.*y.^2)) + 1) - (9.*exp(- 3.*(x + 3./10).^2 - 3.*(y + 1./5).^2).*(10.*x + 3).^2)./25 + (9.*exp(- 3.*(x - 3./10).^2 - 3.*(y - 2./5).^2).*(10.*x - 3).^2)./25 - (36.*exp(- 3.*(x + 3./10).^2 - 3.*(y + 1./5).^2).*(5.*y + 1).^2)./25 + (36.*exp(- 3.*(x - 3./10).^2 - 3.*(y - 2./5).^2).*(5.*y - 2).^2)./25 - 400.*x.^2.*cos(2.*x.^2).^2.*tanh(5.*sin(2.*x.^2) - 5.*sin(2.*y.^2)).*(tanh(5.*sin(2.*x.^2) - 5.*sin(2.*y.^2)).^2 - 1) - 400.*y.^2.*cos(2.*y.^2).^2.*tanh(5.*sin(2.*x.^2) - 5.*sin(2.*y.^2)).*(tanh(5.*sin(2.*x.^2) - 5.*sin(2.*y.^2)).^2 - 1);
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = exp(-3.*((x+0.3).^2+(y+0.2).^2)) - exp(-3.*((x-0.3).^2+(y-0.4).^2)) + 0.5.*tanh(10.*sin(x.^2-y.^2).*cos(x.^2+y.^2));
        end

        function val = exactGradU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [exp(- 3.*(x - 3./10).^2 - 3.*(y - 2./5).^2).*(6.*x - 9./5) - exp(- 3.*(x + 3./10).^2 - 3.*(y + 1./5).^2).*(6.*x + 9./5) - ((20.*x.*cos(x.^2 + y.^2).*cos(x.^2 - y.^2) - 20.*x.*sin(x.^2 + y.^2).*sin(x.^2 - y.^2)).*(tanh(10.*cos(x.^2 + y.^2).*sin(x.^2 - y.^2)).^2 - 1))./2,...
                ((20.*y.*cos(x.^2 + y.^2).*cos(x.^2 - y.^2) + 20.*y.*sin(x.^2 + y.^2).*sin(x.^2 - y.^2)).*(tanh(10.*cos(x.^2 + y.^2).*sin(x.^2 - y.^2)).^2 - 1))./2 - exp(- 3.*(x + 3./10).^2 - 3.*(y + 1./5).^2).*(6.*y + 6./5) + exp(- 3.*(x - 3./10).^2 - 3.*(y - 2./5).^2).*(6.*y - 12./5)];

        end

    end

end