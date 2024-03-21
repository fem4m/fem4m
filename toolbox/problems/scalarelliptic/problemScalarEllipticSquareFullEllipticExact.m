classdef problemScalarEllipticSquareFullEllipticExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticSquareFullEllipticExact - Diffusion-convection-reaction on a square
    %   problemScalarEllipticSquareFullEllipticExact is a full elliptic
    %   (diffusion-convection-reaction) problem on the unit square with
    %   variable convection, constant reaction, a Dirichlet condition on
    %   the left edge x = 0, Neumann conditions on the other three edges
    %   and the smooth exact solution
    %   u = sin(x^3)*cos(y^pi) + x^8 - y^9 + x^6*y^10.
    %
    %   problemScalarEllipticSquareFullEllipticExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %       beta        - Convection field (5*sin(x+y), 6*cos(x+y))
    %       gamma       - Reaction coefficient 7
    %
    %   problemScalarEllipticSquareFullEllipticExact properties:
    %       geom        - Unit square geometry
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry',['N','N','N','D']))
    end


    methods

        function val = f(~,x)
            y = x(:,2);
            x = x(:,1);
            val = -(56+90.*y.^8).*x.^6-(-9.*sin(x.^3).*cos(y.^pi)+30.*y.^10).*x.^4-6.*cos(x.^3)...
                .*x.*cos(y.^pi)+sin(x.^3).*cos(y.^pi).*(y.^pi).^2.*pi.^2./y.^2+sin(x.^3).*sin(y.^pi)...
                .*y.^pi.*pi.^2./y.^2-sin(x.^3).*sin(y.^pi).*y.^pi.*pi./y.^2+72.*y.^7+5.*sin(x+y).*...
                (3.*cos(x.^3).*x.^2.*cos(y.^pi)+8.*x.^7+6.*x.^5.*y.^10)+6.*cos(x+y).*(-sin(x.^3).*...
                sin(y.^pi).*y.^pi.*pi./y-9.*y.^8+10.*x.^6.*y.^9)+7.*sin(x.^3).*cos(y.^pi)+7.*...
                x.^8-7.*y.^9+7.*x.^6.*y.^10;
        end

        function val = exactU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = sin(x.^3).*cos(y.^pi)+x.^8-y.^9+x.^6.*y.^10;
        end

        function val = exactGradU(~,x)
            y = x(:,2);
            x = x(:,1);
            val = [8.*x.^7+6.*x.^5.*y.^10+3.*cos(x.^3).*x.^2.*cos(y.^pi)...
                ,-9.*y.^8+10.*x.^6.*y.^9-sin(x.^3).*pi.*sin(y.^pi).*y.^(pi-1)];
        end

        function val = beta(~,x)
            val = [5*sin(x(:,1)+x(:,2)), 6*cos(x(:,1)+x(:,2))];
        end

        function val = gamma(~,x)
            val = 7*ones(size(x,1),1);
        end

    end

end