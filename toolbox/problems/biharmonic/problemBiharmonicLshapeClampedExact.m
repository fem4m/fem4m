classdef problemBiharmonicLshapeClampedExact < ProblemBiharmonicTemplateExact
    %problemBiharmonicLshapeClampedExact - Clamped L-shape, singular exact solution
    %   problemBiharmonicLshapeClampedExact defines a clamped biharmonic
    %   problem on the L-shape (-1,1)^2 \ [0,1]x[-1,0] with the singular
    %   exact solution of Grisvard type
    %   u = (x^2-1)^2*(y^2-1)^2*r^(1+z)*w(phi) in polar coordinates
    %   (r,phi), phi in (0,3*pi/2), about the re-entrant corner at the
    %   origin, where z = 0.5445 is the first singular exponent of the
    %   clamped plate at the 3*pi/2 corner and w is the corresponding
    %   angular function. Value and normal derivative of u vanish on the
    %   boundary.
    %
    %   problemBiharmonicLshapeClampedExact functions:
    %       f           - Source term from the singular solution
    %       exactU      - Exact (singular) solution
    %       exactGradU  - Exact gradient of the solution
    %       exactD2     - Exact second derivatives
    %
    %   problemBiharmonicLshapeClampedExact properties:
    %       geom        - L-shaped domain geometry
    %
    %   See also ProblemBiharmonicTemplateExact, problemBiharmonicLshapeClamped

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape'))
    end


    methods

        function val = f(~,x,varargin)
            [phi,r] = cart2pol(x(:,1),x(:,2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            z = 0.544483736782464;

            d1 = cos(phi.*3.0+pi.*z.*(3.0./2.0)-z.*phi);
            d2 = cos(phi.*3.0+pi.*z.*(3.0./2.0)+z.*phi);
            d3 = cos(phi.*3.0-pi.*z.*(3.0./2.0)+z.*phi);
            d4 = cos(phi.*3.0-pi.*z.*(3.0./2.0)-z.*phi);
            d5 = cos(phi.*5.0-pi.*z.*(3.0./2.0)+z.*phi);
            d6 = cos(phi.*5.0+pi.*z.*(3.0./2.0)+z.*phi);
            d7 = cos(phi.*5.0-pi.*z.*(3.0./2.0)-z.*phi);
            d8 = cos(phi.*5.0+pi.*z.*(3.0./2.0)-z.*phi);
            d9 = cos(phi+pi.*z.*(3.0./2.0)-z.*phi);
            d10 = cos(phi-pi.*z.*(3.0./2.0)-z.*phi);
            d11 = cos(phi-pi.*z.*(3.0./2.0)+z.*phi);
            d12 = cos(phi+pi.*z.*(3.0./2.0)+z.*phi);
            d13 = cos(phi.*7.0-pi.*z.*(3.0./2.0)-z.*phi);
            d14 = cos(phi.*9.0-pi.*z.*(3.0./2.0)-z.*phi);
            d15 = cos(phi.*7.0+pi.*z.*(3.0./2.0)-z.*phi);
            d16 = cos(phi.*9.0+pi.*z.*(3.0./2.0)-z.*phi);

            val = -(r.^(z-1.0).*( ...
                r.^4.*d3.*1.6e2-r.^4.*d2.*1.6e2-r.^4.*d4.*1.6e2+r.^4.*d1.*1.6e2-r.^6.*d3.*9.0e1+ ...
                r.^6.*d2.*9.0e1+r.^6.*d4.*9.0e1-r.^6.*d1.*9.0e1-r.^6.*d5.*4.2e1+r.^6.*d6.*4.2e1+ ...
                r.^6.*d7.*4.2e1-r.^6.*d8.*4.2e1+z.*d10.*1.28e2-z.*d9.*1.28e2+z.^2.*d10.*1.28e2- ...
                z.^2.*d9.*3.84e2-z.^3.*d9.*2.56e2+r.^2.*d11.*4.8e2-r.^2.*d12.*4.8e2- ...
                r.^2.*d10.*4.8e2+r.^2.*d9.*4.8e2-r.^4.*d11.*5.76e2+r.^4.*d12.*5.76e2+ ...
                r.^4.*d10.*5.76e2-r.^4.*d9.*5.76e2+r.^6.*d11.*1.8e2-r.^6.*d12.*1.8e2- ...
                r.^6.*d10.*1.8e2+r.^6.*d9.*1.8e2-z.*r.^2.*d11.*5.6e2-z.*r.^2.*d12.*4.0e2- ...
                z.*r.^2.*d10.*7.2e2+z.*r.^2.*d9.*1.68e3+z.*r.^4.*d11.*8.16e2+z.*r.^4.*d12.*3.36e2+ ...
                z.*r.^4.*d10.*4.8e2-z.*r.^4.*d9.*1.632e3-z.*r.^6.*d11.*2.79e2-z.*r.^6.*d12.*8.1e1- ...
                z.*r.^6.*d10.*1.05e2+z.*r.^6.*d9.*4.65e2+z.*r.^2.*d4.*4.8e1-z.*r.^2.*d1.*4.8e1+ ...
                z.*r.^4.*d3.*3.92e2-z.*r.^4.*d2.*7.2e1-z.*r.^2.*d7.*8.0e1+z.*r.^2.*d8.*8.0e1- ...
                z.*r.^4.*d4.*2.4e2-z.*r.^4.*d1.*8.0e1-z.*r.^6.*d3.*2.13e2+z.*r.^6.*d2.*3.3e1- ...
                z.*r.^4.*d7.*1.2e2+z.*r.^4.*d8.*1.2e2+z.*r.^6.*d4.*7.5e1+z.*r.^6.*d1.*1.05e2+ ...
                z.*r.^6.*d5.*7.1e1+z.*r.^6.*d6.*1.3e1+z.*r.^6.*d7.*6.3e1-z.*r.^6.*d8.*1.47e2- ...
                z.*r.^6.*d13.*7.0+z.*r.^6.*d15.*7.0+z.*r.^6.*d14.*9.0-z.*r.^6.*d16.*9.0- ...
                z.^2.*r.^2.*d11.*7.2e2-z.^2.*r.^2.*d12.*8.0e1-z.^2.*r.^2.*d10.*2.4e2+ ...
                z.^2.*r.^2.*d9.*1.68e3-z.^3.*r.^2.*d11.*1.6e2+z.^2.*r.^4.*d11.*6.24e2+ ...
                z.^2.*r.^4.*d12.*4.8e1+z.^3.*r.^2.*d9.*4.8e2+z.^2.*r.^4.*d10.*9.6e1- ...
                z.^2.*r.^4.*d9.*1.056e3+z.^3.*r.^4.*d11.*9.6e1-z.^2.*r.^6.*d11.*1.53e2- ...
                z.^2.*r.^6.*d12.*9.0-z.^3.*r.^4.*d9.*1.92e2-z.^2.*r.^6.*d10.*1.5e1+ ...
                z.^2.*r.^6.*d9.*2.25e2-z.^3.*r.^6.*d11.*1.8e1+z.^3.*r.^6.*d9.*3.0e1+ ...
                z.^2.*r.^2.*d4.*4.8e1+z.^2.*r.^2.*d1.*4.8e1+z.^2.*r.^4.*d3.*1.52e2- ...
                z.^2.*r.^4.*d2.*8.0+z.^3.*r.^2.*d1.*9.6e1+z.^2.*r.^2.*d7.*8.0e1+ ...
                z.^2.*r.^2.*d8.*8.0e1-z.^2.*r.^4.*d4.*8.0e1-z.^2.*r.^4.*d1.*4.0e2+ ...
                z.^3.*r.^4.*d3.*1.6e1-z.^2.*r.^6.*d3.*6.9e1+z.^2.*r.^6.*d2.*3.0- ...
                z.^3.*r.^2.*d8.*1.6e2-z.^3.*r.^4.*d1.*1.6e2-z.^2.*r.^4.*d7.*1.2e2+ ...
                z.^2.*r.^4.*d8.*3.6e2+z.^2.*r.^6.*d4.*1.5e1+z.^2.*r.^6.*d1.*1.35e2- ...
                z.^3.*r.^6.*d3.*6.0+z.^2.*r.^6.*d5.*2.5e1+z.^2.*r.^6.*d6+z.^3.*r.^4.*d8.*2.4e2+ ...
                z.^3.*r.^6.*d1.*3.0e1+z.^2.*r.^6.*d7.*2.1e1-z.^2.*r.^6.*d8.*1.47e2+ ...
                z.^3.*r.^6.*d5.*2.0-z.^3.*r.^6.*d8.*4.2e1-z.^2.*r.^6.*d13.*7.0-z.^2.*r.^6.*d15.*7.0- ...
                z.^3.*r.^6.*d15.*1.4e1-z.^2.*r.^6.*d14.*9.0-z.^2.*r.^6.*d16.*9.0+z.^3.*r.^6.*d16.*1.8e1 ...
                ).*(-1.0./2.0))./(z.^2-1.0);
        end


        function val = exactU(obj,x)
            [phi,r] = cart2pol(x(:,1),x(:,2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            z = 0.544483736782464;
            val =  (r.^2.*(cos(phi)).^2-1).^2.*(r.^2.*(sin(phi)).^2-1).^2.*r.^(1+z).*obj.gfunc(phi);
        end

        function val = exactGradU(obj,x)
            [phi,r] = cart2pol(x(:,1),x(:,2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            z = 0.544483736782464;
            u_dr = (-r.^z.*((r.^4.*cos(4*phi))/8 + r.^2 - r.^4/8 - 1).*...
                (z - (9*r.^4.*cos(4*phi))/8 - r.^2.*z + (r.^4.*z)/8 - 5*r.^2 ...
                + (9*r.^4)/8 - (r.^4.*z.*cos(4*phi))/8 + 1) ).* obj.gfunc(phi);
            u_dphi = (-r.^(z + 5).*sin(4*phi).*(r.^2.*sin(phi).^2 - 1).*((r.^2.*(2*sin(phi).^2 - 1))/2 - r.^2/2 + 1)).*obj.gfunc(phi) + ...
                (r.^2.*cos(phi).^2-1).^2.*(r.^2.*sin(phi).^2-1).^2.*r.^(1+z).*obj.Dgfunc(phi);
            val = [ cos(phi).*u_dr - 1./r.*sin(phi).*u_dphi,...
                sin(phi).*u_dr + 1./r.*cos(phi).*u_dphi];
        end

        function val = exactD2(~,x)
            [phi,r] = cart2pol(x(:,1),x(:,2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            z = 0.544483736782464;

            d1 = r.^2.*sin(phi).^2-1.0;
            d2 = r.^2.*cos(phi).^2-1.0;
            d3 = (sin(pi.*(z-1.0).*(3.0./2.0))./(z-1.0)-sin(pi.*(z+1.0).*(3.0./2.0))./(z+1.0));
            d4 = (cos(pi.*(z-1.0).*(3.0./2.0))-cos(pi.*(z+1.0).*(3.0./2.0)));
            d5 = (sin(phi.*(z-1.0)).*(z-1.0)-sin(phi.*(z+1.0)).*(z+1.0));
            d6 = (cos(phi.*(z-1.0))-cos(phi.*(z+1.0)));
            d7 = sin(phi.*(z-1.0));
            d8 = sin(phi.*(z+1.0));
            d9 = cos(phi.*(z-1.0));
            d10 = cos(phi.*(z+1.0));
            d11 = d7./(z-1.0)-d8./(z+1.0);

            x1 = (d3.*d6-d4.*d11);
            x2 = (d3.*d5+d4.*d6);
            x3 = r.^z.*cos(phi).*sin(phi).*(z+1.0);
            x5 = r.^(z+1.0).*cos(phi).^2;
            x4 = x5.*sin(phi).^2;
            x6 = r.^(z+1.0).*sin(phi).^2;
            x7 = r.^(z+1.0).*cos(phi).*sin(phi).^3;
            x8 = r.^(z+1.0).*cos(phi).*sin(phi);
            x9 = r.^(z+1.0).*cos(phi).^3.*sin(phi);
            x10 = r.^(z+1.0).*sin(phi).^4;
            x11 = r.^(z+1.0).*cos(phi).^4;

            d12 = (x5.*d2.*x1.*d1.^2.*4.0+r.^2.*x11.*x1.*d1.^2.*8.0+r.^2.*x10.*d2.^2.*x1.*8.0+x6.*d2.^2.*x1.*d1.*4.0 ...
                +z.*r.^(z-1.0).*(z+1.0).*d2.^2.*x1.*d1.^2+r.^2.*x4.*d2.*x1.*d1.*3.2e1 ...
                +r.*r.^z.*cos(phi).^2.*(z+1.0).*d2.*x1.*d1.^2.*8.0 ...
                +r.*r.^z.*sin(phi).^2.*(z+1.0).*d2.^2.*x1.*d1.*8.0);
            d13 = (r.^(z+1.0).*x2.*d2.^2.*d1.^2-r.^2.*x8.*d2.^2.*x1.*d1.*4.0+r.^2.*x8.*d2.*x1.*d1.^2.*4.0);
            d14 = (r.^z.*x2.*(z+1.0).*d2.^2.*d1.^2+r.*x5.*x2.*d2.*d1.^2.*4.0+r.*x6.*x2.*d2.^2.*d1.*4.0 ...
                -r.^3.*x7.*d2.^2.*x1.*8.0+r.^3.*x9.*x1.*d1.^2.*8.0+r.^3.*x7.*d2.*x1.*d1.*1.6e1-r.^3.*x9.*d2.*x1.*d1.*1.6e1 ...
                -r.*x8.*d2.^2.*x1.*d1.*8.0+r.*x8.*d2.*x1.*d1.^2.*8.0-r.^2.*x3.*d2.^2.*x1.*d1.*4.0+r.^2.*x3.*d2.*x1.*d1.^2.*4.0);
            d15 = (r.^z.*(z+1.0).*d2.^2.*x1.*d1.^2+r.*x5.*d2.*x1.*d1.^2.*4.0+r.*x6.*d2.^2.*x1.*d1.*4.0);
            d16 = (r.^(z+1.0).*(d4.*d5-d3.*(d9.*(z-1.0).^2-d10.*(z+1.0).^2)).*d2.^2.*d1.^2+r.^4.*x4.*d2.^2.*x1.*8.0 ...
                +r.^4.*x4.*x1.*d1.^2.*8.0+r.^2.*x5.*d2.^2.*x1.*d1.*4.0-r.^2.*x5.*d2.*x1.*d1.^2.*4.0-r.^2.*x6.*d2.^2.*x1.*d1.*4.0 ...
                +r.^2.*x6.*d2.*x1.*d1.^2.*4.0-r.^2.*x8.*x2.*d2.^2.*d1.*8.0+r.^2.*x8.*x2.*d2.*d1.^2.*8.0-r.^4.*x4.*d2.*x1.*d1.*3.2e1);

            val(:,1) = cos(phi).*(cos(phi).*d12-1.0./r.^2.*sin(phi).*d13+(sin(phi).*d14)./r) ...
                +(sin(phi).*( cos(phi).*d14+sin(phi).*d15-(cos(phi).*d13)./r+(sin(phi).*d16)./r))./r;
            val(:,2) = sin(phi).*(sin(phi).*d12+1.0./r.^2.*cos(phi).*d13-(cos(phi).*d14)./r) ...
                +(cos(phi).*(-sin(phi).*d14+cos(phi).*d15+(sin(phi).*d13)./r+(cos(phi).*d16)./r))./r;
            val(:,3) = cos(phi).*(sin(phi).*d12+1.0./r.^2.*cos(phi).*d13-(cos(phi).*d14)./r) ...
                -(sin(phi).*(-sin(phi).*d14+cos(phi).*d15+(sin(phi).*d13)./r+(cos(phi).*d16)./r))./r;

        end

    end

    methods(Static,Access=private)

        function val = gfunc(phi)
            z = 0.544483736782464;
            omega = 3/2*pi;
            val = ( 1./(z-1).*sin((z-1).*omega) - 1./(z+1).*sin((z+1).*omega) )...
                * ( cos((z-1)*phi) - cos((z+1)*phi) )...
                - ( 1./(z-1).*sin((z-1).*phi) - 1./(z+1).*sin((z+1).*phi) )...
                * ( cos((z-1)*omega) - cos((z+1)*omega) );
        end

        function val = Dgfunc(phi)
            z = 0.544483736782464;
            omega = 3/2*pi;
            val = ( 1./(z-1).*sin((z-1).*omega) - 1./(z+1).*sin((z+1).*omega) )...
                * ( -(z-1).*sin((z-1).*phi) + (z+1).*sin((z+1).*phi) )...
                - ( cos((z-1).*phi) - cos((z+1).*phi) )...
                * ( cos((z-1).*omega) - cos((z+1).*omega) );
        end

    end
end