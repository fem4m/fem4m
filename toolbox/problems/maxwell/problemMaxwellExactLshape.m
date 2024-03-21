classdef problemMaxwellExactLshape
    %problemMaxwellExactLshape - Maxwell on an L-shape, singular solution
    %   problemMaxwellExactLshape is a Maxwell source problem on the
    %   L-shape (-1,1)^2 \ [0,1]^2 (re-entrant corner at the origin, first
    %   quadrant removed) with a singular exact field for convergence
    %   studies. The field behaves like r^(-1/3) at the re-entrant corner
    %   (the rotated gradient of an r^(2/3) function), and its tangential
    %   trace vanishes on the boundary.
    %
    %   OBJ = problemMaxwellExactLshape() uses ALPHA = 1.
    %
    %   OBJ = problemMaxwellExactLshape(ALPHA) sets the coefficient ALPHA
    %   of the zero-order term in curl(rot(u)) + ALPHA*u = f.
    %
    %   problemMaxwellExactLshape functions:
    %       f           - Source current
    %       rotf        - Rot of the source current
    %       u4Db        - Zero Dirichlet field
    %       g           - Zero Neumann-type boundary data
    %       exactU      - Exact singular solution
    %       exactRotU   - Rot of the exact solution
    %
    %   problemMaxwellExactLshape properties:
    %       geom        - L-shaped geometry
    %       m           - Number of holes of the domain (0)
    %       alpha       - Coefficient of the zero-order term
    %
    %   See also problemMaxwellExactSquare, problemMaxwellExactDblConnected,
    %   problemMaxwellLshape3dExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape','rotation',pi/2))
        m = 0
        alpha
    end


    methods

        function obj = problemMaxwellExactLshape(alpha)
            if nargin < 1 || isempty(alpha)
                alpha = 1;
            end
            obj.alpha = alpha;
        end

        function val = f(obj,x)
            alpha = obj.alpha;
            [phi,r] = cart2pol(x(:,1), x(:,2));
            phi(phi<=0) = phi(phi<=0)+2.*pi;

            val(:,1) =  ((8.*sin(pi./3 + (7.*phi)./3))./27 + (386.*r.^4.*cos(pi./6 + (5.*phi)./3))./27 - (473.*r.^6.*cos(pi./6 + (5.*phi)./3))./108 + (511.*r.^8.*cos(pi./6 + (5.*phi)./3))./432 - (101.*r.^6.*cos(pi./6 + (11.*phi)./3))./12 + (143.*r.^8.*cos(pi./6 + (11.*phi)./3))./48 + (955.*r.^8.*cos(pi./6 + (17.*phi)./3))./864 + (20.*r.^2.*sin(pi./3 + phi./3))./3 - (91.*r.^4.*sin(pi./3 + phi./3))./3 + (29.*r.^6.*sin(pi./3 + phi./3))./2 + (20.*r.^2.*sin(pi./3 + (7.*phi)./3))./27 - (47.*r.^8.*sin(pi./3 + phi./3))./16 + (19.*r.^4.*sin(pi./3 + (7.*phi)./3))./27 + (455.*r.^6.*sin(pi./3 + (7.*phi)./3))./54 - (193.*r.^8.*sin(pi./3 + (7.*phi)./3))./72 - (r.^4.*sin(pi./3 + (13.*phi)./3))./3 - (25.*r.^6.*sin(pi./3 + (13.*phi)./3))./12 - (r.^8.*sin(pi./3 + (13.*phi)./3))./16 - (r.^4.*sin(pi./3 + (19.*phi)./3))./27 - (5.*r.^6.*sin(pi./3 + (19.*phi)./3))./108 + (169.*r.^8.*sin(pi./3 + (19.*phi)./3))./432 + (r.^8.*sin(pi./3 + (25.*phi)./3))./48 + (r.^8.*sin(pi./3 + (31.*phi)./3))./864)./r.^(7./3) - ((8.*sin(pi./3 + (7.*phi)./3))./27 - (334.*r.^4.*cos(pi./6 + (5.*phi)./3))./27 + (1903.*r.^6.*cos(pi./6 + (5.*phi)./3))./108 - (1757.*r.^8.*cos(pi./6 + (5.*phi)./3))./432 + (173.*r.^6.*cos(pi./6 + (11.*phi)./3))./36 - (251.*r.^8.*cos(pi./6 + (11.*phi)./3))./144 + (1675.*r.^8.*cos(pi./6 + (17.*phi)./3))./864 - (20.*r.^2.*sin(pi./3 + phi./3))./9 + (127.*r.^4.*sin(pi./3 + phi./3))./9 - (89.*r.^6.*sin(pi./3 + phi./3))./6 + (20.*r.^2.*sin(pi./3 + (7.*phi)./3))./27 + (167.*r.^8.*sin(pi./3 + phi./3))./48 + (127.*r.^4.*sin(pi./3 + (7.*phi)./3))./27 - (445.*r.^6.*sin(pi./3 + (7.*phi)./3))./54 + (167.*r.^8.*sin(pi./3 + (7.*phi)./3))./72 + (r.^4.*sin(pi./3 + (13.*phi)./3))./9 + (25.*r.^6.*sin(pi./3 + (13.*phi)./3))./36 - (83.*r.^8.*sin(pi./3 + (13.*phi)./3))./48 - (r.^4.*sin(pi./3 + (19.*phi)./3))./27 - (5.*r.^6.*sin(pi./3 + (19.*phi)./3))./108 - (83.*r.^8.*sin(pi./3 + (19.*phi)./3))./432 - (r.^8.*sin(pi./3 + (25.*phi)./3))./144 + (r.^8.*sin(pi./3 + (31.*phi)./3))./864)./r.^(7./3) + (alpha.*(r.^2.*cos(2.*phi) + r.^2 - 2).^2.*((r.^2.*cos(2.*phi))./2 - r.^2./2 + 1).*(sin(pi./3 + phi./3) + (13.*r.^2.*cos(pi./6 + (5.*phi)./3))./4 - (7.*r.^2.*sin(pi./3 + phi./3))./2 + (r.^2.*sin(pi./3 + (7.*phi)./3))./4))./(6.*r.^(1./3));
            val(:,2) =  ((8.*cos(pi./3 + (7.*phi)./3))./27 - (20.*r.^2.*cos(pi./3 + phi./3))./3 + (91.*r.^4.*cos(pi./3 + phi./3))./3 - (29.*r.^6.*cos(pi./3 + phi./3))./2 + (20.*r.^2.*cos(pi./3 + (7.*phi)./3))./27 + (47.*r.^8.*cos(pi./3 + phi./3))./16 + (19.*r.^4.*cos(pi./3 + (7.*phi)./3))./27 + (455.*r.^6.*cos(pi./3 + (7.*phi)./3))./54 - (193.*r.^8.*cos(pi./3 + (7.*phi)./3))./72 + (r.^4.*cos(pi./3 + (13.*phi)./3))./3 + (25.*r.^6.*cos(pi./3 + (13.*phi)./3))./12 + (r.^8.*cos(pi./3 + (13.*phi)./3))./16 - (r.^4.*cos(pi./3 + (19.*phi)./3))./27 - (5.*r.^6.*cos(pi./3 + (19.*phi)./3))./108 + (169.*r.^8.*cos(pi./3 + (19.*phi)./3))./432 - (r.^8.*cos(pi./3 + (25.*phi)./3))./48 + (r.^8.*cos(pi./3 + (31.*phi)./3))./864 + (386.*r.^4.*sin(pi./6 + (5.*phi)./3))./27 - (473.*r.^6.*sin(pi./6 + (5.*phi)./3))./108 + (511.*r.^8.*sin(pi./6 + (5.*phi)./3))./432 + (101.*r.^6.*sin(pi./6 + (11.*phi)./3))./12 - (143.*r.^8.*sin(pi./6 + (11.*phi)./3))./48 + (955.*r.^8.*sin(pi./6 + (17.*phi)./3))./864)./r.^(7./3) - ((8.*cos(pi./3 + (7.*phi)./3))./27 + (20.*r.^2.*cos(pi./3 + phi./3))./9 - (127.*r.^4.*cos(pi./3 + phi./3))./9 + (89.*r.^6.*cos(pi./3 + phi./3))./6 + (20.*r.^2.*cos(pi./3 + (7.*phi)./3))./27 - (167.*r.^8.*cos(pi./3 + phi./3))./48 + (127.*r.^4.*cos(pi./3 + (7.*phi)./3))./27 - (445.*r.^6.*cos(pi./3 + (7.*phi)./3))./54 + (167.*r.^8.*cos(pi./3 + (7.*phi)./3))./72 - (r.^4.*cos(pi./3 + (13.*phi)./3))./9 - (25.*r.^6.*cos(pi./3 + (13.*phi)./3))./36 + (83.*r.^8.*cos(pi./3 + (13.*phi)./3))./48 - (r.^4.*cos(pi./3 + (19.*phi)./3))./27 - (5.*r.^6.*cos(pi./3 + (19.*phi)./3))./108 - (83.*r.^8.*cos(pi./3 + (19.*phi)./3))./432 + (r.^8.*cos(pi./3 + (25.*phi)./3))./144 + (r.^8.*cos(pi./3 + (31.*phi)./3))./864 - (334.*r.^4.*sin(pi./6 + (5.*phi)./3))./27 + (1903.*r.^6.*sin(pi./6 + (5.*phi)./3))./108 - (1757.*r.^8.*sin(pi./6 + (5.*phi)./3))./432 - (173.*r.^6.*sin(pi./6 + (11.*phi)./3))./36 + (251.*r.^8.*sin(pi./6 + (11.*phi)./3))./144 + (1675.*r.^8.*sin(pi./6 + (17.*phi)./3))./864)./r.^(7./3) - (alpha.*(r.^2.*cos(2.*phi) - r.^2 + 2).^2.*((r.^2.*cos(2.*phi))./2 + r.^2./2 - 1).*((7.*r.^2.*cos(pi./3 + phi./3))./2 - cos(1./3.*pi + 1./3.*phi) + (r.^2.*cos(pi./3 + (7.*phi)./3))./4 + (13.*r.^2.*sin(pi./6 + (5.*phi)./3))./4))./(6.*r.^(1./3));

        end

        function val = rotf(obj,x)
            alpha = obj.alpha;
            [phi,r] = cart2pol(x(:,1), x(:,2));
            phi(phi<=0) = phi(phi<=0)+2.*pi;

            val = r.^(2./3).*((16.*cos(pi./3 + (10.*phi)./3))./3 + (1600.*sin(pi./6 + (2.*phi)./3))./9 + (40.*alpha.*sin(pi./6 + (2.*phi)./3))./3 + (400.*r.^2.*cos(pi./3 + (10.*phi)./3))./9 - (100.*r.^4.*cos(pi./3 + (10.*phi)./3))./3 - (7.*r.^4.*cos(pi./3 + (22.*phi)./3))./9 - (704.*r.^2.*sin(pi./6 + (2.*phi)./3))./3 + 77.*r.^4.*sin(pi./6 + (2.*phi)./3) - (170.*r.^4.*sin(pi./6 + (14.*phi)./3))./9 + (4.*alpha.*r.^2.*cos(pi./3 + (10.*phi)./3))./3 - (25.*alpha.*r.^4.*cos(pi./3 + (10.*phi)./3))./6 + alpha.*r.^6.*cos(pi./3 + (10.*phi)./3) - (alpha.*r.^6.*cos(pi./3 + (22.*phi)./3))./12 - (80.*alpha.*r.^2.*sin(pi./6 + (2.*phi)./3))./3 + 11.*alpha.*r.^4.*sin(pi./6 + (2.*phi)./3) - (7.*alpha.*r.^6.*sin(pi./6 + (2.*phi)./3))./4 - (17.*alpha.*r.^4.*sin(pi./6 + (14.*phi)./3))./6 + (5.*alpha.*r.^6.*sin(pi./6 + (14.*phi)./3))./6);

        end

        function val = u4Db(~,x)
            %u4Db - Zero Dirichlet field
            val = zeros(size(x,1),2);
        end

        function val = g(~,x)
            val = zeros(size(x,1),2);
        end

        function val = exactU(~,x)
            [phi,r] = cart2pol(x(:,1), x(:,2));
            phi(phi<=0) = phi(phi<=0)+2.*pi;

            dphi_u =  -(2.*r.^(2./3).*((r.^2.*cos(2.*phi))./2 - r.^2./2 + 1).*((r.^2.*cos(2.*phi))./2 + r.^2./2 - 1).*(cos(pi./6 + (2.*phi)./3) - r.^2.*cos(pi./6 + (2.*phi)./3) + (r.^4.*cos(pi./6 + (2.*phi)./3))./8 - (13.*r.^4.*cos(pi./6 + (14.*phi)./3))./16 + (11.*r.^4.*sin(pi./3 + (10.*phi)./3))./16))./3;
            dr_u = (2.*sin(pi./6 + (2.*phi)./3).*((r.^2.*cos(2.*phi))./2 - r.^2./2 + 1).*((r.^2.*cos(2.*phi))./2 + r.^2./2 - 1).*((13.*r.^4.*cos(4.*phi))./8 + 7.*r.^2 - (13.*r.^4)./8 - 1))./(3.*r.^(1./3));

            grad(:,1) = cos(phi).*dr_u - 1./r.*sin(phi).*dphi_u;
            grad(:,2) = sin(phi).*dr_u + 1./r.*cos(phi).*dphi_u;
            val = [grad(:,2),-grad(:,1)];
        end

        function val = exactRotU(~,x)

            [phi,r] = cart2pol(x(:,1), x(:,2));
            phi(phi<=0) = phi(phi<=0)+2.*pi;

            val= (cos(phi).*((r.^(5./3).*sin(2.*phi).*(r.^2.*cos(2.*phi) + r.^2 - 2).^2.*(sin(pi./3 + phi./3) + (13.*r.^2.*cos(pi./6 + (5.*phi)./3))./4 - (7.*r.^2.*sin(pi./3 + phi./3))./2 + (r.^2.*sin(pi./3 + (7.*phi)./3))./4))./6 - ((r.^2.*cos(2.*phi) + r.^2 - 2).^2.*((r.^2.*cos(2.*phi))./2 - r.^2./2 + 1).*(cos(pi./3 + phi./3)./3 - (7.*r.^2.*cos(pi./3 + phi./3))./6 + (7.*r.^2.*cos(pi./3 + (7.*phi)./3))./12 - (65.*r.^2.*sin(pi./6 + (5.*phi)./3))./12))./(6.*r.^(1./3)) + (2.*r.^(5./3).*sin(2.*phi).*(r.^2.*cos(2.*phi) + r.^2 - 2).*((r.^2.*cos(2.*phi))./2 - r.^2./2 + 1).*(sin(pi./3 + phi./3) + (13.*r.^2.*cos(pi./6 + (5.*phi)./3))./4 - (7.*r.^2.*sin(pi./3 + phi./3))./2 + (r.^2.*sin(pi./3 + (7.*phi)./3))./4))./3))./r - (sin(phi).*((r.^(5./3).*sin(2.*phi).*(r.^2.*cos(2.*phi) - r.^2 + 2).^2.*((7.*r.^2.*cos(pi./3 + phi./3))./2 - cos(1./3.*pi + 1./3.*phi) + (r.^2.*cos(pi./3 + (7.*phi)./3))./4 + (13.*r.^2.*sin(pi./6 + (5.*phi)./3))./4))./6 - ((r.^2.*cos(2.*phi) - r.^2 + 2).^2.*((r.^2.*cos(2.*phi))./2 + r.^2./2 - 1).*(sin(pi./3 + phi./3)./3 + (65.*r.^2.*cos(pi./6 + (5.*phi)./3))./12 - (7.*r.^2.*sin(pi./3 + phi./3))./6 - (7.*r.^2.*sin(pi./3 + (7.*phi)./3))./12))./(6.*r.^(1./3)) + (2.*r.^(5./3).*sin(2.*phi).*(r.^2.*cos(2.*phi) - r.^2 + 2).*((r.^2.*cos(2.*phi))./2 + r.^2./2 - 1).*((7.*r.^2.*cos(pi./3 + phi./3))./2 - cos(1./3.*pi + 1./3.*phi) + (r.^2.*cos(pi./3 + (7.*phi)./3))./4 + (13.*r.^2.*sin(pi./6 + (5.*phi)./3))./4))./3))./r - (sin(phi).*(r.^2.*cos(2.*phi) + r.^2 - 2).*(16.*sin(pi./3 + phi./3) - 236.*r.^2.*cos(pi./6 + (5.*phi)./3) + 644.*r.^4.*cos(pi./6 + (5.*phi)./3) - (575.*r.^6.*cos(pi./6 + (5.*phi)./3))./4 + 95.*r.^4.*cos(pi./6 + (11.*phi)./3) - (161.*r.^6.*cos(pi./6 + (11.*phi)./3))./2 + (299.*r.^6.*cos(pi./6 + (17.*phi)./3))./4 + 408.*r.^2.*sin(pi./3 + phi./3) - 734.*r.^4.*sin(pi./3 + phi./3) + 161.*r.^6.*sin(pi./3 + phi./3) + 4.*r.^2.*sin(pi./3 + (7.*phi)./3) - 28.*r.^4.*sin(pi./3 + (7.*phi)./3) + (253.*r.^6.*sin(pi./3 + (7.*phi)./3))./4 + 23.*r.^4.*sin(pi./3 + (13.*phi)./3) - (161.*r.^6.*sin(pi./3 + (13.*phi)./3))./2 + (23.*r.^6.*sin(pi./3 + (19.*phi)./3))./4))./(144.*r.^(4./3)) + (cos(phi).*(r.^2.*cos(2.*phi) - r.^2 + 2).*(16.*cos(pi./3 + phi./3) + 408.*r.^2.*cos(pi./3 + phi./3) - 734.*r.^4.*cos(pi./3 + phi./3) + 161.*r.^6.*cos(pi./3 + phi./3) - 4.*r.^2.*cos(pi./3 + (7.*phi)./3) + 28.*r.^4.*cos(pi./3 + (7.*phi)./3) - (253.*r.^6.*cos(pi./3 + (7.*phi)./3))./4 + 23.*r.^4.*cos(pi./3 + (13.*phi)./3) - (161.*r.^6.*cos(pi./3 + (13.*phi)./3))./2 - (23.*r.^6.*cos(pi./3 + (19.*phi)./3))./4 + 236.*r.^2.*sin(pi./6 + (5.*phi)./3) - 644.*r.^4.*sin(pi./6 + (5.*phi)./3) + (575.*r.^6.*sin(pi./6 + (5.*phi)./3))./4 + 95.*r.^4.*sin(pi./6 + (11.*phi)./3) - (161.*r.^6.*sin(pi./6 + (11.*phi)./3))./2 - (299.*r.^6.*sin(pi./6 + (17.*phi)./3))./4))./(144.*r.^(4./3));
        end
    end
end