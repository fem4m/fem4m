classdef problemStokesExactSlit < ProblemStokesTemplateExact
    %problemStokesExactSlit - Stokes flow on a slit domain
    %   problemStokesExactSlit is a Stokes problem on the slit domain
    %   (-1,1)^2 minus [0,1]x{0} with an exact singular velocity field
    %   (u ~ sqrt(r), r the distance to the origin) and matching stress.
    %   Pressure and stress scale with the viscosity NU, so the body force
    %   vanishes for every NU.
    %
    %   OBJ = problemStokesExactSlit() creates the problem with the viscosity
    %   NU = 1.
    %
    %   OBJ = problemStokesExactSlit(NU) sets the viscosity NU.
    %
    %   problemStokesExactSlit functions:
    %       f                - Zero body force
    %       exactU           - Exact velocity
    %       exactGradU       - Velocity gradient
    %       exactP           - Exact pressure
    %       exactSigma       - Exact stress 2*NU*eps(u) - p*I as
    %                          [s11 s22 s12]
    %       exactRotEpsilon  - Row-wise rot of the strain
    %
    %   problemStokesExactSlit properties:
    %       geom             - Slit domain geometry
    %
    %   See also ProblemStokesTemplateExact, problemStokesExactCollidingFlow

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Slit'))
    end


    methods

        function val = f(~,x)
            val = zeros(size(x));
        end

        function val = exactU(~,x)
            [phi,r] = cart2pol(x(:,1), x(:,2));
            phi(phi<0) = phi(phi<0)+2*pi;
            ur = 3/2*sqrt(r).*(cos(phi/2)-cos(3/2*phi));
            uphi = 3/2*sqrt(r).*(3*sin(phi/2)-sin(3/2*phi));
            val = [ur,uphi];
        end

        function val = exactGradU(~,x)
            [phi,r] = cart2pol(x(:,1), x(:,2));
            phi(phi<0) = phi(phi<0)+2*pi;
            dr_ur =  3/4./sqrt(r).*(cos(phi/2)-cos(3/2*phi));
            dr_uphi = 3/4./sqrt(r).*(3*sin(phi/2)-sin(3/2*phi));
            dphi_ur = 3/2*sqrt(r).*(-sin(phi/2)/2+3/2*sin(3/2*phi));
            dphi_uphi = 3/2*sqrt(r).*(3/2*cos(phi/2)-3/2*cos(3/2*phi));
            val(:,1,1) = cos(phi).*dr_ur - 1./r.*sin(phi).*dphi_ur;
            val(:,1,2) = sin(phi).*dr_ur + 1./r.*cos(phi).*dphi_ur;
            val(:,2,1) = cos(phi).*dr_uphi - 1./r.*sin(phi).*dphi_uphi;
            val(:,2,2) = sin(phi).*dr_uphi + 1./r.*cos(phi).*dphi_uphi;
        end

        function val = exactP(obj,x)
            [phi,r] = cart2pol(x(:,1), x(:,2));
            phi(phi<0) = phi(phi<0)+2*pi;
            val = -6*obj.nu./sqrt(r).*cos(phi/2);
        end

        function val = exactSigma(obj,x)
            gradU = obj.exactGradU(x);
            p = obj.exactP(x);
            val(:,1) = 2*obj.nu*gradU(:,1,1)-p;
            val(:,2) = 2*obj.nu*gradU(:,2,2)-p;
            val(:,3) = obj.nu*(gradU(:,1,2)+gradU(:,2,1));
        end

        function val = exactRotEpsilon(~,x)
            [phi,r] = cart2pol(x(:,1), x(:,2));
            phi(phi<0) = phi(phi<0)+2*pi;
            dr_ur =  3/4./sqrt(r).*(cos(phi/2)-cos(3/2*phi));
            dr_uphi = 3/4./sqrt(r).*(3*sin(phi/2)-sin(3/2*phi));
            drdr_ur = -3/8./r./sqrt(r).*(cos(phi/2)-cos(3/2*phi));
            drdr_uphi = -3/8./r./sqrt(r).*(3*sin(phi/2)-sin(3/2*phi));
            dphidr_ur = 3/4./sqrt(r).*(-1/2*sin(phi/2)+3/2*sin(3/2*phi));
            dphidr_uphi = 3/4./sqrt(r).*(3/2*cos(phi/2)-3/2*cos(3/2*phi));
            dphi_ur = 3/2*sqrt(r).*(-sin(phi/2)/2+3/2*sin(3/2*phi));
            dphi_uphi = 3/2*sqrt(r).*(3/2*cos(phi/2)-3/2*cos(3/2*phi));
            dphidphi_ur = 3/2*sqrt(r).*(-cos(phi/2)/4+9/4*cos(3/2*phi));
            dphidphi_uphi = 3/2*sqrt(r).*(-3/4*sin(phi/2)+9/4*sin(3/2*phi));
            u2xx = cos(phi).*(cos(phi).*drdr_uphi+1./r.^2.*sin(phi).*dphi_uphi-1./r.*sin(phi).*dphidr_uphi)...
                - 1./r.*sin(phi).*(-sin(phi).*dr_uphi+cos(phi).*dphidr_uphi-1./r.*cos(phi).*dphi_uphi-1./r.*sin(phi).*dphidphi_uphi);
            u1yy = sin(phi).*(sin(phi).*drdr_ur-1./r.^2.*cos(phi).*dphi_ur+1./r.*cos(phi).*dphidr_ur)...
                + 1./r.*cos(phi).*(cos(phi).*dr_ur-1./r.*sin(phi).*dphi_ur+sin(phi).*dphidr_ur+1./r.*cos(phi).*dphidphi_ur);
            u1xy = cos(phi).*(sin(phi).*drdr_ur - 1./r.^2.*cos(phi).*dphi_ur+1./r.*cos(phi).*dphidr_ur) ...
                - 1./r.*sin(phi).*(cos(phi).*dr_ur+sin(phi).*dphidr_ur-1./r.*sin(phi).*dphi_ur+1./r.*cos(phi).*dphidphi_ur);
            u2xy = cos(phi).*(sin(phi).*drdr_uphi - 1./r.^2.*cos(phi).*dphi_uphi+1./r.*cos(phi).*dphidr_uphi) ...
                - 1./r.*sin(phi).*(cos(phi).*dr_uphi+sin(phi).*dphidr_uphi-1./r.*sin(phi).*dphi_uphi+1./r.*cos(phi).*dphidphi_uphi);
            val = 1/2*[u2xx-u1xy,u2xy-u1yy];
        end
    end
end