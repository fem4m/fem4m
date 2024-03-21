classdef problemObstacleLshapeExact < ProblemScalarEllipticTemplateExact & ProblemObstacleTemplate
    %problemObstacleLshapeExact - Obstacle problem on an L-shape
    %   problemObstacleLshapeExact couples a scalar elliptic operator on the
    %   L-shape (-2,2)^2 \ [0,2]x[-2,0] with the zero obstacle. The exact
    %   solution r^(2/3)*gamma(r)*sin(2*phi/3) has a singular gradient at
    %   the re-entrant corner (the origin); the C^2 cutoff gamma is 1 for
    %   r <= 1/4 and 0 for r >= 3/4, so the contact set is r >= 3/4.
    %
    %   problemObstacleLshapeExact functions:
    %       lower       - Zero obstacle
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemObstacleLshapeExact properties:
    %       geom        - L-shaped domain geometry
    %
    %   See also ProblemObstacleTemplate, ProblemScalarEllipticTemplateExact,
    %   problemObstacleSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape','bdry','D','length',4))
    end


    methods

        function val = lower(~,x)
            val = zeros(size(x,1),1);
        end

        function val = f(~,x)
            [phi,r] = cart2pol(x(:,1),x(:,2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            gamma2 = ones(size(x,1),1);
            gamma2(r<=5/4) = 0;
            R = 2*(r-1/4);
            dgamma1 = 2*(-30*R.^4 + 60*R.^3 - 30*R.^2);
            dgamma1(R<0) = 0;
            dgamma1(R>=1) = 0;
            ddgamma1 = 4*(-120*R.^3 + 180*R.^2 - 60*R);
            ddgamma1(R<0) = 0;
            ddgamma1(R>=1) = 0;
            val = -r.^(2/3).*sin(2/3*phi).*(dgamma1./r+ddgamma1)...
                -4/3.*r.^(-1/3).*sin(2/3*phi).*dgamma1 ...
                -gamma2;
        end

        function val = exactU(~,x)
            [phi,r] = cart2pol(x(:,1),x(:,2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            R = 2*(r-1/4);
            gamma = -6*R.^5 + 15*R.^4 - 10*R.^3 +1;
            gamma(R<0) = 1;
            gamma(R>=1) = 0;
            val =  r.^(2/3).*gamma.*sin(2/3*phi);
        end

        function val = exactGradU(~,x)
            [phi,r] = cart2pol(x(:,1),x(:,2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            R = 2*(r-1/4);
            gamma = -6*R.^5 + 15*R.^4 - 10*R.^3 +1;
            gamma(R<0) = 1;
            gamma(R>=1) = 0;
            dgamma = 2*(-30*R.^4 + 60*R.^3 - 30*R.^2);
            dgamma(R<0) = 0;
            dgamma(R>=1) = 0;
            dr = 2/3.*r.^(-1/3).*gamma.*sin(2/3*phi)+r.^(2/3).*dgamma.*sin(2/3*phi);
            dphi = 2/3.*r.^(2/3).*gamma.*cos(2/3*phi);
            val = [cos(phi).*dr - 1./r.*sin(phi).*dphi , sin(phi).*dr + 1./r.*cos(phi).*dphi];
        end

    end

end