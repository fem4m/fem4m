classdef problemScalarEllipticLbrickExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticLbrickExact - Elliptic problem on an L-brick
    %   problemScalarEllipticLbrickExact is a Poisson problem on the
    %   L-brick ((-1,1)^2 \ [0,1]x[-1,0]) x (0,1) with the harmonic exact
    %   solution u = r^mu*sin(omega+mu*phi) in the polar coordinates
    %   (r,phi), phi in [0,2*pi), of (x,y). The solution does not depend
    %   on z, its gradient is singular along the re-entrant edge
    %   x = y = 0 for mu < 1, and the Dirichlet data are taken from it on
    %   the whole boundary. The defaults mu = 2/3 and omega = 0 give u = 0
    %   on the two faces at the re-entrant edge.
    %
    %   problemScalarEllipticLbrickExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticLbrickExact properties:
    %       geom        - L-shaped 3d geometry
    %       mu          - Singularity exponent
    %       omega       - Phase shift of the angular part sin(omega+mu*phi)
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticLbrick, problemScalarEllipticLshapeExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape3d'))
        mu = 2/3
        omega = 0
    end


    methods

        function val = f(~,x)
            val = zeros(size(x,1),1);
        end

        function val = exactU(obj,x)
            [phi,r] = cart2pol(x(:,1),x(:,2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            val =  r.^(obj.mu).*sin(obj.omega+obj.mu*phi);
        end

        function val = exactGradU(obj,x)
            [phi,r] = cart2pol(x(:,1),x(:,2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;

            gradPolar = [obj.mu*r.^(obj.mu-1).*sin(obj.omega+obj.mu*phi),...
                obj.mu*r.^(obj.mu).*cos(obj.omega+obj.mu*phi)];

            gradientTransform = reshape([cos(phi),sin(phi),-1./r.*sin(phi),1./r.*cos(phi)],[],2,2);
            val = matMul(gradientTransform,gradPolar);
            val(:,3) = 0;
        end

    end
end