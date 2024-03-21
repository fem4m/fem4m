classdef problemScalarEllipticSlitExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticSlitExact - Elliptic problem on a slit domain
    %   problemScalarEllipticSlitExact is a Poisson problem on the slit
    %   domain (-1,1)^2 \ [0,1]x{0} with the harmonic exact solution
    %   u = r^mu*sin(omega+mu*phi) in polar coordinates, phi in [0,2*pi),
    %   whose gradient is singular at the tip of the slit (the origin) for
    %   mu < 1. The geometry has a Neumann condition on the lower side of
    %   the slit and Dirichlet conditions on its upper side and on the
    %   outer boundary; the defaults mu = 1/4 and omega = 0 give u = 0 on
    %   the upper side and a zero normal derivative on the lower side.
    %
    %   problemScalarEllipticSlitExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticSlitExact properties:
    %       geom        - Slit domain geometry
    %       mu          - Singularity exponent
    %       omega       - Phase shift of the angular part sin(omega+mu*phi)
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticLshape, problemScalarEllipticLshapeExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Slit','bdry',['N','D','D'],'slitdelta',eps))
        mu = 1/4
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
        end

    end

end