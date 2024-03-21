classdef problemScalarEllipticLshapeExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticLshapeExact - Elliptic problem on an L-shape
    %   problemScalarEllipticLshapeExact is a Poisson problem on the
    %   L-shape (-1,1)^2 \ [0,1]x[-1,0] with the harmonic exact solution
    %   u = r^mu*sin(omega+mu*phi) in polar coordinates, phi in [0,2*pi),
    %   whose gradient is singular at the re-entrant corner (the origin)
    %   for mu < 1. The defaults mu = 2/3 and omega = 0 give u = 0 on the
    %   two edges at the re-entrant corner.
    %
    %   OBJ = problemScalarEllipticLshapeExact() uses Dirichlet boundary
    %   conditions on the two edges at the re-entrant corner and Neumann
    %   boundary conditions on the remaining edges (BDRY = ['D','N']).
    %
    %   OBJ = problemScalarEllipticLshapeExact(BDRY) sets the boundary
    %   conditions BDRY, one letter 'D' (Dirichlet) or 'N' (Neumann) per
    %   boundary part: part 1 is the two edges at the re-entrant corner,
    %   part 2 the remaining edges (see Geometry).
    %
    %   problemScalarEllipticLshapeExact functions:
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticLshapeExact properties:
    %       geom        - L-shaped domain geometry
    %       mu          - Singularity exponent
    %       omega       - Phase shift of the angular part sin(omega+mu*phi)
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticLshape, problemScalarEllipticLbrickExact,
    %   problemScalarEllipticSlitExact, Geometry

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape','bdry',['D','N']))
        mu = 2/3
        omega = 0
    end


    methods

        function obj = problemScalarEllipticLshapeExact(bdry)
            if nargin > 0 && ~isempty(bdry)
                obj.geom = Geometry(struct('name','Lshape','bdry',bdry));
            end
        end

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