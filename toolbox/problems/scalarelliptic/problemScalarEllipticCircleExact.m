classdef problemScalarEllipticCircleExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticCircleExact - Dirichlet problem on a disk
    %   problemScalarEllipticCircleExact is a Poisson problem on the unit
    %   disk with the radial exact solution u = cos(pi*r/2) and
    %   homogeneous Dirichlet data.
    %
    %   problemScalarEllipticCircleExact functions:
    %       u4Db        - Dirichlet data
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticCircleExact properties:
    %       geom        - Unit disk geometry
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticSquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Circle','bdry','D'))
    end


    methods

        function val = u4Db(~,x)
            val = zeros(size(x,1),1);
        end

        function val = f(~,x)
            [~,r] = cart2pol(x(:,1),x(:,2));
            val = pi^2/4*cos(pi/2*r)+pi/2./r.*sin(pi/2*r);
        end

        function val = exactU(~,x)
            [~,r] = cart2pol(x(:,1),x(:,2));
            val = cos(pi/2*r);
        end

        function val = exactGradU(~,x)
            [theta,r] = cart2pol(x(:,1),x(:,2));
            val = [ cos(theta).*(-pi/2*sin(pi/2*r)),...
                sin(theta).*(-pi/2*sin(pi/2*r))];
        end

    end

end