classdef problemObstacleSquareExact < ProblemScalarEllipticTemplateExact & ProblemObstacleTemplate
    %problemObstacleSquareExact - Obstacle problem with a radial exact solution
    %   problemObstacleSquareExact is an obstacle problem on the square
    %   (-1,1)^2 with the zero obstacle and the radial exact solution
    %   u = (r^2-r0^2)^2 for r > r0 and u = 0 on the contact set r <= r0.
    %
    %   problemObstacleSquareExact functions:
    %       lower       - Zero obstacle
    %       f           - Source term
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemObstacleSquareExact properties:
    %       geom        - Square geometry
    %       r0          - Radius of the contact set (0 < r0 < 1)
    %       d           - Dimension in the formula of f (must be 2)
    %
    %   See also ProblemObstacleTemplate, problemObstacleSquareCusp

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','midpoint',[0,0],'length',2))
        r0 = 0.7 % Radius of the contact set (0 < r0 < 1)
        d = 2 % Dimension in the formula of f (must be 2)
    end


    methods

        function val = lower(~,x)
            val = zeros(size(x,1),1);
        end

        function val = f(obj,x)
            [~,r] = cart2pol(x(:,1),x(:,2));
            val = -4*(2*r.^2 + obj.d*(r.^2-obj.r0^2));
            val(r<=obj.r0) = -8*obj.r0^2*(1 -( r(r<=obj.r0).^2 - obj.r0^2));
        end

        function val = exactU(obj,x)
            [~,r] = cart2pol(x(:,1),x(:,2));
            val = (r.^2-obj.r0^2).^2;
            val(r<=obj.r0) = 0;
        end

        function val = exactGradU(obj,x)
            [phi,r] = cart2pol(x(:,1),x(:,2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            % dx = cos(phi) dr - 1/r sin(phi) dphi
            % dy = sin(phi) dr + 1/r cos(phi) dphi
            val = (4*(r.^2-obj.r0^2).*r).*[cos(phi),sin(phi)];
            val(r<=obj.r0,:) = 0;
        end

    end

end