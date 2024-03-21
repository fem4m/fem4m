classdef problemScalarEllipticSquareKelloggExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticSquareKelloggExact - Kellogg counterexample
    %   problemScalarEllipticSquareKelloggExact is a scalar elliptic problem
    %   on the unit square with the piecewise constant diffusion
    %   coefficient alpha = R*I in the first and third quadrants around
    %   the singular point X0 and alpha = I in the other two, f = 0, and
    %   Kellogg's exact solution u = r^GAMMA*m(phi) in polar coordinates
    %   around X0 with a piecewise trigonometric angular part m. The
    %   Dirichlet data are taken from u on the whole boundary. R, SIGMA
    %   and RHO solve Kellogg's nonlinear system for the exponent GAMMA,
    %   so they cannot be changed one at a time; the defaults are the
    %   values of Bonito et al. (GAMMA = 1/4).
    %
    %   OBJ = problemScalarEllipticSquareKelloggExact() uses the singular
    %   point X0 = [1/2,1/2].
    %
    %   OBJ = problemScalarEllipticSquareKelloggExact(X0) sets the singular
    %   point X0 (1-by-2), where the coefficient interfaces x = X0(1) and
    %   y = X0(2) cross.
    %
    %   problemScalarEllipticSquareKelloggExact functions:
    %       f           - Zero source term
    %       exactU      - Exact piecewise solution
    %       exactGradU  - Exact gradient of the solution
    %       alpha       - Piecewise-constant diffusion tensor
    %
    %   problemScalarEllipticSquareKelloggExact properties:
    %       geom   - Unit square geometry
    %       x0     - Singular point
    %       R      - Diffusion ratio across the interfaces
    %       SIGMA  - Angular parameter sigma of Kellogg's solution
    %       GAMMA  - Singularity exponent
    %       RHO    - Angular parameter rho of Kellogg's solution
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticLshapeExact
    %
    %   Note: Example 5.3 in Morin et al., "Data oscillation and
    %   convergence of adaptive FEM" (GAMMA = 0.1, values in a comment of
    %   the source); Test 2 in Bonito et al., "Adaptive FEM for elliptic
    %   problems with discontinuous coefficients" (the defaults).

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        x0 = [1/2,1/2]
        % Morin et al. values (GAMMA = 0.1):
        %         R = 161.4476387975881
        %         SIGMA = -14.92256510455152
        %         GAMMA = 0.1
        %         RHO = pi/4

        % Bonito et al. values (GAMMA = 1/4, the defaults):
        R = 25.27414236908818
        SIGMA = -5.49778714378214
        GAMMA = 1/4
        RHO = pi/4
    end


    methods

        function obj = problemScalarEllipticSquareKelloggExact(x0)
            if nargin >= 1 && ~isempty(x0)
                obj.x0 = x0;
            end
        end

        function val = f(~,x)
            val = zeros(size(x,1),1);
        end

        function val = exactU(obj,x)
            [phi,r] = cart2pol(x(:,1)-obj.x0(1),x(:,2)-obj.x0(2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            val = zeros(size(x,1),1);
            q1 = phi < pi/2;
            q2 = (phi>= pi/2) & (phi < pi);
            q3 = (phi>=pi) & (phi < 3/2*pi);
            q4 = (phi>=3/2*pi);
            val(q1) = cos((pi/2-obj.SIGMA)*obj.GAMMA)*cos((phi(q1)-pi/2+obj.RHO)*obj.GAMMA);
            val(q2) = cos(obj.RHO*obj.GAMMA)*cos((phi(q2)-pi+obj.SIGMA)*obj.GAMMA);
            val(q3) = cos(obj.SIGMA*obj.GAMMA)*cos((phi(q3)-pi-obj.RHO)*obj.GAMMA);
            val(q4) = cos((pi/2-obj.RHO)*obj.GAMMA)*cos((phi(q4)-3/2*pi-obj.SIGMA)*obj.GAMMA);
            val = r.^(obj.GAMMA).*val;
        end

        function val = exactGradU(obj,x)
            [phi,r] = cart2pol(x(:,1)-obj.x0(1),x(:,2)-obj.x0(2));
            phi( phi<0 ) = phi( phi<0 ) + 2*pi;
            val = zeros(size(x,1),2);
            q1 = phi < pi/2;
            q2 = (phi>= pi/2) & (phi < pi);
            q3 = (phi>=pi) & (phi < 3/2*pi);
            q4 = (phi>=3/2*pi);
            val(q1,:) = [ obj.GAMMA*r(q1).^(obj.GAMMA - 1).*cos(phi(q1) - obj.GAMMA*phi(q1) + (pi*obj.GAMMA)/2 - obj.GAMMA*obj.RHO)*cos(obj.GAMMA*(obj.SIGMA - pi/2)), ...
                obj.GAMMA*r(q1).^(obj.GAMMA - 1).*sin(phi(q1) - obj.GAMMA*phi(q1) + (pi*obj.GAMMA)/2 - obj.GAMMA*obj.RHO)*cos(obj.GAMMA*(obj.SIGMA - pi/2))];
            val(q2,:) = [ obj.GAMMA*r(q2).^(obj.GAMMA - 1).*cos(phi(q2) - obj.GAMMA*phi(q2) + pi*obj.GAMMA - obj.GAMMA*obj.SIGMA)*cos(obj.GAMMA*obj.RHO), ...
                obj.GAMMA*r(q2).^(obj.GAMMA - 1).*sin(phi(q2) - obj.GAMMA*phi(q2) + pi*obj.GAMMA - obj.GAMMA*obj.SIGMA)*cos(obj.GAMMA*obj.RHO)];
            val(q3,:) = [ obj.GAMMA*r(q3).^(obj.GAMMA - 1).*cos(phi(q3) - obj.GAMMA*phi(q3) + pi*obj.GAMMA + obj.GAMMA*obj.RHO)*cos(obj.GAMMA*obj.SIGMA), ...
                obj.GAMMA*r(q3).^(obj.GAMMA - 1).*sin(phi(q3) - obj.GAMMA*phi(q3) + pi*obj.GAMMA + obj.GAMMA*obj.RHO)*cos(obj.GAMMA*obj.SIGMA)];
            val(q4,:) = [ obj.GAMMA*r(q4).^(obj.GAMMA - 1).*cos(phi(q4) - obj.GAMMA*phi(q4) + (3*pi*obj.GAMMA)/2 + obj.GAMMA*obj.SIGMA)*cos(obj.GAMMA*(obj.RHO - pi/2)),...
                obj.GAMMA*r(q4).^(obj.GAMMA - 1).*sin(phi(q4) - obj.GAMMA*phi(q4) + (3*pi*obj.GAMMA)/2 + obj.GAMMA*obj.SIGMA)*cos(obj.GAMMA*(obj.RHO - pi/2))];
        end

        function val = alpha(obj,x)
            y = x(:,2);
            x = x(:,1);
            val = zeros(size(x,1),2,2);
            val(:,1,1) = 1;
            val(:,2,2) = 1;
            indx = (x-obj.x0(1)).*(y-obj.x0(2)) >= 0;
            val(indx,1,1) = obj.R;
            val(indx,2,2) = obj.R;
        end

    end

end