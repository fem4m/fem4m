classdef problemHelmholtzSquareReflection < ProblemHelmholtzTemplateExact
    %problemHelmholtzSquareReflection - Scattering with reflection/transmission
    %   problemHelmholtzSquareReflection is a Helmholtz problem on the
    %   square (-1,1)^2 with the interface x2 = 0 between the medium of
    %   refractive index n1 = 2 (x2 <= 0) and the medium of refractive index
    %   n2 = 1 (x2 > 0). The exact solution is the total field of a plane
    %   wave in the direction (cos(THETA),sin(THETA)) that comes from the
    %   medium n1: the incident plus the reflected wave for x2 <= 0 and the
    %   transmitted wave for x2 > 0, with the reflection coefficient R and
    %   the transmission coefficient T = 1 + R; f = 0. For
    %   n1*cos(THETA) > n2 (THETA < pi/3 for the default indices) the wave is
    %   totally reflected and the field for x2 > 0 is evanescent.
    %
    %   OBJ = problemHelmholtzSquareReflection(K,THETA) sets the wavenumber
    %   K and the angle THETA between the incident direction and the
    %   interface (default deg2rad(29)). A THETA larger than pi/2 is taken
    %   in degrees and converted to radians. K is required.
    %
    %   problemHelmholtzSquareReflection functions:
    %       refraction  - Piecewise-constant refractive index
    %       exactU      - Exact total field
    %       exactGradU  - Gradient of the exact field
    %
    %   problemHelmholtzSquareReflection properties:
    %       geom        - Square geometry
    %       k           - Wavenumber (equal to omega)
    %       theta       - Angle (radians) of the incident wave to the interface
    %       n1          - Refractive index for x2 <= 0 (default 2)
    %       n2          - Refractive index for x2 > 0 (default 1)
    %
    %   See also ProblemHelmholtzTemplateExact, problemHelmholtzPlaneWave

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','midpoint',[0,0],'length',2))
        k
        theta
        n1 = 2
        n2 = 1
    end

    methods

        function obj = problemHelmholtzSquareReflection(k,theta)
            if nargin < 1
                k = [];
            end
            obj = obj@ProblemHelmholtzTemplateExact(k); % for compatibility
            obj.k = k;
            if nargin < 2 || isempty(theta)
                theta = deg2rad(29);
            end
            if theta > pi/2
                obj.theta = deg2rad(theta);
            else
                obj.theta = theta;
            end
        end

        function val = refraction(obj,x)
            val = ones(size(x,1),1)*obj.n2;
            val(x(:,2)<=0) = obj.n1;
        end

        function val = exactU(obj,x)
            y = x(:,2);
            x = x(:,1);
            K1 = obj.k*obj.n1*cos(obj.theta);
            K2 = obj.k*sqrt(obj.n2^2-obj.n1^2*(cos(obj.theta))^2);
            R = -(K2-obj.k*obj.n1*sin(obj.theta))/(K2+obj.k*obj.n1*sin(obj.theta));
            T = 1+R;
            val = T*exp(1i*(K1*x+K2*y));
            val2 = exp(1i*obj.k*obj.n1*(x*cos(obj.theta)+y*sin(obj.theta)))...
                +R*exp(1i*obj.k*obj.n1*(x*cos(obj.theta)-y*sin(obj.theta)));
            val(y<=0) = val2(y<=0);
        end

        function val = exactGradU(obj,x)
            y = x(:,2);
            x = x(:,1);
            K1 = obj.k*obj.n1*cos(obj.theta);
            K2 = obj.k*sqrt(obj.n2^2-obj.n1^2*(cos(obj.theta))^2);
            R = -(K2-obj.k*obj.n1*sin(obj.theta))/(K2+obj.k*obj.n1*sin(obj.theta));
            T = 1+R;
            val = [1i*K1*T*exp(1i*(K1*x+K2*y)),1i*K2*T*exp(1i*(K1*x+K2*y))];
            val2(:,1) = 1i*obj.k*obj.n1*cos(obj.theta)*exp(1i*obj.k*obj.n1*(x*cos(obj.theta)+y*sin(obj.theta)))...
                + 1i*obj.k*obj.n1*cos(obj.theta)*R*exp(1i*obj.k*obj.n1*(x*cos(obj.theta)-y*sin(obj.theta)));
            val2(:,2) = 1i*obj.k*obj.n1*sin(obj.theta)*exp(1i*obj.k*obj.n1*(x*cos(obj.theta)+y*sin(obj.theta)))...
                -1i*obj.k*obj.n1*sin(obj.theta)*R*exp(1i*obj.k*obj.n1*(x*cos(obj.theta)-y*sin(obj.theta))) ;
            val(y<=0,:) = val2(y<=0,:);
        end

    end
end