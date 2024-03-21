classdef problemElasticityLshapeExact < ProblemElasticityTemplateExact
    %problemElasticityLshapeExact - Elastic L-shape with singular solution
    %   OBJ = problemElasticityLshapeExact(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies a linear elasticity problem on the L-shape
    %   (-1,1)^2 \ [0,1]x[-1,0], rotated by -3*pi/4 about the re-entrant
    %   corner at the origin, so that its two re-entrant edges lie on the
    %   rays phi = -3*pi/4 and phi = 3*pi/4. The re-entrant edges are
    %   traction-free (Neumann, g = 0), the other edges carry the exact
    %   displacement as Dirichlet data, and f = 0. The exact displacement
    %   is the corner singularity r^alpha*Phi(phi)/(2*mu) with the
    %   material-independent exponent alpha = 0.544483736782; the angular
    %   part Phi depends on lambda and mu through
    %   2*(lambda+2*mu)/(lambda+mu). exactGradU returns 0 at the corner
    %   r = 0, where the gradient is singular; exactD2U and exactRotEpsilon
    %   return Inf or NaN there. The material parameters
    %   NAME1,VALUE1,NAME2,VALUE2 are given as in ElasticityParameters and
    %   are required.
    %
    %   problemElasticityLshapeExact functions:
    %       f               - Zero body force
    %       g               - Zero traction data
    %       exactU          - Exact displacement
    %       exactSigma      - Exact stress
    %       exactGradU      - Exact displacement gradient
    %       exactRotEpsilon - Row-wise rot of the strain
    %       exactD2U        - Second derivatives of the solution
    %
    %   problemElasticityLshapeExact properties:
    %       geom            - Rotated L-shape geometry
    %       scale           - Amplification of the displacement plots
    %
    %   See also ProblemElasticityTemplateExact, problemElasticityCrack,
    %   problemElasticitySquareExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape','rotation',-3/4*pi,'bdry',['N','D']))
        scale = 100
    end


    methods

        function val = f(~,x)
            val = zeros(size(x,1),2);
        end

        function val = g(~,x,~)
            val = zeros(size(x,1),2);
        end

        function val = exactU(obj,x)
            alpha = 0.544483736782;
            c1 = -cos((alpha+1)*3/4*pi)/cos((alpha-1)*3/4*pi);
            c2 = 2*(obj.lambda+2*obj.mu)/(obj.lambda+obj.mu);
            [phi,r] = cart2pol(x(:,1), x(:,2));
            ur = r.^alpha/(2*obj.mu).*(-(alpha+1)*cos((alpha+1)*phi) + (c2-alpha-1)*c1*cos((alpha-1)*phi));
            uphi = r.^alpha/(2*obj.mu).*((alpha+1)*sin((alpha+1)*phi) + (c2+alpha-1)*c1*sin((alpha-1)*phi));
            val = [ur .* cos(phi) - uphi .* sin(phi), ur .* sin(phi) + uphi .* cos(phi)];
        end

        function val = exactSigma(obj,x)
            gradU = exactGradU(obj,x);
            val = ([2*obj.mu+obj.lambda,obj.lambda,0;obj.lambda,2*obj.mu+obj.lambda,0;0,0,obj.mu]...
                *[gradU(:,1,1),gradU(:,2,2),gradU(:,1,2)+gradU(:,2,1)]')';
        end

        function val = exactGradU(obj,x)
            alpha = 0.544483736782;
            c1 = -cos((alpha+1)*3/4*pi)/cos((alpha-1)*3/4*pi);
            c2 = 2*(obj.lambda+2*obj.mu)/(obj.lambda+obj.mu);
            [phi,r] = cart2pol(x(:,1), x(:,2));
            val = zeros(length(r),1,1);
            ind0 = find(r==0);
            ind1 = find(r~=0);
            r(ind0)=[];
            phi(ind0)=[];
            cos_phi = cos(phi);
            sin_phi = sin(phi);
            ur = r.^alpha/(2*obj.mu).*(-(alpha+1)*cos((alpha+1)*phi) + (c2-alpha-1)*c1*cos((alpha-1)*phi));
            uphi = r.^alpha/(2*obj.mu).*((alpha+1)*sin((alpha+1)*phi) + (c2+alpha-1)*c1*sin((alpha-1)*phi));
            dr_ur = alpha*r.^(alpha-1)/(2*obj.mu).*(-(alpha+1)*cos((alpha+1)*phi)+(c2-alpha-1)*c1*cos((alpha-1)*phi));
            dr_uphi = alpha*r.^(alpha-1)/(2*obj.mu).*((alpha+1)*sin((alpha+1)*phi)+(c2+alpha-1)*c1*sin((alpha-1)*phi));
            dphi_ur = r.^(alpha)/(2*obj.mu).*((alpha+1)^2*sin((alpha+1)*phi)-(c2-alpha-1)*c1*(alpha-1)*sin((alpha-1)*phi));
            dphi_uphi = r.^(alpha)/(2*obj.mu).*((alpha+1)^2*cos((alpha+1)*phi)+(c2+alpha-1)*c1*(alpha-1)*cos((alpha-1)*phi));
            val(ind1,1,1) = cos_phi.*(cos_phi.*dr_ur-sin_phi.*dr_uphi) - 1./r.*sin_phi.*(cos_phi.*dphi_ur-sin_phi.*ur-sin_phi.*dphi_uphi-cos_phi.*uphi);
            val(ind1,1,2) = sin_phi.*(cos_phi.*dr_ur-sin_phi.*dr_uphi) + 1./r.*cos_phi.*(cos_phi.*dphi_ur-sin_phi.*ur-sin_phi.*dphi_uphi-cos_phi.*uphi);
            val(ind1,2,1) = cos_phi.*(sin_phi.*dr_ur+cos_phi.*dr_uphi) - 1./r.*sin_phi.*(sin_phi.*dphi_ur+cos_phi.*ur+cos_phi.*dphi_uphi-sin_phi.*uphi);
            val(ind1,2,2) = sin_phi.*(sin_phi.*dr_ur+cos_phi.*dr_uphi) + 1./r.*cos_phi.*(sin_phi.*dphi_ur+cos_phi.*ur+cos_phi.*dphi_uphi-sin_phi.*uphi);
        end

        function val = exactRotEpsilon(obj,x)
            d2 = obj.exactD2U(x);
            val(:,1) = 0.5*(d2(:,1,3)+d2(:,2,1))-d2(:,1,3);
            val(:,2) = d2(:,2,3) - 0.5*(d2(:,1,2)+d2(:,2,3));
        end

        function val = exactD2U(obj,x)
            alpha = 0.544483736782;
            C1 = -cos((alpha+1)*3/4*pi)/cos((alpha-1)*3/4*pi);
            C2 = 2*(obj.lambda+2*obj.mu)/(obj.lambda+obj.mu);
            [t,r] = cart2pol(x(:,1), x(:,2));
            val(:,1,1) = (alpha.*r.^(alpha - 2).*(4.*cos(t.*(alpha - 2)) - 4.*alpha.^2.*cos(t.*(alpha - 2)) + 10.*C1.*cos(2.*t - alpha.*t) - 10.*C1.*cos(4.*t - alpha.*t) + 2.*C1.*cos(t.*(alpha - 2)) + 2.*C1.*cos(t.*(alpha - 4)) + 3.*C1.*alpha.^2.*cos(2.*t - alpha.*t) - 3.*C1.*alpha.^2.*cos(4.*t - alpha.*t) - 3.*C1.*alpha.^2.*cos(t.*(alpha - 2)) - C1.*alpha.^2.*cos(t.*(alpha - 4)) - 2.*C1.*C2.*cos(2.*t - alpha.*t) + 2.*C1.*C2.*cos(4.*t - alpha.*t) - 11.*C1.*alpha.*cos(2.*t - alpha.*t) + 11.*C1.*alpha.*cos(4.*t - alpha.*t) - 2.*C1.*C2.*cos(t.*(alpha - 2)) - 2.*C1.*C2.*cos(t.*(alpha - 4)) - C1.*alpha.*cos(t.*(alpha - 2)) + C1.*alpha.*cos(t.*(alpha - 4)) + C1.*C2.*alpha.*cos(2.*t - alpha.*t) - C1.*C2.*alpha.*cos(4.*t - alpha.*t) + 3.*C1.*C2.*alpha.*cos(t.*(alpha - 2)) + C1.*C2.*alpha.*cos(t.*(alpha - 4))))./(8.*obj.mu);
            val(:,1,2) = (alpha.*r.^(alpha - 2).*(alpha - 1).*(cos(t.*(alpha - 2)) - C1.*cos(t.*(alpha - 2)) - 2.*C1.*cos(t.*(alpha - 4)) + alpha.*cos(t.*(alpha - 2)) - C1.*C2.*cos(t.*(alpha - 2)) + C1.*alpha.*cos(t.*(alpha - 4))))./(2.*obj.mu);
            val(:,1,3) = (alpha.*r.^(alpha - 2).*(alpha - 1).*(sin(t.*(alpha - 2)) + C1.*sin(t.*(alpha - 2)) - 2.*C1.*sin(t.*(alpha - 4)) + alpha.*sin(t.*(alpha - 2)) + C1.*alpha.*sin(t.*(alpha - 4)) - C1.*C2.*sin(t.*(alpha - 2))))./(2.*obj.mu);
            val(:,2,1) = (alpha.*r.^(alpha - 2).*(alpha - 1).*(sin(t.*(alpha - 2)) + C1.*sin(t.*(alpha - 2)) - 2.*C1.*sin(t.*(alpha - 4)) + alpha.*sin(t.*(alpha - 2)) + C1.*alpha.*sin(t.*(alpha - 4)) + C1.*C2.*sin(t.*(alpha - 2))))./(2.*obj.mu);
            val(:,2,2) = -(alpha.*r.^(alpha - 2).*(alpha - 1).*(sin(t.*(alpha - 2)) - 3.*C1.*sin(t.*(alpha - 2)) - 2.*C1.*sin(t.*(alpha - 4)) + alpha.*sin(t.*(alpha - 2)) + C1.*alpha.*sin(t.*(alpha - 4)) + C1.*C2.*sin(t.*(alpha - 2))))./(2.*obj.mu);
            val(:,2,3) = (alpha.*r.^(alpha - 2).*(alpha - 1).*(cos(t.*(alpha - 2)) - C1.*cos(t.*(alpha - 2)) - 2.*C1.*cos(t.*(alpha - 4)) + alpha.*cos(t.*(alpha - 2)) + C1.*C2.*cos(t.*(alpha - 2)) + C1.*alpha.*cos(t.*(alpha - 4))))./(2.*obj.mu);
        end

    end
end