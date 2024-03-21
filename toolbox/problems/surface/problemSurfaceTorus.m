classdef problemSurfaceTorus < ProblemSurfaceTemplate
    %problemSurfaceTorus - Surface of a torus
    %   problemSurfaceTorus describes the torus
    %   (r0 - sqrt(x^2 + y^2))^2 + z^2 = r1^2 with major radius r0
    %   (property, default 1) and tube (minor) radius r1 (property,
    %   default 0.25) via an implicit level-set function, its derivatives,
    %   and a mesh correction.
    %
    %   problemSurfaceTorus functions:
    %       surfaceFunction     - Level-set function of the torus
    %       gradSurfaceFunction - Gradient of the level-set function
    %       d2SurfaceFunction   - Hessian of the level-set function
    %       correct             - Snap the nodes to a 12x8 grid on the torus
    %
    %   problemSurfaceTorus properties:
    %       geom  - Torus geometry
    %       r0    - Major radius (default 1)
    %       r1    - Minor (tube) radius (default 0.25)
    %
    %   See also ProblemSurfaceTemplate, problemSurfaceTorusExact,
    %   problemSurfaceSphere

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','TorusFine'))
        r0 = 1
        r1 = 0.25
    end


    methods

        function val = surfaceFunction(obj,x)
            %surfaceFunction - Level-set function of the torus
            % sqrt((r0 - sqrt(x(:,1).^2+x(:,2).^2)).^2 + x(:,3).^2)-r1
            c1 = sqrt(x(:,1).^2+x(:,2).^2);
            c2 = obj.r0 - c1;
            c3 = sqrt(x(:,3).^2 + c2.^2);
            val = c3-obj.r1;
        end

        function val = gradSurfaceFunction(obj,x)
            c1 = sqrt(x(:,1).^2+x(:,2).^2);
            c2 = obj.r0 - c1;
            c3 = sqrt(x(:,3).^2 + c2.^2);
            val = [-x(:,1).*c2./c1./c3,-x(:,2).*c2./c1./c3,x(:,3)./c3];
        end

        function val = d2SurfaceFunction(obj,x)
            c1 = sqrt(x(:,1).^2+x(:,2).^2);
            c2 = obj.r0 - c1;
            c3 = sqrt(x(:,3).^2 + c2.^2);
            y = x(:,2);
            z = x(:,3);
            x = x(:,1);
            val = zeros(size(x,1),3,3);
            val(:,1,1) = x.^2./(c1.^2.*c3) - c2./(c1.*c3) - (x.^2.*c2.^2)./(c1.^2.*(z.^2 + c2.^2).^(3/2)) + (x.^2.*c2)./(c1.^3.*c3);
            val(:,1,2) = (x.*y)./(c1.^2.*c3) - (x.*y.*c2.^2)./(c1.^2.*(z.^2 + c2.^2).^(3/2)) + (x.*y.*c2)./(c1.^3.*c3);
            val(:,1,3) = (x.*z.*c2)./(c1.*(z.^2 + c2.^2).^(3/2));
            val(:,2,1) = val(:,1,2);
            val(:,2,2) = y.^2./(c1.^2.*c3) - c2./(c1.*c3) - (y.^2.*c2.^2)./(c1.^2.*(z.^2 + c2.^2).^(3/2)) + (y.^2.*c2)./(c1.^3.*c3);
            val(:,2,3) = (y.*z.*c2)./(c1.*(z.^2 + c2.^2).^(3/2));
            val(:,3,1) = val(:,1,3);
            val(:,3,2) = val(:,2,3);
            val(:,3,3) = 1./c3 - z.^2./(z.^2 + c2.^2).^(3/2);
        end

        function mesh = correct(obj,mesh)

            c4n = mesh.c4n;
            indx = abs(c4n(:,3))>0.8*obj.r1;
            c4n(indx,3) = c4n(indx,3)./abs(c4n(indx,3)).*8*obj.r1;
            for i = 1 : mesh.nrNodes
                alpha = (0 : 1/6 : 2-1/6)*pi;
                beta = (0: 1/4 : 2-1/4)*pi;
                normal = -obj.r1*[cos(alpha);sin(alpha);zeros(1,numel(alpha))]';
                dist = 10;
                for curBeta = beta
                    curPoints = obj.r0*[cos(alpha);sin(alpha);zeros(1,numel(alpha))]' + ...
                        normal*cos(curBeta) + ([0,0,obj.r1*sin(curBeta)]'*ones(1,numel(alpha)))';
                    [curDist,indx] = min(sum(((c4n(i,:)'*ones(1,numel(alpha)))' - curPoints).^2,2));
                    if dist > curDist
                        tmp = curPoints(indx,:);
                        tmp(abs(tmp)<=10*eps)=0;
                        dist = curDist;
                    end
                end
                c4n(i,:) = tmp;
            end

            mesh.c4n = c4n;
        end

    end

end