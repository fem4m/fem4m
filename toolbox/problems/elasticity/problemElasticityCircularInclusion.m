classdef problemElasticityCircularInclusion < ProblemElasticityTemplateExact
    %problemElasticityCircularInclusion - Plate with a rigid circular inclusion
    %   OBJ = problemElasticityCircularInclusion(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies a linear elasticity problem with known solution: a rigid
    %   circular inclusion of radius 1/4 at the origin in an infinite
    %   plate in plane strain under unit uniaxial tension in the x
    %   direction (the stress tends to [s11 s22 s12] = [1 0 0] far away),
    %   restricted to the unit square minus the quarter disk. The displacement
    %   is zero on the arc (Dirichlet), the sides x = 0 and y = 0 are symmetry
    %   boundaries (u1 = 0 and u2 = 0), and the sides x = 1 and y = 1 carry the
    %   traction of the exact stress (Neumann). The body force is zero. The
    %   material parameters NAME1,VALUE1,NAME2,VALUE2 are given as in
    %   ElasticityParameters and are required.
    %
    %   problemElasticityCircularInclusion functions:
    %       f            - Zero body force
    %       u4Db         - Zero displacement data
    %       exactU       - Exact displacement
    %       exactSigma   - Exact stress
    %       exactGradU   - Exact displacement gradient
    %       correctNodes - Project the Dirichlet nodes onto the circle of
    %                      radius 1/4
    %       generateMesh - Mesh arrays of the unit square minus the quarter
    %                      disk of radius 1/4
    %
    %   problemElasticityCircularInclusion properties:
    %       geom         - Unit square minus the quarter disk
    %                      (CircularInclusion_*.dat); correctNodes moves
    %                      the new arc nodes of refined meshes
    %       scale        - Amplification of the displacement plots
    %
    %   See also ProblemElasticityTemplateExact, ElasticityParameters,
    %   Geometry

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','CircularInclusion','projectBoundary',@(mesh)(problemElasticityCircularInclusion.correctNodes(mesh))))
        scale = 10000
    end


    methods

        function val = f(~,x)
            val = zeros(size(x,1),2);
        end

        function val = u4Db(~,x)
            val = zeros(size(x,1),2);
        end

        function val = exactU(obj,x)
            a = 0.25;
            nu = obj.lambda/(2*(obj.lambda+obj.mu));
            k = 3 - 4*nu;
            delta = 1/k;
            beta = -2*delta;
            gamma = 2*nu-1;
            [phi,r] = cart2pol(x(:,1), x(:,2));
            ur = 1./(8.*obj.mu.*r).*((k-1).*r.^2 + 2*gamma*a^2 + (beta*(k+1)*a^2 + 2*r.^2 + 2*delta*a^4./r.^2).*cos(2*phi));
            uphi = -1./(8*obj.mu*r).*(beta*(k-1)*a^2 + 2*r.^2 - 2*delta*a^4./r.^2).*sin(2*phi);
            val = [ur .* cos(phi) - uphi .* sin(phi), ur .* sin(phi) + uphi .*cos(phi)];
        end

        function val = exactSigma(obj,x)
            a = 0.25;
            nu = obj.lambda/(2*(obj.lambda+obj.mu));
            k = 3 - 4*nu;
            delta = 1/k;
            beta = -2*delta;
            gamma = 2*nu-1;
            [phi,r] = cart2pol(x(:,1), x(:,2));

            ur = 1/2*(1-gamma*a^2./r.^2+(1-2*beta*a^2./r.^2-3*delta*a^4./r.^4).*cos(2*phi));
            uphi = 1/2*(1+gamma*a^2./r.^2-(1-3*delta*a^4./r.^4).*cos(2*phi));
            urphi = -1/2*(1+beta*a^2./r.^2+3*delta*a^4./r.^4).*sin(2*phi);
            v11 = ur .* cos(phi) - urphi .* sin(phi);
            v21 = urphi .* cos(phi) - uphi .* sin(phi);
            v12 = ur .* sin(phi) + urphi .* cos(phi);
            v22 = urphi .* sin(phi) + uphi .* cos(phi);
            val1 = v11.*cos(phi) - v21.*sin(phi);
            val2 = v12.*sin(phi) + v22.*cos(phi);
            val3 = v12.*cos(phi) - v22.*sin(phi);
            val = [val1,val2,val3];
        end

        function val = exactGradU(obj,x)
            a = 0.25;
            nu = obj.lambda/(2*(obj.lambda+obj.mu));
            k = 3 - 4*nu;
            delta = 1/k;
            beta = -2*delta;
            gamma = 2*nu-1;
            [phi,r] = cart2pol(x(:,1), x(:,2));
            ur = 1./(8.*obj.mu.*r).*((k-1).*r.^2 + 2*gamma*a^2 + (beta*(k+1)*a^2 + 2*r.^2 + 2*delta*a^4./r.^2).*cos(2*phi));
            uphi = -1./(8*obj.mu*r).*(beta*(k-1)*a^2 + 2*r.^2 - 2*delta*a^4./r.^2).*sin(2*phi);

            dr_ur = (cos(2.*phi).*(4.*r - (4.*a.^4.*delta)./r.^3) + 2.*r.*(k - 1))./(8.*obj.mu.*r) - (2.*a.^2.*gamma + cos(2.*phi).*(2.*r.^2 + a.^2.*beta.*(k + 1) + (2.*a.^4.*delta)./r.^2) + r.^2.*(k - 1))./(8.*obj.mu.*r.^2);
            dphi_ur = -(sin(2.*phi).*(2.*r.^2 + a.^2.*beta.*(k + 1) + (2.*a.^4.*delta)./r.^2))./(4.*obj.mu.*r);
            dr_uphi = (sin(2.*phi).*(2.*r.^2 + a.^2.*beta.*(k - 1) - (2.*a.^4.*delta)./r.^2))./(8.*obj.mu.*r.^2) - (sin(2.*phi).*(4.*r + (4.*a.^4.*delta)./r.^3))./(8.*obj.mu.*r);
            dphi_uphi = -(cos(2.*phi).*(2.*r.^2 + a.^2.*beta.*(k - 1) - (2.*a.^4.*delta)./r.^2))./(4.*obj.mu.*r);

            cos_phi = cos(phi);
            sin_phi = sin(phi);
            val(:,1,1) = cos_phi.*(cos_phi.*dr_ur-sin_phi.*dr_uphi) - 1./r.*sin_phi.*(cos_phi.*dphi_ur-sin_phi.*ur-sin_phi.*dphi_uphi-cos_phi.*uphi);
            val(:,1,2) = sin_phi.*(cos_phi.*dr_ur-sin_phi.*dr_uphi) + 1./r.*cos_phi.*(cos_phi.*dphi_ur-sin_phi.*ur-sin_phi.*dphi_uphi-cos_phi.*uphi);
            val(:,2,1) = cos_phi.*(sin_phi.*dr_ur+cos_phi.*dr_uphi) - 1./r.*sin_phi.*(sin_phi.*dphi_ur+cos_phi.*ur+cos_phi.*dphi_uphi-sin_phi.*uphi);
            val(:,2,2) = sin_phi.*(sin_phi.*dr_ur+cos_phi.*dr_uphi) + 1./r.*cos_phi.*(sin_phi.*dphi_ur+cos_phi.*ur+cos_phi.*dphi_uphi-sin_phi.*uphi);

        end

    end

    methods(Static)

        function mesh = correctNodes(mesh)
            %correctNodes - Project the Dirichlet nodes onto the circle of radius 1/4
            %   MESH = correctNodes(MESH) moves the first node of every
            %   Dirichlet edge radially onto the circle of radius 1/4,
            %   unless that node lies on x = 0, x = 1, y = 0 or y = 1.
            c4n = mesh.c4n;
            n4Db = mesh.n4Db;
            for curDbEdge = 1 : size(n4Db,1)
                x = c4n(n4Db(curDbEdge,1),1);
                y = c4n(n4Db(curDbEdge,1),2);
                if x~=0 && x~=1 && y~=0 && y~=1
                    c4n(n4Db(curDbEdge,1),:) = c4n(n4Db(curDbEdge,1),:).*(0.25/sqrt(x^2+y^2));
                end
            end
            mesh.c4n = c4n;
        end

        function [c4n,n4e,n4Db,n4Nb,n4Sb] = generateMesh(level)
            %generateMesh - Mesh arrays of the unit square minus the quarter disk of radius 1/4
            %   [C4N,N4E,N4DB,N4NB,N4SB] = generateMesh(LEVEL) refines the
            %   unit square LEVEL times uniformly and once more near the
            %   disk, removes the elements with a vertex inside the disk
            %   and a few further elements next to the arc (selected by
            %   hard-coded coordinates such as 1/32 and 1/64, which match
            %   only some LEVEL), moves the arc nodes onto the circle,
            %   makes the longest edge of every element its first edge and
            %   returns the mesh arrays: N4DB on the arc, N4NB on x = 1 and
            %   y = 1, and N4SB on x = 0 and y = 0 with the fixed axis
            %   (1 on x = 0, 2 on y = 0) in column 3.

            % Init
            mesh = Mesh([0,0;1,0;0,1;1,1],[4 1 2;1 4 3],[1 2;2 4;4 3;3 1],[]);
            % for k = 1: 5
            %     indx = find((sqrt(sum(mesh.c4n(mesh.n4e(:,1),:).^2,2)) < 0.25) |(sqrt(sum(mesh.c4n(mesh.n4e(:,2),:).^2,2)) < 0.25)|(sqrt(sum(mesh.c4n(mesh.n4e(:,3),:).^2,2)) < 0.25));
            %     eta4e = zeros(mesh.nrElems,1);eta4e(indx) = 1;
            %     eta4ed = zeros(mesh.nrEdges,1);eta4ed(mesh.ed4e(eta4e==1)) = 1;
            %     mesh = mesh.refine(eta4ed,'rgb');
            % end
            mesh = mesh.refine('uniform',level);

            indx = find((sqrt(sum(mesh.c4n(mesh.n4e(:,1),:).^2,2)) < 0.25) |(sqrt(sum(mesh.c4n(mesh.n4e(:,2),:).^2,2)) < 0.25)|(sqrt(sum(mesh.c4n(mesh.n4e(:,3),:).^2,2)) < 0.25));
            eta4ed = zeros(mesh.nrEdges,1);
            eta4ed(mesh.ed4e(indx,1)) = 1;
            mesh = mesh.refine(eta4ed,'rgb');

            % clf;plot(mesh)

            c4n = mesh.c4n;
            n4e = mesh.n4e;
            indx1 = find(sqrt(sum(c4n(n4e(:,1),:).^2,2))<0.25);
            indx2 = find(sqrt(sum(c4n(n4e(:,2),:).^2,2))<0.25);
            indx3 = find(sqrt(sum(c4n(n4e(:,3),:).^2,2))<0.25);
            % indx = intersect(indx1,indx2);
            % indx = intersect(indx,indx3);
            % n4e(indx,:)=[];
            n4e([indx1(:);indx2(:);indx3(:)],:) = [];

            indx = c4n(n4e(:,1),1)==0.25 & c4n(n4e(:,1),2)==0 & c4n(n4e(:,3),1)<0.25;
            n4e(indx,:) = [];
            indx = c4n(n4e(:,1),1)==0.25 & c4n(n4e(:,1),2)==1/32& c4n(n4e(:,3),1)<0.25;
            n4e(indx,:) = [];
            indx = c4n(n4e(:,2),2)==0.25 & c4n(n4e(:,2),1)==0 & c4n(n4e(:,3),2)<0.25;
            n4e(indx,:) = [];
            indx = c4n(n4e(:,2),2)==0.25 & c4n(n4e(:,2),1)==1/32& c4n(n4e(:,3),2)<0.25;
            n4e(indx,:) = [];
            indx = c4n(n4e(:,3),1)==(0.25-1/32-1/64) & c4n(n4e(:,3),2)==(3/32+1/64);
            n4e(indx,:) = [];
            indx = c4n(n4e(:,3),2)==(0.25-1/32-1/64) & c4n(n4e(:,3),1)==(3/32+1/64);
            n4e(indx,:) = [];

            % delete coordinates
            while( ~isempty(setdiff(1:size(c4n,1),unique(n4e(:)))))
                I =  setdiff(1:size(c4n,1),unique(n4e(:)));
                I = I(1);
                c4n(I,:) = [];
                n4e(n4e>I) = n4e(n4e>I) - 1;
            end


            %find the boundary edges (split into Db, Nb and Sb below)
            allEdges = [n4e(:,[1 2]);n4e(:,[2 3]);n4e(:,[3 1])];
            [allEdges,indx] = sort(allEdges,2);
            [allEdges,indy] = sortrows(allEdges);
            indx = indx(indy,:);

            n4Db = [];
            if ~(allEdges(1,1) == allEdges(2,1) && allEdges(1,2) == allEdges(2,2))
                n4Db = [n4Db;allEdges(1,indx(1,:))];
            end
            for curEdge = 2:size(allEdges,1)-1
                if ~(allEdges(curEdge,1) == allEdges(curEdge+1,1) && allEdges(curEdge,2) == allEdges(curEdge+1,2)) && ...
                        ~(allEdges(curEdge,1) == allEdges(curEdge-1,1) && allEdges(curEdge,2) == allEdges(curEdge-1,2))
                    n4Db = [n4Db;allEdges(curEdge,indx(curEdge,:))]; %#ok<AGROW>
                end
            end
            if ~(allEdges(end-1,1) == allEdges(end,1) && allEdges(end-1,2) == allEdges(end,2))
                n4Db = [n4Db;allEdges(end,indx(end,:))];
            end

            % correct nodes
            for curDbEdge = 1 : size(n4Db,1)
                x = c4n(n4Db(curDbEdge,1),1);
                y = c4n(n4Db(curDbEdge,1),2);
                if x~=0 && x~=1 && y~=0 && y~=1
                    c4n(n4Db(curDbEdge,1),:) = c4n(n4Db(curDbEdge,1),:).*(0.25/sqrt(x^2+y^2));
                end
            end

            % longest edge is refedge
            mesh = Mesh(c4n,n4e,n4Db,[]);
            length4ed = mesh.length4ed;
            len4e = length4ed(mesh.ed4e);
            [~,indx] = sort(len4e,2,'descend');
            for curElem = 1: size(n4e,1)
                I = mod([indx(curElem,1)-1,indx(curElem,1),indx(curElem,1)+1],3)+1;
                n4e(curElem,:) = n4e(curElem,I);
            end

            % NB on x == 1 or y == 1, SB on x == 0 or y == 0
            indx1 = find((c4n(n4Db(:,1),1) == 1 & c4n(n4Db(:,2),1) == 1) | (c4n(n4Db(:,1),2) == 1 & c4n(n4Db(:,2),2) == 1));
            indx2 = find((c4n(n4Db(:,1),1) == 0 & c4n(n4Db(:,2),1) == 0) | (c4n(n4Db(:,1),2) == 0 & c4n(n4Db(:,2),2) == 0));

            n4Nb = n4Db(indx1,:);
            n4Sb = n4Db(indx2,:);
            n4Db([indx1;indx2],:) = [];

            indx1 = c4n(n4Sb(:,1),1) == 0 & c4n(n4Sb(:,2),1) == 0;
            indx2 = c4n(n4Sb(:,1),2) == 0 & c4n(n4Sb(:,2),2) == 0;
            n4Sb(indx1,3) = 1;
            n4Sb(indx2,3) = 2;
        end

    end
end