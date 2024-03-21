classdef HZ1 < AW1
    %HZ1 - Hu-Zhang element for elasticity
    %   OBJ = HZ1(DEGREE,MESH) creates the cubic Hu-Zhang symmetric
    %   H(div)-conforming stress finite element on the triangular mesh
    %   MESH with 21 degrees of freedom per element: the values of
    %   sigma at the vertices (3 per vertex), the zeroth order moments
    %   of n*sigma*n and t*sigma*n and the first order moment of
    %   n*sigma*n on every edge (3 per edge), and the zeroth order
    %   moments of t*sigma*t on the three edges of the element (3 per
    %   element). Unlike AW1, the tangential traction t*sigma*n is
    %   quadratic on every edge, so its first order moment is dropped,
    %   and the element degrees of freedom are the t*sigma*t moments
    %   instead of the element means of sigma. The basis is obtained
    %   from the 30-dimensional space of symmetric cubic matrix fields
    %   by the constraints t*sigma*n in P2 on each edge and div(sigma)
    %   in P1. DEGREE must be 3.
    %
    %   HZ1 functions:
    %       prepareCoefficients - Compute the basis expansion coefficients
    %       boundaryMoments     - Assemble the Neumann boundary moments
    %       intMomentNN         - Evaluate the normal-normal moment integrand
    %       intMomentNT         - Evaluate the normal-tangent moment integrand
    %       intMomentNN1        - Evaluate the linear normal-normal integrand
    %
    %   See also AW1, AW0, AWk, FEMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = HZ1(degree,mesh)
            %HZ1 - Create a Hu-Zhang finite element object
            if nargin < 2
                [degree,mesh] = deal([]);
            end
            obj = obj@AW1(degree,mesh);
        end

        %% Calculation of the Hu & Zhang FE transformation
        function obj = prepareCoefficients(obj)
            %prepareCoefficients - Compute the basis expansion coefficients
            %   OBJ = prepareCoefficients(OBJ) solves on every element the
            %   local 30x30 Vandermonde system given by the 9 nodal
            %   values, per edge the zeroth order moments of n*sigma*n and
            %   t*sigma*n and the first order moment of n*sigma*n, the 3
            %   zeroth order edge moments of t*sigma*t (the element
            %   dofs), and the 9 constraints t*sigma*n in P2 on each edge
            %   and div(sigma) in P1 to obtain the expansion of the
            %   Hu-Zhang basis in the P3 basis (basisCoefficients).
            normals4e = obj.mesh.normal4e;
            nrElems = obj.mesh.nrElems;
            gradP14e = obj.gradP1;
            gradP14e = permute(matMul(obj.mesh.transformation,permute(gradP14e,[3 2 1])),[1 3 2]);
            ed4e = obj.mesh.ed4e;
            normal4ed = obj.mesh.normal4ed;
            length4ed = obj.mesh.length4ed;
            % Initialisation
            obj.basisCoefficients = zeros(nrElems,21,30);
            C = zeros(30,30);
            % Nodal degrees of freedom
            C(1:9,1:9) = 60*eye(9,9);

            for curElem = 1 : nrElems

                O = zeros(1,3);
                curNormal = normal4ed(ed4e(curElem,:),:);
                sign = sum(permute(normals4e(curElem,:,:),[2 3 1]).*curNormal,2);

                % Zeroth order edge moments n*sigma*n
                N1 = [curNormal(1,:).^2,2*prod(curNormal(1,:))];
                N2 = [curNormal(2,:).^2,2*prod(curNormal(2,:))];
                N3 = [curNormal(3,:).^2,2*prod(curNormal(3,:))];
                R = [30*N1,30*N1,O;
                    O,30*N2,30*N2;
                    30*N3,O,30*N3];
                S = [10*N1,O,O;
                    O,10*N2,O;
                    O,O,10*N3];
                C(10:12,1:18) = [R,S];

                % Zeroth order edge moments t*sigma*n
                N1 = [-prod(curNormal(1,:)),prod(curNormal(1,:)),curNormal(1,1).^2-curNormal(1,2).^2];
                N2 = [-prod(curNormal(2,:)),prod(curNormal(2,:)),curNormal(2,1).^2-curNormal(2,2).^2];
                N3 = [-prod(curNormal(3,:)),prod(curNormal(3,:)),curNormal(3,1).^2-curNormal(3,2).^2];
                R = [30*N1,30*N1,O;
                    O,30*N2,30*N2;
                    30*N3,O,30*N3];
                S = [10*N1,O,O;
                    O,10*N2,O;
                    O,O,10*N3];
                C(13:15,1:18) = [R,S];

                % First order edge moments n*sigma*n
                M1 = [curNormal(1,:).^2,2*prod(curNormal(1,:))]*sign(1);
                M2 = [curNormal(2,:).^2,2*prod(curNormal(2,:))]*sign(2);
                M3 = [curNormal(3,:).^2,2*prod(curNormal(3,:))]*sign(3);
                R = [-5*M1,5*M1,O;
                    O,-5*M2,5*M2;
                    5*M3,O,-5*M3];
                S1 = [O,O,O;
                    O,O,O;
                    O,O,O];
                S2 = [-M1,O,O;
                    O,-M2,O;
                    O,O,-M3];
                C(16:18,1:27) = [R,S1,S2];

                % Element dofs: zeroth order edge moments t*sigma*t
                N1 = [curNormal(1,2).^2,curNormal(1,1).^2,-2*prod(curNormal(1,:))];
                N2 = [curNormal(2,2).^2,curNormal(2,1).^2,-2*prod(curNormal(2,:))];
                N3 = [curNormal(3,2).^2,curNormal(3,1).^2,-2*prod(curNormal(3,:))];
                R = [30*N1,30*N1,O;
                    O,30*N2,30*N2;
                    30*N3,O,30*N3];
                S = [10*N1,O,O;
                    O,10*N2,O;
                    O,O,10*N3];
                C(19:21,1:18) = [R,S];

                % Constraint t*sigma*n in P_2 (difference to AW1)
                M1 = [-prod(curNormal(1,:)),prod(curNormal(1,:)),curNormal(1,1).^2-curNormal(1,2).^2] *sign(1);
                M2 = [-prod(curNormal(2,:)),prod(curNormal(2,:)),curNormal(2,1).^2-curNormal(2,2).^2] *sign(2);
                M3 = [-prod(curNormal(3,:)),prod(curNormal(3,:)),curNormal(3,1).^2-curNormal(3,2).^2] *sign(3);
                R  = [-M1,O,O;
                    O,-M2,O;
                    O,O,-M3];
                C(22:24,19:27) = R;

                % Constraint div(sigma) in P1
                gradP1X = gradP14e(curElem,:,1);
                gradP1Y = gradP14e(curElem,:,2);
                dxW = [3*(gradP1X(1)-gradP1X(2)), -gradP1X(3), gradP1X(3), gradP1X(3);
                    gradP1X(1), 3*(gradP1X(2)-gradP1X(3)), -gradP1X(1), gradP1X(1);
                    -gradP1X(2), gradP1X(2), 3*(gradP1X(3)-gradP1X(1)), gradP1X(2)];
                dyW = [3*(gradP1Y(1)-gradP1Y(2)), -gradP1Y(3), gradP1Y(3), gradP1Y(3);
                    gradP1Y(1), 3*(gradP1Y(2)-gradP1Y(3)), -gradP1Y(1), gradP1Y(1);
                    -gradP1Y(2), gradP1Y(2), 3*(gradP1Y(3)-gradP1Y(1)), gradP1Y(2)];
                O = zeros(3,1);
                Q = [dxW(:,1),O,dyW(:,1),dxW(:,2),O,dyW(:,2),dxW(:,3),O,dyW(:,3),dxW(:,4),O,dyW(:,4);
                    O,dyW(:,1),dxW(:,1),O,dyW(:,2),dxW(:,2),O,dyW(:,3),dxW(:,3),O,dyW(:,4),dxW(:,4)];
                C(25:30,19:30) = 60*max(length4ed(ed4e(curElem,:)))*Q;

                % Assembly
                x = (C/60)\[eye(21);zeros(9,21)];
                obj.basisCoefficients(curElem,:,:) = x';

            end
        end

        function x = boundaryMoments(obj,g)
            %boundaryMoments - Assemble the Neumann boundary moments
            %   X = boundaryMoments(OBJ,G) computes the zeroth order
            %   normal-normal and normal-tangent and the first order
            %   normal-normal edge moments of the traction function
            %   handle G(PTS,NORMALS) on all Neumann boundary edges
            %   (Gauss rules exact for degree 8 for the zeroth and 12 for
            %   the first order moments) and returns them at the Neumann
            %   edge dofs of the vector X (nrDoF-by-1, zero elsewhere).
            x = zeros(obj.nrDoF,1);
            if ~isempty(obj.mesh.n4Nb)
                % moments
                ed4Nb = obj.mesh.ed4Nb;
                nrNodes = obj.mesh.nrNodes;
                nrEdges = obj.mesh.nrEdges;
                length4ed = obj.mesh.length4ed;
                n4Nb = obj.mesh.n4Nb;
                % order 0 NN
                momentsNN = (1./length4ed(ed4Nb)).*integrate(...
                    @(points,refp,ind)obj.intMomentNN(points,refp,ind,g,obj.mesh),n4Nb,ed4Nb,obj.mesh,8);
                fixNodes0 = 3*nrNodes+ed4Nb;
                x(fixNodes0(:)) = momentsNN(:);
                % order 0 NT
                momentsNT = (1./length4ed(ed4Nb)).*integrate(...
                    @(points,refp,ind)obj.intMomentNT(points,refp,ind,g,obj.mesh),n4Nb,ed4Nb,obj.mesh,8);
                fixNodes1 = 3*nrNodes+nrEdges+ed4Nb;
                x(fixNodes1(:)) = momentsNT(:);
                % order 1 NN
                momentsNN1 = (1./length4ed(ed4Nb)).^2.*integrate(...
                    @(points,refp,ind)obj.intMomentNN1(points,refp,ind,g,obj.mesh),n4Nb,ed4Nb,obj.mesh,12);
                fixNodes2 = 3*nrNodes+2*nrEdges+ed4Nb;
                x(fixNodes2(:)) = momentsNN1(:);
            end
        end

        function val = intMomentNN(~,points,~,~,g,mesh)
            %intMomentNN - Evaluate the normal-normal moment integrand
            %   VAL = intMomentNN(OBJ,POINTS,REFP,IND,G,MESH) evaluates
            %   the normal component of the traction G(POINTS,NORMALS)
            %   on the Neumann edges of MESH, the integrand of the zeroth
            %   order normal-normal moments. REFP and IND are ignored.
            normals = mesh.normal4ed;
            normals = normals(mesh.ed4Nb,:);
            g = g(points,normals);
            g = g(:,1).*normals(:,1) + g(:,2).*normals(:,2);
            val = g;
        end

        function val = intMomentNT(~,points,~,~,g,mesh)
            %intMomentNT - Evaluate the normal-tangent moment integrand
            %   VAL = intMomentNT(OBJ,POINTS,REFP,IND,G,MESH) evaluates
            %   the tangential component of the traction G(POINTS,NORMALS)
            %   on the Neumann edges of MESH, the integrand of the zeroth
            %   order normal-tangent moments. REFP and IND are ignored.
            tangents = mesh.tangent4ed;
            tangents = tangents(mesh.ed4Nb,:);
            normals = mesh.normal4ed;
            normals = normals(mesh.ed4Nb,:);
            g = g(points,normals);
            g = g(:,1).*tangents(:,1) + g(:,2).*tangents(:,2);
            val = g;
        end

        function val = intMomentNN1(~,points,~,~,g,mesh)
            %intMomentNN1 - Evaluate the linear normal-normal integrand
            %   VAL = intMomentNN1(OBJ,POINTS,REFP,IND,G,MESH) evaluates
            %   the normal component of the traction G(POINTS,NORMALS) on
            %   the Neumann edges of MESH weighted by the signed
            %   tangential distance (POINTS-MID)*T of POINTS from the edge
            %   midpoints MID, the integrand of the first order
            %   normal-normal moments. REFP and IND are ignored.
            tangents = mesh.tangent4ed;
            tangents = tangents(mesh.ed4Nb,:);
            normals = mesh.normal4ed;
            normals = normals(mesh.ed4Nb,:);
            mid = mesh.point4ed(1/2);
            mid = mid(mesh.ed4Nb,:);
            g = g(points,normals);
            g = g(:,1).*normals(:,1) + g(:,2).*normals(:,2);
            val = (sum((points-mid).*tangents,2)).*g;
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of the 3 nodal stress values, 3
            %   edge moments, and 3 element moments.
            obj.nrBasis4n = 3;
            obj.nrBasis4ed = 3;
            obj.nrBasis4e = 3;
        end

        function byEntity = dofOrder(~)
            %dofOrder - Local order of the dofs per entity type
            %   BYENTITY = dofOrder(OBJ) is [true false false]: the nodal
            %   dofs are listed node by node, the edge moments index by
            %   index.
            byEntity = [true false false];
        end

        function dofs = symmetryDoFs(obj,edgeDoF)
            %symmetryDoFs - Edge moments fixed on the symmetry boundary
            %   DOFS = symmetryDoFs(OBJ,EDGEDOF) returns the second edge
            %   moments EDGEDOF(ED,2), the zeroth order normal-tangent
            %   moments, of all symmetry boundary edges ED.
            ed4Sb = obj.mesh.ed4Sb;
            flag = obj.mesh.n4Sb(:,3);
            dofs = [edgeDoF(ed4Sb(flag==2),2);edgeDoF(ed4Sb(flag==1),2)];
        end

    end
end