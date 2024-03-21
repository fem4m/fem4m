classdef AWNC1 < FE
    %AWNC1 - Non-conforming Arnold-Winther element for elasticity
    %   OBJ = AWNC1(DEGREE,MESH) creates the non-conforming quadratic
    %   Arnold-Winther-type stress finite element on the triangular mesh
    %   MESH. Its space consists of the symmetric quadratic matrix fields
    %   whose normal-normal component n*sigma*n is linear on every edge.
    %   The 15 degrees of freedom per element are the zeroth and first
    %   order moments of the traction sigma*n on each edge (4 per edge)
    %   and the mean values of the three stress components on the
    %   element; there are no nodal degrees of freedom. DEGREE must be 2;
    %   the 1 in the name is the degree of the displacements. The basis
    %   coefficients depend on the element geometry (isAffineEquivalent
    %   is false).
    %
    %   AWNC1 functions:
    %       prepareCoefficients - Compute the basis expansion coefficients
    %       getBasis            - Evaluate the non-conforming basis functions
    %       getDivBasis         - Evaluate the divergence of the basis
    %       basisP2             - Evaluate the quadratic reference basis
    %       gradP2              - Evaluate the gradients of the quadratic basis
    %       basisP1             - Evaluate the linear reference basis
    %       gradP1              - Evaluate the gradients of the linear basis
    %       interpolateBoundary - Interpolate Neumann data into dofs
    %       intMoment0          - Evaluate the zeroth order moment integrand
    %       intMoment1          - Evaluate the first order moment integrand
    %
    %   AWNC1 properties:
    %       mesh              - Mesh on which the element is defined
    %       degree            - Polynomial degree of the element
    %       nrBasisFunc       - Number of basis functions per element
    %       dof4e             - Local-to-global dofs per element
    %       nrDoF             - Total number of degrees of freedom
    %       basisCoefficients - Expansion coefficients in the P2 basis
    %       fixedNodes        - Indices of the fixed degrees of freedom
    %       freeNodes         - Indices of the free degrees of freedom
    %       scaling           - Scaling factors of the degrees of freedom
    %       dof4ed            - Local-to-global dofs per edge
    %       nrBasis4n         - Number of nodal basis functions
    %       nrBasis4ed        - Number of edge basis functions
    %       nrBasis4e         - Number of interior basis functions
    %       isRefElemBased    - Basis lives on the reference element
    %
    %   See also AWNC0, GGNC1, AW1, FEMatrix

    % Copyright (c) 2026 by the fem4m authors

    properties
        mesh              % mesh on which the element is defined
        degree            % polynomial degree of the element
        nrBasisFunc       % number of basis functions per element
        dof4e             % local-to-global dofs per element
        nrDoF             % total number of degrees of freedom
        basisCoefficients % expansion coefficients in the P2 basis
        fixedNodes        % indices of the fixed degrees of freedom
        freeNodes         % indices of the free degrees of freedom
        scaling           % scaling factors of the degrees of freedom
        dof4ed            % local-to-global dofs per edge
        nrBasis4n         % number of nodal basis functions
        nrBasis4ed        % number of edge basis functions
        nrBasis4e         % number of interior basis functions
        isRefElemBased = true % basis lives on the reference element
    end

    methods

        function obj = AWNC1(degree,mesh)
            %AWNC1 - Create a non-conforming Arnold-Winther element object
            if nargin < 2 || ~isequal(degree,2)
                error(['fem4m:',class(obj),':invalidDegree'], ...
                    '%s degree must be 2, the stresses are quadratic',class(obj));
            end
            obj.isAffineEquivalent = false;
            obj.mesh = mesh;
            obj.degree = degree;
            obj = obj.prepareCoefficients();
            obj = obj.dofCounts();
            obj = obj.dofManager();
        end

        %% Calculation of the Arnold Winther transformation
        function obj = prepareCoefficients(obj)
            %prepareCoefficients - Compute the basis expansion coefficients
            %   OBJ = prepareCoefficients(OBJ) solves, on every element,
            %   the 18x18 Vandermonde system given by the 15 edge and
            %   element moments and the three constraints that the
            %   quadratic part of n*sigma*n vanishes on each edge, and
            %   stores the expansion of the non-conforming
            %   Arnold-Winther basis in the P2 basis of basisP2 in
            %   basisCoefficients.
            normals4e = obj.mesh.normal4e;
            nrElems = obj.mesh.nrElems;
            ed4e = obj.mesh.ed4e;
            normal4ed = obj.mesh.normal4ed;
            % Initialisation
            obj.basisCoefficients = zeros(nrElems,15,18);
            C = zeros(18,18);
            % Element degrees of freedom
            I = [eye(3,3),eye(3,3),eye(3,3)];
            C(13:15,:) = [20*I,5*I];
            % Element specific transformation
            for curElem = 1 : nrElems
                % Edge degrees of freedom
                curNormal = normal4ed(ed4e(curElem,:),:);
                N1 = [diag(curNormal(1,:)),curNormal(1,[2 1])'];
                N2 = [diag(curNormal(2,:)),curNormal(2,[2 1])'];
                N3 = [diag(curNormal(3,:)),curNormal(3,[2 1])'];
                curNormal = permute(normals4e(curElem,:,:),[2 3 1]);
                M1 = [diag(curNormal(1,:)),curNormal(1,[2 1])'];
                M2 = [diag(curNormal(2,:)),curNormal(2,[2 1])'];
                M3 = [diag(curNormal(3,:)),curNormal(3,[2 1])'];
                O = zeros(2,3);
                R = [30*N1,30*N1,O;-5*M1,5*M1,O;
                    O,30*N2,30*N2;O,-5*M2,5*M2;
                    30*N3,O,30*N3;5*M3,O,-5*M3];
                S = [10*N1,O,O;O,O,O;
                    O,10*N2,O;O,O,O;
                    O,O,10*N3;O,O,O];
                C(1:12,:) = [R,S];
                % Constraints n*sigma*n in P1 on E1,E2,E3
                N1 = [curNormal(1,:).^2,2*curNormal(1,1)*curNormal(1,2)];
                N2 = [curNormal(2,:).^2,2*curNormal(2,1)*curNormal(2,2)];
                N3 = [curNormal(3,:).^2,2*curNormal(3,1)*curNormal(3,2)];
                O = zeros(size(N1));
                C(16:18,10:18) = 60*[N1, O, O; O, N2, O; O, O, N3];
                % Assembly
                x = (C/60)\[eye(15);zeros(3,15)];
                obj.basisCoefficients(curElem,:,:) = x';
            end
        end

        %% Stress basis functions
        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the non-conforming basis functions
            %   VAL = getBasis(OBJ,PTS,X,PARTS) evaluates the symmetric
            %   quadratic non-conforming Arnold-Winther stress basis
            %   functions in all elements at the reference point X
            %   (1-by-2) in Voigt notation [s11 s22 s12], an
            %   nrElems-by-nrBasisFunc-by-3 array. PTS and PARTS are
            %   ignored.
            basisPk = obj.basisP2(x(1),x(2))';
            a = obj.basisCoefficients(:,:,1:3:end);
            b = obj.basisCoefficients(:,:,2:3:end);
            c = obj.basisCoefficients(:,:,3:3:end);
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,3);
            val(:,:,1) = matMul(a,basisPk);
            val(:,:,2) = matMul(b,basisPk);
            val(:,:,3) = matMul(c,basisPk);
        end

        function val = getDivBasis(obj,~,x,~)
            %getDivBasis - Evaluate the divergence of the basis
            %   VAL = getDivBasis(OBJ,PTS,X,PARTS) evaluates the row-wise
            %   divergences of the non-conforming Arnold-Winther stress
            %   basis functions in all elements at the reference point X
            %   (1-by-2), an nrElems-by-nrBasisFunc-by-2 array. PTS and
            %   PARTS are ignored.
            gradBasisPk = obj.gradP2(x(1),x(2));
            gradBasisPk = permute(matMul(obj.mesh.transformation,permute(gradBasisPk,[3 2 1])),[1 3 2]);
            a = obj.basisCoefficients(:,:,1:3:end);
            b = obj.basisCoefficients(:,:,2:3:end);
            c = obj.basisCoefficients(:,:,3:3:end);
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,1);
            val(:,:,1) = matMul(a,gradBasisPk(:,:,1)) + matMul(c,gradBasisPk(:,:,2));
            val(:,:,2) = matMul(c,gradBasisPk(:,:,1)) + matMul(b,gradBasisPk(:,:,2));
        end

        %% P2 basis functions
        function val = basisP2(obj,x,y)
            %basisP2 - Evaluate the quadratic reference basis
            %   VAL = basisP2(OBJ,X,Y) evaluates the six hierarchical
            %   quadratic basis functions lambda_1, lambda_2, lambda_3,
            %   lambda_1*lambda_2, lambda_2*lambda_3, lambda_3*lambda_1
            %   of the reference triangle (lambda_i the barycentric
            %   coordinates of basisP1) at the point (X,Y) and returns
            %   them as a 6-by-1 vector.
            basisP1 = obj.basisP1(x,y);

            val = zeros(6,1);
            val( 1) = basisP1(1);
            val( 2) = basisP1(2);
            val( 3) = basisP1(3);
            val( 4) = basisP1(1).*basisP1(2);
            val( 5) = basisP1(2).*basisP1(3);
            val( 6) = basisP1(3).*basisP1(1);
        end

        function val = gradP2(obj,x,y)
            %gradP2 - Evaluate the gradients of the quadratic basis
            %   VAL = gradP2(OBJ,X,Y) evaluates the gradients with respect
            %   to the reference coordinates of the six quadratic basis
            %   functions of basisP2 at the point (X,Y) and returns them
            %   as a 6-by-2 array.
            gradP1 = obj.gradP1;
            basisP1 = obj.basisP1(x,y);

            val = zeros(6,2);
            val(1,:) = gradP1(1,:);
            val(2,:) = gradP1(2,:);
            val(3,:) = gradP1(3,:);
            val(4,:) = gradP1(1,:)*basisP1(2) + gradP1(2,:)*basisP1(1);
            val(5,:) = gradP1(2,:)*basisP1(3) + gradP1(3,:)*basisP1(2);
            val(6,:) = gradP1(3,:)*basisP1(1) + gradP1(1,:)*basisP1(3);
        end

        %% P1 basis functions
        function val = basisP1(~,x,y)
            %basisP1 - Evaluate the linear reference basis
            %   VAL = basisP1(OBJ,X,Y) evaluates the three linear nodal
            %   basis functions 1-X-Y, X, Y of the reference triangle at
            %   the point (X,Y) and returns them as a 3-by-1 vector. OBJ
            %   is ignored.
            val = zeros(3,1);
            val(1) = 1 - x - y;
            val(2) = x;
            val(3) = y;
        end

        function val = gradP1(~)
            %gradP1 - Evaluate the gradients of the linear basis
            %   VAL = gradP1(OBJ) returns the constant gradients of the
            %   three linear nodal basis functions of the reference
            %   triangle as a 3-by-2 array. OBJ is ignored.
            val = zeros(3,2);
            val(1,:) = [-1,-1];
            val(2,:) = [1,0];
            val(3,:) = [0,1];
        end

        function x = interpolateBoundary(obj,g)
            %interpolateBoundary - Interpolate Neumann data into dofs
            %   X = interpolateBoundary(OBJ,G) returns the vector X
            %   (nrDoF-by-1) whose dofs on the Neumann boundary edges E
            %   are the zeroth and first order moments of the traction G,
            %   the mean value of G and (1/|E|^2)*int_E s*G ds with the
            %   signed tangential distance s from the edge midpoint, and
            %   whose other entries are zero. G is a function handle
            %   G(PTS,NORMALS) that returns the traction sigma*n (N-by-2)
            %   at the points PTS (N-by-2) with the outward unit normals
            %   NORMALS (N-by-2). The moments are computed with a Gauss
            %   rule exact for polynomials of degree 12.
            x = zeros(obj.nrDoF,1);
            if ~isempty(obj.mesh.n4Nb)
                % moments
                ed4Nb = obj.mesh.ed4Nb;
                length4ed = obj.mesh.length4ed;
                % order 0
                moments0 = (1./length4ed(obj.mesh.ed4Nb)*ones(1,2)).*integrate(...
                    @(points,refp,ind)obj.intMoment0(points,refp,ind,g,obj.mesh),obj.mesh.n4Nb,obj.mesh.ed4Nb,obj.mesh,12);
                fixNodes0 = [ed4Nb,obj.mesh.nrEdges+ed4Nb];
                x(fixNodes0(:)) = moments0(:);
                % order 1
                moments1 = (1./length4ed(obj.mesh.ed4Nb)).^2*ones(1,2).*integrate(...
                    @(points,refp,ind)obj.intMoment1(points,refp,ind,g,obj.mesh),obj.mesh.n4Nb,obj.mesh.ed4Nb,obj.mesh,12);
                fixNodes1 = [2*obj.mesh.nrEdges+ed4Nb,3*obj.mesh.nrEdges+ed4Nb];
                x(fixNodes1(:)) = moments1(:);
            end
        end

        function val = intMoment0(~,points,~,~,g,mesh)
            %intMoment0 - Evaluate the zeroth order moment integrand
            %   VAL = intMoment0(OBJ,POINTS,REFP,IND,G,MESH) evaluates
            %   the traction G(POINTS,NORMALS) with the unit normals
            %   NORMALS of the Neumann edges MESH.ed4Nb, one row of POINTS
            %   per Neumann edge, the integrand of the zeroth order edge
            %   moments. OBJ, REFP and IND are ignored.
            normals = mesh.normal4ed;
            normals = normals(mesh.ed4Nb,:);
            val = g(points,normals);
        end

        function val = intMoment1(~,points,~,~,g,mesh)
            %intMoment1 - Evaluate the first order moment integrand
            %   VAL = intMoment1(OBJ,POINTS,REFP,IND,G,MESH) evaluates
            %   the traction G(POINTS,NORMALS) on the Neumann edges
            %   MESH.ed4Nb, one row of POINTS per Neumann edge, weighted
            %   by the signed tangential distance (POINTS-M)*t of POINTS
            %   from the edge midpoints M, with the unit edge tangents t of
            %   tangent4ed, the integrand of the first order edge moments.
            %   OBJ, REFP and IND are ignored.
            tangents = mesh.tangent4ed;
            tangents = tangents(mesh.ed4Nb,:);
            normals = mesh.normal4ed;
            normals = normals(mesh.ed4Nb,:);
            mid = mesh.point4ed(1/2);
            mid = mid(mesh.ed4Nb,:);
            val = (sum((points-mid).*tangents,2)*ones(1,2)).*g(points,normals);
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of no nodal, 4 edge moments, and
            %   3 element moments.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 4;
            obj.nrBasis4e = 3;
        end

        %% Organisation of the degrees of freedom
        function obj = dofManager(obj)
            %dofManager - Assemble the stress degrees of freedom
            %   OBJ = dofManager(OBJ) enumerates the edge-moment and
            %   element stress degrees of freedom of dofCounts with
            %   FE.dofManager (the edge dofs of each element listed edge
            %   by edge), fixes all edge moments on the Neumann boundary
            %   and those of symmetryDoFs on the symmetry boundary, and
            %   sets the scaling of the edge dofs to the square root of
            %   area4ed (the area of the one or two elements adjacent to
            %   the edge) and of the element dofs to the square root of
            %   the element area.
            obj.nrBasisFunc = 3*obj.nrBasis4ed + obj.nrBasis4e;
            obj = dofManager@FE(obj,false,false,[false true false]);
            % edge moments of each edge
            edgeDoF = obj.dof4ed(:,2*obj.nrBasis4n+1:end);
            obj.fixedNodes = reshape(edgeDoF(obj.mesh.ed4Nb,:),[],1);
            if ~isempty(obj.mesh.n4Sb)
                obj.fixedNodes = [obj.fixedNodes;obj.symmetryDoFs(edgeDoF)];
            end
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
            % scale dofs properly
            obj.scaling = [repmat(sqrt(obj.mesh.area4ed),obj.nrBasis4ed,1);...
                zeros(obj.nrBasis4e*obj.mesh.nrElems,1)];
            obj.scaling(obj.dof4e(:,end-obj.nrBasis4e+1:end)) = sqrt(obj.mesh.area4e)*ones(1,obj.nrBasis4e);
        end

        function dofs = symmetryDoFs(obj,edgeDoF)
            %symmetryDoFs - Edge moments fixed on the symmetry boundary
            %   DOFS = symmetryDoFs(OBJ,EDGEDOF) returns the edge moments
            %   EDGEDOF(ED,J) that vanish on the symmetry boundary edges
            %   ED: the odd moments J on the edges fixed in y, where
            %   (sigma*n)_x = 0, and the even moments on the edges fixed
            %   in x, where (sigma*n)_y = 0.
            ed4Sb = obj.mesh.ed4Sb;
            flag = obj.mesh.n4Sb(:,3);
            dofs = [];
            for j = 1 : size(edgeDoF,2)
                if mod(j,2) == 1
                    dofs = [dofs;edgeDoF(ed4Sb(flag==2),j)]; %#ok<AGROW>
                else
                    dofs = [dofs;edgeDoF(ed4Sb(flag==1),j)]; %#ok<AGROW>
                end
            end
        end

    end
end