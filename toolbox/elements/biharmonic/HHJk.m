classdef HHJk < FE
    %HHJk - Hellan-Herrmann-Johnson element of degree k
    %   OBJ = HHJk(DEGREE,MESH) creates the Hellan-Herrmann-Johnson
    %   mixed finite element object of polynomial degree DEGREE >= 0 on
    %   the triangular mesh MESH. Its functions are symmetric 2-by-2
    %   matrix fields of degree DEGREE with continuous normal-normal
    %   component n'*sigma*n across the edges, stored in the component
    %   order [sigma_xx, sigma_yy, sigma_xy]. The degrees of freedom are
    %   the coefficients of a hierarchical basis with DEGREE+1 edge
    %   basis functions per edge and 3*DEGREE*(DEGREE+1)/2 interior
    %   basis functions per element. fixedNodes holds the edge degrees
    %   of freedom on the Dirichlet and the Neumann boundary.
    %
    %   This is the HHJ element for plate bending problems. HHJB adds the
    %   divdiv-free bubbles of the TDNNS method for linear elasticity.
    %
    %   HHJk functions:
    %       getBasis               - Evaluate the basis functions
    %       getDivBasis            - Evaluate the divergence of the basis
    %       getDivDivBasis         - Evaluate double divergences of the basis
    %       projectBoundary        - Project data onto boundary edges
    %       intBoundaryMatrix      - Integrate the boundary mass matrix
    %       intBoundaryRHS         - Integrate the boundary load vector
    %       getEdgeRefBasis        - Evaluate the edge reference basis
    %       getEdgeDivRefBasis     - Evaluate the edge div reference basis
    %       getEdgeDivDivRefBasis  - Evaluate the edge divdiv ref basis
    %       getElemRefBasis        - Evaluate the element reference basis
    %       getElemDivRefBasis     - Evaluate the element div ref basis
    %       getElemDivDivRefBasis  - Evaluate the element divdiv ref basis
    %       p1RefBasis             - Evaluate the P1 reference basis
    %
    %   HHJk properties:
    %       dof4e          - Local-to-global dofs per element
    %       dof4ed         - Local-to-global degrees of freedom per edge
    %       nrDoF          - Total number of degrees of freedom
    %       mesh           - Mesh on which the element is defined
    %       freeNodes      - Indices of the free degrees of freedom
    %       fixedNodes     - Indices of the fixed degrees of freedom
    %       degree         - Polynomial degree of the element
    %       nrBasisFunc    - Number of basis functions per element
    %       jacobian       - Transposed Jacobians of the element maps
    %       negativeEdges  - Orientation sign of each edge per element
    %       degree4DoF     - Polynomial degree of each local dof
    %       dof4Db         - Dofs on the Dirichlet boundary
    %       dof4Nb         - Dofs on the Neumann boundary
    %       nrBasis4n      - Number of nodal basis functions
    %       nrBasis4ed     - Number of edge basis functions
    %       nrBasis4e      - Number of interior basis functions
    %       bubbles        - True for the divdiv-free bubbles of HHJB
    %       isRefElemBased - True if the basis lives on the reference element
    %
    %   See also HHJB, TUBAk, Argyris, Pk, FE

    % Copyright (c) 2026 by the fem4m authors

    properties
        dof4e          % local-to-global dofs per element
        dof4ed         % local-to-global degrees of freedom per edge
        nrDoF          % total number of degrees of freedom
        mesh           % mesh on which the element is defined
        freeNodes      % indices of the free degrees of freedom
        fixedNodes     % indices of the fixed degrees of freedom
        degree         % polynomial degree of the element
        nrBasisFunc    % number of basis functions per element
        jacobian % transposed Jacobians of the element maps
        negativeEdges % orientation sign of each edge per element
        degree4DoF     % polynomial degree of each local dof
        dof4Db         % dofs on the Dirichlet boundary
        dof4Nb         % dofs on the Neumann boundary
        nrBasis4n      % number of nodal basis functions
        nrBasis4ed     % number of edge basis functions
        nrBasis4e      % number of interior basis functions
        isRefElemBased = true % true if the basis lives on the reference element
    end

    properties (SetAccess = protected)
        bubbles = false % true for the divdiv-free bubbles of HHJB
    end

    methods

        %% Constructor
        function obj = HHJk(degree,mesh)
            %HHJk - Create a Hellan-Herrmann-Johnson finite element object
            if degree < 0
                error('fem4m:HHJk:invalidDegree','HHJ degree must be >= 0');
            end
            obj.degree = degree;
            obj.mesh = mesh;
            obj = obj.dofCounts();
            obj = dofManager(obj);
            obj.jacobian = obj.mesh.jacobian;
            % orientation
            obj.negativeEdges = zeros(obj.mesh.nrElems,3);
            innerEdge = obj.mesh.e4ed(:,2)>0;
            obj.negativeEdges(:,1) = (obj.mesh.n4e(:,1) > obj.mesh.n4e(:,2)) & innerEdge(obj.mesh.ed4e(:,1));
            obj.negativeEdges(:,2) = (obj.mesh.n4e(:,2) > obj.mesh.n4e(:,3)) & innerEdge(obj.mesh.ed4e(:,2));
            obj.negativeEdges(:,3) = (obj.mesh.n4e(:,3) > obj.mesh.n4e(:,1)) & innerEdge(obj.mesh.ed4e(:,3));
            obj.negativeEdges(obj.negativeEdges==0) = -1;
            obj.negativeEdges = -obj.negativeEdges;
        end

        %% basis functions

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   Piola-mapped HHJ basis functions B*S*B'/det(B)^2 of the
            %   reference basis functions S at the reference point PTS_REF
            %   (1-by-2) in all elements, with the Jacobian B of the element
            %   map. VAL is nrElems-by-nrBasisFunc-by-3 in the component order
            %   [sigma_xx, sigma_yy, sigma_xy]. PTS and PARTS are ignored.
            basis = zeros(obj.mesh.nrElems,obj.nrBasisFunc,3);
            basis(:,1:3*(obj.degree+1),:) = ones(obj.mesh.nrElems,1).*obj.getEdgeRefBasis(obj.degree,x);
            for k = 2 : 2 : obj.degree+1
                basis(:,3*(k-1)+[1,2,3],1) = basis(:,3*(k-1)+[1,2,3],1).*obj.negativeEdges;
                basis(:,3*(k-1)+[1,2,3],2) = basis(:,3*(k-1)+[1,2,3],2).*obj.negativeEdges;
                basis(:,3*(k-1)+[1,2,3],3) = basis(:,3*(k-1)+[1,2,3],3).*obj.negativeEdges;
            end
            basis(:,3*(obj.degree+1)+1 : obj.nrBasisFunc,:) = ones(obj.mesh.nrElems,1).*obj.getElemRefBasis(obj.degree,x,obj.nrBasisFunc,obj.bubbles);

            % Piola transformation
            % J'*sigma*J
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,3);
            val(:,:,1) = basis(:,:,1).*obj.jacobian(:,1,1).^2 + 2*basis(:,:,3).*obj.jacobian(:,1,1).*obj.jacobian(:,2,1) + basis(:,:,2).*obj.jacobian(:,2,1).^2;
            val(:,:,2) = basis(:,:,1).*obj.jacobian(:,1,2).^2 + 2*basis(:,:,3).*obj.jacobian(:,1,2).*obj.jacobian(:,2,2) + basis(:,:,2).*obj.jacobian(:,2,2).^2;
            val(:,:,3) = obj.jacobian(:,1,2).*(basis(:,:,1).*obj.jacobian(:,1,1) + basis(:,:,3).*obj.jacobian(:,2,1)) ...
                + obj.jacobian(:,2,2).*(basis(:,:,2).*obj.jacobian(:,2,1) + basis(:,:,3).*obj.jacobian(:,1,1));
            val = (0.5./obj.mesh.area4e).^2.*val;
        end

        function val = getDivBasis(obj,~,x,~)
            %getDivBasis - Evaluate the divergence of the basis
            %   VAL = getDivBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   Piola-mapped divergences of the HHJ basis functions at the
            %   reference point PTS_REF (1-by-2) in all elements. VAL is
            %   nrElems-by-nrBasisFunc-by-2. PTS and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,2);
            val(:,1:3*(obj.degree+1),:) = ones(obj.mesh.nrElems,1).*obj.getEdgeDivRefBasis(obj.degree,x);
            for k = 2 : 2 : obj.degree+1
                val(:,3*(k-1)+[1,2,3],1) = val(:,3*(k-1)+[1,2,3],1).*obj.negativeEdges;
                val(:,3*(k-1)+[1,2,3],2) = val(:,3*(k-1)+[1,2,3],2).*obj.negativeEdges;
            end
            val(:,3*(obj.degree+1)+1 : obj.nrBasisFunc,:) = ones(obj.mesh.nrElems,1).*obj.getElemDivRefBasis(obj.degree,x,obj.nrBasisFunc,obj.bubbles);

            % Piola transformation
            val = matMul(val,obj.jacobian);
            val = (0.5./obj.mesh.area4e).^2.*val;
        end

        function val = getDivDivBasis(obj,~,x,~)
            %getDivDivBasis - Evaluate double divergences of the basis
            %   VAL = getDivDivBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   double divergences of the HHJ basis functions at the
            %   reference point PTS_REF (1-by-2) in all elements. VAL is
            %   nrElems-by-nrBasisFunc. PTS and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc);
            val(:,1:3*(obj.degree+1)) = ones(obj.mesh.nrElems,1).*obj.getEdgeDivDivRefBasis(obj.degree,x);
            for k = 2 : 2 : obj.degree+1
                val(:,3*(k-1)+[1,2,3]) = val(:,3*(k-1)+[1,2,3],1).*obj.negativeEdges;
            end
            val(:,3*(obj.degree+1)+1 : obj.nrBasisFunc) = ones(obj.mesh.nrElems,1).*obj.getElemDivDivRefBasis(obj.degree,x,obj.nrBasisFunc,obj.bubbles);

            % Piola transformation
            val = (0.5./obj.mesh.area4e).^2.*val;
        end

        %% FE values

        function val = projectBoundary(obj,ed4Bd,g,gDegree)
            %projectBoundary - Project data onto boundary edges
            %   VAL = projectBoundary(OBJ,ed4Bd,G) computes the L2
            %   projection of the normal-normal component (normal
            %   traction) n'*sigma*n given by the function handle G onto
            %   the edge degrees of freedom of the boundary edges ed4Bd,
            %   solving the local boundary mass matrix systems. G is
            %   called as G(X,NORMALS) or, if it takes four arguments, as
            %   G(PTS,PTS_REF,INDICES,NORMALS) with the outer unit normals
            %   NORMALS of the edges, and returns one scalar per point.
            %   VAL is the nrDoF-by-1 coefficient vector, zero away from
            %   the edges ed4Bd.
            %
            %   VAL = projectBoundary(OBJ,ed4Bd,G,gDegree) also specifies
            %   the polynomial degree gDegree of G, which sets the
            %   quadrature degree OBJ.degree+gDegree of the right-hand
            %   side; it defaults to OBJ.degree.
            if nargin < 4 || isempty(gDegree)
                gDegree = obj.degree;
            end
            normals = obj.mesh.normal4ed;
            normals = normals(ed4Bd,:);
            tangents = obj.mesh.tangent4ed;
            tangents = tangents(ed4Bd,:);
            M = integrate(@(pts,pts_ref,indices)obj.intBoundaryMatrix(pts,pts_ref,indices,obj,obj.mesh,normals,tangents),...
                obj.mesh.n4ed(ed4Bd,:),ed4Bd,obj.mesh,2*obj.degree);
            e4ed = obj.mesh.e4ed(ed4Bd,1);
            dofs = obj.dof4e(e4ed,:);
            [I,J] = localDoFtoGlobalDoF(dofs);
            M = sparse(I(:),J(:),M(:),obj.nrDoF,obj.nrDoF);
            b = integrate(@(pts,pts_ref,indices)obj.intBoundaryRHS(pts,pts_ref,indices,obj,obj.mesh,normals,tangents,g),...
                obj.mesh.n4ed(ed4Bd,:),ed4Bd,obj.mesh,obj.degree+gDegree);
            b = accumarray(dofs(:),b(:),[obj.nrDoF 1]);
            freeNodesBb = unique(obj.dof4ed(ed4Bd,:));
            val = zeros(obj.nrDoF,1);
            perm = symamd(M(freeNodesBb,freeNodesBb));
            freeNodesBb = freeNodesBb(perm);
            val(freeNodesBb) = M(freeNodesBb,freeNodesBb)\b(freeNodesBb);
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the number nrBasisFunc of basis
            %   functions, with 2*DEGREE more for the bubbles of HHJB,
            %   and the numbers nrBasis4n, nrBasis4ed, and nrBasis4e of
            %   nodal, edge, and interior degrees of freedom. The
            %   constructor calls it right before dofManager, so that
            %   subclasses may override it.
            k = obj.degree;
            obj.nrBasisFunc = max(3,3*(3*k+(k-1)*max(0,k-2)/2)) + 2*k*obj.bubbles;
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = obj.degree+1;
            obj.nrBasis4e = obj.nrBasisFunc-3*obj.nrBasis4ed;
        end

        %% Organisation of the degrees of freedom
        function obj = dofManager(obj)
            %dofManager - Number the degrees of freedom of the element
            %   OBJ = dofManager(OBJ) enumerates the edge and element
            %   degrees of freedom of the HHJ element with
            %   FE.dofManager, fixes them on the Dirichlet and Neumann
            %   boundaries, and records the polynomial degree of each
            %   local degree of freedom.
            obj = dofManager@FE(obj);
            obj.fixedNodes = [obj.dof4Db;obj.dof4Nb];
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
            % determine degree4DoF (useful for hierarchical fem), the
            % polynomial degree of the basis functions (elemFunctions):
            % p for the edge functions -Airy(u_(p+2)), i+j-2 for the
            % types I-III, j for type IV and DEGREE+1 for the bubbles
            obj.degree4DoF = zeros(1,size(obj.dof4e,2));
            for p = 0 : obj.degree
                obj.degree4DoF(3*p+(1:3)) = p;
            end
            indx = 3*(obj.degree+1)+1;
            for i = 2 : obj.degree + 1
                for j = 1 : obj.degree - i + 2
                    obj.degree4DoF(indx) = i+j-2;
                    indx = indx + 1;
                end
            end
            for i = 2 : obj.degree + 1
                for j = 1 : obj.degree - i + 2
                    obj.degree4DoF(indx) = i+j-2;
                    indx = indx + 1;
                end
            end
            for i = 2 : obj.degree + 1
                for j = 2 : obj.degree - i + 2
                    obj.degree4DoF(indx) = i+j-2;
                    indx = indx + 1;
                end
            end
            for j = 1 : obj.degree
                obj.degree4DoF(indx) = j;
                indx = indx + 1;
            end
            % divdiv-free bubbles of HHJB
            obj.degree4DoF(indx:end) = obj.degree+1;
        end

    end

    methods(Static)

        function val = intBoundaryMatrix(~,pts_ref,indices,fe,mesh,normals,~)
            %intBoundaryMatrix - Integrate the boundary mass matrix
            %   VAL = intBoundaryMatrix(PTS,PTS_REF,INDICES,FE,MESH,NORMALS,TANGENTS)
            %   evaluates all basis functions of FE in the element
            %   e4ed(:,1) of each edge INDICES at the quadrature point
            %   PTS_REF of the edge, takes their normal-normal components
            %   n'*sigma_i*n with the normals NORMALS (one row per edge)
            %   and returns their products, the integrand of the boundary
            %   mass matrix, of size
            %   numel(INDICES)-by-nrBasisFunc-by-nrBasisFunc. PTS and
            %   TANGENTS are ignored.
            basis = getEdgeValues4RefElemFunktion(pts_ref,@(x)fe.getBasis([],x,[]),mesh);
            basis = basis(indices,:,:,1);
            basis = matMul(basis,[normals(:,[1 2]).^2,2*prod(normals,2)]);
            val = matMul(basis,permute(basis,[1 3 2]));
        end

        function val = intBoundaryRHS(pts,pts_ref,indices,fe,mesh,normals,~,g)
            %intBoundaryRHS - Integrate the boundary load vector
            %   VAL = intBoundaryRHS(PTS,PTS_REF,INDICES,FE,MESH,NORMALS,TANGENTS,G)
            %   evaluates all basis functions of FE in the element
            %   e4ed(:,1) of each edge INDICES at the quadrature point
            %   PTS_REF of the edge, takes their normal-normal components
            %   n'*sigma_i*n with the normals NORMALS (one row per edge)
            %   and multiplies them by the scalar data G(PTS,NORMALS), or
            %   G(PTS,PTS_REF,INDICES,NORMALS) if G takes four arguments.
            %   VAL is numel(INDICES)-by-nrBasisFunc, the integrand of
            %   the boundary load vector. TANGENTS is ignored.
            basis = getEdgeValues4RefElemFunktion(pts_ref,@(x)fe.getBasis([],x,[]),mesh);
            basis = basis(indices,:,:,1);
            basis = matMul(basis,[normals(:,[1 2]).^2,2*prod(normals,2)]);
            if nargin(g) == 4
                curG = g(pts,pts_ref,indices,normals);
            else
                curG = g(pts,normals);
            end
            val = curG.*basis;
        end

        % The basis functions are built from the scaled integrated Jacobi
        % factors u_i = L_i(lambda_b,lambda_a+lambda_b) and
        % v_j = L_j^alpha(lambda_3,1) of BarycentricJacobi with the Airy
        % function Airy(w) = [w_yy, w_xx, -w_xy] in the component order
        % [sigma_xx, sigma_yy, sigma_xy].

        %% Edge Basis helper functions
        function val = getEdgeRefBasis(p,x)
            %getEdgeRefBasis - Evaluate the edge reference basis
            %   VAL = getEdgeRefBasis(P,X) evaluates the integrated
            %   Jacobi-based HHJ edge basis functions -Airy(u_i),
            %   i = 2,...,P+2, of the edges [a b] on the reference
            %   triangle at the point X (1-by-2). VAL is
            %   1-by-3*(P+1)-by-3, with column 3*(i-2)+e for edge e.
            edgeNodes = [1 2; 2 3; 3 1];
            lambda = BarycentricJacobi.coordinates(x);
            val = zeros(1,3*(p+1),3);
            for e = 1 : 3
                U = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(edgeNodes(e,2),edgeNodes(e,:),2,2),p+2,0);
                val(1,e:3:end,:) = -U(3:end,[5 4 6]).*[1 1 -1];
            end
        end

        function val = getEdgeDivRefBasis(p,x)
            %getEdgeDivRefBasis - Evaluate the edge div reference basis
            %   VAL = getEdgeDivRefBasis(P,X) returns zeros since all HHJ
            %   edge basis functions are divergence free on the reference
            %   triangle. VAL is size(X,1)-by-3*(P+1)-by-2.
            val = zeros(size(x,1),3*(p+1),2);
        end

        function val = getEdgeDivDivRefBasis(p,x)
            %getEdgeDivDivRefBasis - Evaluate the edge divdiv ref basis
            %   VAL = getEdgeDivDivRefBasis(P,X) returns zeros since all
            %   HHJ edge basis functions have vanishing double divergence
            %   on the reference triangle. VAL is size(X,1)-by-3*(P+1).
            val = zeros(size(x,1),3*(p+1));
        end

        %% Elem Basis helper functions
        function val = getElemRefBasis(p,x,nrBasisFunc,bubbles)
            %getElemRefBasis - Evaluate the element reference basis
            %   VAL = getElemRefBasis(P,X,nrBasisFunc,BUBBLES) evaluates
            %   the interior HHJ basis functions of degree P on the
            %   reference triangle at the point X (1-by-2), including the
            %   divdiv-free bubbles of HHJB if BUBBLES is true. VAL is
            %   1-by-(nrBasisFunc-3*(P+1))-by-3.
            val = HHJk.elemFunctions(p,x,nrBasisFunc,bubbles,0);
        end

        function val = getElemDivRefBasis(p,x,nrBasisFunc,bubbles)
            %getElemDivRefBasis - Evaluate the element div ref basis
            %   VAL = getElemDivRefBasis(P,X,nrBasisFunc,BUBBLES)
            %   evaluates the divergences of the interior HHJ basis
            %   functions of degree P, with the bubbles of HHJB if BUBBLES
            %   is true, at the point X (1-by-2) of the reference triangle.
            %   VAL is 1-by-(nrBasisFunc-3*(P+1))-by-2.
            val = HHJk.elemFunctions(p,x,nrBasisFunc,bubbles,1);
        end

        function val = getElemDivDivRefBasis(p,x,nrBasisFunc,bubbles)
            %getElemDivDivRefBasis - Evaluate the element divdiv ref basis
            %   VAL = getElemDivDivRefBasis(P,X,nrBasisFunc,BUBBLES)
            %   evaluates the double divergences of the interior HHJ basis
            %   functions of degree P, with the bubbles of HHJB if BUBBLES
            %   is true, at the point X (1-by-2) of the reference triangle.
            %   VAL is 1-by-(nrBasisFunc-3*(P+1)).
            val = HHJk.elemFunctions(p,x,nrBasisFunc,bubbles,2);
        end

        %% P1 basis functions
        function val = p1RefBasis(x)
            %p1RefBasis - Evaluate the P1 reference basis
            %   VAL = p1RefBasis(X) evaluates the three linear nodal
            %   basis functions of the reference triangle at the points X
            %   (N-by-2). VAL is N-by-3.
            val = zeros(size(x,1),3);
            val(:,2) = x(:,1);
            val(:,3) = x(:,2);
            val(:,1) = 1 - x(:,1) - x(:,2);
        end

    end

    methods (Static, Access = private)

        function val = elemFunctions(p,x,nrBasisFunc,bubbles,m)
            %elemFunctions - Evaluate the cell basis or its divergences
            %   VAL = elemFunctions(P,X,nrBasisFunc,BUBBLES,M) evaluates
            %   the cell basis functions at the point X (1-by-2) for M = 0,
            %   their divergences for M = 1, or their double divergences
            %   for M = 2, with u_i = L_i(lambda_2,lambda_1+lambda_2) and
            %   v_j = L_j^(2i-1)(lambda_3,1), i = 2,...,P+1: -Airy(u_i v_j),
            %   j = 1,...,P-i+2 (type I, div-free), -Airy(u_i v_j) +
            %   2 cof(grad(u_i) grad(v_j)' + grad(v_j) grad(u_i)'),
            %   j = 1,...,P-i+2 (type II), u_i Airy(v_j) - v_j Airy(u_i),
            %   j = 2,...,P-i+2 (type III), d_y(v_j) [x, 0, (y-1)/2] with
            %   alpha = 1, j = 1,...,P (type IV), and for BUBBLES the
            %   divdiv-free sym(curl(w) grad(lambda_3)') and
            %   -sym(curl(w) grad(lambda_2)') with w = u_i v_(P-i+2)
            %   (type V). VAL is 1-by-(nrBasisFunc-3*(P+1))-by-3 for
            %   M = 0, 1-by-(nrBasisFunc-3*(P+1))-by-2 for M = 1 and
            %   1-by-(nrBasisFunc-3*(P+1)) for M = 2.
            lambda = BarycentricJacobi.coordinates(x);
            M = 2+m;
            termsM = BarycentricJacobi.leibnizTerms(BarycentricJacobi.multiIndices(2,M),M);
            terms = BarycentricJacobi.leibnizTerms(BarycentricJacobi.multiIndices(2,m),m);
            U = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(2,[1 2],M,2),p+1,0);
            chainV = BarycentricJacobi.chainRule(3,1:3,M,2);
            d = @(F,alpha) HHJk.shift(F,alpha,M,m);
            % Leibniz rule for the product of the derivatives alpha and beta
            uv = @(u,alpha,V,beta) BarycentricJacobi.product(d(u,alpha),d(V,beta),terms);
            S = cell(5,1);
            for i = 2 : p+1
                V = BarycentricJacobi.factor(lambda,chainV,p-i+2,2*i-1);
                u = U(i+1,:);
                W = BarycentricJacobi.product(u,V,termsM);
                airy = {-d(W,[0 2]), -d(W,[2 0]), d(W,[1 1])};
                S{1} = [S{1}; HHJk.rows(airy,2:size(V,1))];
                cof = {4*uv(u,[0 1],V,[0 1]), 4*uv(u,[1 0],V,[1 0]), -2*(uv(u,[1 0],V,[0 1])+uv(u,[0 1],V,[1 0]))};
                S{2} = [S{2}; HHJk.rows(cellfun(@plus,airy,cof,'UniformOutput',false),2:size(V,1))];
                diff = {uv(u,[0 0],V,[0 2])-uv(u,[0 2],V,[0 0]), uv(u,[0 0],V,[2 0])-uv(u,[2 0],V,[0 0]), ...
                    uv(u,[1 1],V,[0 0])-uv(u,[0 0],V,[1 1])};
                S{3} = [S{3}; HHJk.rows(diff,3:size(V,1))];
                if bubbles
                    n = size(V,1);
                    S{5} = [S{5}; HHJk.rows({0*d(W,[0 0]), -d(W,[1 0]), d(W,[0 1])/2},n); ...
                        HHJk.rows({-d(W,[0 1]), 0*d(W,[0 0]), d(W,[1 0])/2},n)];
                end
            end
            % type IV
            V = BarycentricJacobi.factor(lambda,chainV,p,1);
            n = (m+1)*(m+2)/2;
            X = [x(1) 1 0 0 0 0];
            Y = [(x(2)-1)/2 0 1/2 0 0 0];
            S{4} = HHJk.rows({BarycentricJacobi.product(X(1:n),d(V,[0 1]),terms), 0*d(V,[0 0]), ...
                BarycentricJacobi.product(Y(1:n),d(V,[0 1]),terms)},2:size(V,1));
            S = [zeros(0,3*n); vertcat(S{:})];
            S11 = S(:,1:n); S22 = S(:,n+(1:n)); S12 = S(:,2*n+(1:n));
            switch m
                case 0
                    val = [S11, S22, S12];
                case 1
                    val = [S11(:,2)+S12(:,3), S12(:,2)+S22(:,3)];
                otherwise
                    val = S11(:,4)+2*S12(:,6)+S22(:,5);
            end
            val = reshape(val,1,max(0,nrBasisFunc-3*(p+1)),size(val,2));
        end

        function G = shift(F,alpha,M,m)
            %SHIFT - Derivatives of a derivative
            %   G = SHIFT(F,ALPHA,M,M0) returns the derivatives up to
            %   order M0 of the derivative d^ALPHA f, where the columns of
            %   F hold the derivatives of f up to order M as in
            %   BarycentricJacobi.multiIndices(2,M).
            D = BarycentricJacobi.multiIndices(2,M);
            [~,cols] = ismember(BarycentricJacobi.multiIndices(2,m)+alpha,D,'rows');
            G = F(:,cols);
        end

        function S = rows(components,rows)
            %ROWS - Stack the derivatives of the three components
            %   S = ROWS(COMPONENTS,ROWS) returns [S11, S22, S12] for the
            %   rows ROWS of the derivatives of the three components.
            S = [components{1}(rows,:), components{2}(rows,:), components{3}(rows,:)];
        end
    end

end