classdef HPk3d < FE
    %HPk3d - Hierarchical high-order H1 element on tetrahedra
    %   OBJ = HPk3d(DEGREE,MESH) creates the tetrahedral hierarchical
    %   high-order scalar H1-conforming finite element of polynomial
    %   degree DEGREE >= 1 on the tetrahedral mesh MESH with nodal,
    %   edge, face and cell bubble degrees of freedom. DEGREE is one
    %   scalar for all elements. The basis functions are products of
    %   scaled integrated Jacobi polynomials in the barycentric
    %   coordinates (BarycentricJacobi), oriented by the global node
    %   numbers of the edges and faces (sign4e, perm4e).
    %
    %   HPk3d functions:
    %       getBasis            - Evaluate the basis functions
    %       getGradBasis        - Evaluate the gradients of the basis
    %       getD2Basis          - Evaluate the Hessians of the basis
    %       getD3Basis          - Evaluate third derivatives of the basis
    %       getD4Basis          - Evaluate fourth derivatives of the basis
    %       getEdgeRefBasis     - Evaluate the edge reference basis
    %       getEdgeGradRefBasis - Evaluate the edge gradient basis
    %       getEdgeD2RefBasis   - Evaluate the edge D2 reference basis
    %       getEdgeD3RefBasis   - Evaluate the edge D3 reference basis
    %       getEdgeD4RefBasis   - Evaluate the edge D4 reference basis
    %       getFaceRefBasis     - Evaluate the face reference basis
    %       getFaceGradRefBasis - Evaluate the face gradient basis
    %       getFaceD2RefBasis   - Evaluate the face D2 reference basis
    %       getFaceD3RefBasis   - Evaluate the face D3 reference basis
    %       getFaceD4RefBasis   - Evaluate the face D4 reference basis
    %       getElemRefBasis     - Evaluate the element reference basis
    %       getElemGradRefBasis - Evaluate the element gradient basis
    %       getElemD2RefBasis   - Evaluate the element D2 reference basis
    %       getElemD3RefBasis   - Evaluate the element D3 reference basis
    %       getElemD4RefBasis   - Evaluate the element D4 reference basis
    %       p1RefBasis          - Evaluate the P1 reference basis
    %       p1GradRefBasis      - Evaluate the gradients of the P1 basis
    %
    %   HPk3d properties:
    %       nrDoF             - Total number of degrees of freedom
    %       dof4e             - Local-to-global dofs per element
    %       freeNodes         - Indices of the free degrees of freedom
    %       fixedNodes        - Indices of the fixed degrees of freedom
    %       dof4Db            - Dofs on the Dirichlet boundary
    %       degree            - Polynomial degree of the element
    %       nrBasis4n         - Number of nodal basis functions
    %       nrBasis4ed        - Number of edge basis functions
    %       nrBasis4f         - Number of face basis functions
    %       nrBasis4e         - Number of element bubble functions
    %       nrBasisFunc       - Number of basis functions per element
    %       sign4e            - Orientation signs of the edges and faces
    %       perm4e            - Face permutation index per element
    %       facePolyIndices   - Index pairs of the face polynomials
    %       transformation    - Inverse Jacobians of the element maps
    %       mesh              - Mesh on which the element is defined
    %       basisCoefficients - Orientation coefficients of the basis
    %       isRefElemBased    - Basis lives on the reference element
    %
    %   See also HPk, Pk3d, DGDk3d, Mesh3d, BarycentricJacobi, FE

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = public)
        nrDoF             % total number of degrees of freedom
        dof4e             % local-to-global dofs per element
        freeNodes         % indices of the free degrees of freedom
        fixedNodes        % indices of the fixed degrees of freedom
        dof4Db            % dofs on the Dirichlet boundary
        degree            % polynomial degree of the element
        nrBasis4n         % number of nodal basis functions
        nrBasis4ed        % number of edge basis functions
        nrBasis4f         % number of face basis functions
        nrBasis4e         % number of element bubble functions
        nrBasisFunc       % number of basis functions per element
        sign4e            % orientation signs of the edges and faces
        perm4e            % face permutation index per element
        facePolyIndices   % index pairs of the face polynomials
        transformation    % inverse Jacobians of the element maps
        mesh              % mesh on which the element is defined
        basisCoefficients % orientation coefficients of the basis
        isRefElemBased = true % basis lives on the reference element
    end

    methods

        function obj = HPk3d(degree,mesh)
            %HPk3d - Create a tetrahedral hierarchical Pk element object
            if degree < 1
                error('fem4m:HPk3d:invalidDegree','HPk3d degree must be at least 1');
            end
            obj.degree = degree;
            obj = obj.dofCounts();
            obj.nrBasisFunc = 4*obj.nrBasis4n + 6*obj.nrBasis4ed + 4*obj.nrBasis4f + obj.nrBasis4e;
            obj = obj.dofManager(mesh);
            obj = obj.getBasisCoefficients(mesh);
            obj.transformation = mesh.transformation();
            obj.mesh = mesh;
        end

        function val = get.basisCoefficients(obj)
            %basisCoefficients - Build the orientation coefficients from sign4e
            nrElems = size(obj.dof4e,1);
            val = ones(nrElems,obj.nrBasisFunc);
            for k = 2 : 2 : obj.nrBasis4ed
                val(:,4*obj.nrBasis4n+6*(k-1)+(1:6)) = obj.sign4e(:,1:6);
            end
            for k = 1 : obj.nrBasis4f
                for curFace = 1 : 4
                    for perm = 1 : 3
                        indx = obj.perm4e(:,curFace)==perm;
                        if mod(obj.facePolyIndices(k,1),2)==1
                            val(indx,4*obj.nrBasis4n+6*obj.nrBasis4ed + 4*(k-1)+curFace) = obj.sign4e(indx,6+curFace);
                        end
                    end
                end
            end
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the nodal,
            %   edge, face and cell bubble basis functions, oriented per
            %   element, in all elements at the reference point PTS_REF
            %   (1-by-3). VAL is nrElems-by-nrBasisFunc. PTS and PARTS
            %   are ignored.
            val = obj.refBasis(x,0);
        end

        function val = getGradBasis(obj,~,x,~)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   gradients of the basis functions at the reference point
            %   PTS_REF (1-by-3) and maps them to the physical elements
            %   with transformation. VAL is nrElems-by-nrBasisFunc-by-3.
            %   PTS and PARTS are ignored.
            val = FE.mapDerivatives(obj.transformation,obj.refBasis(x,1),1);
        end

        function val = getD2Basis(obj,~,x,~)
            %getD2Basis - Evaluate the Hessians of the basis
            %   VAL = getD2Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   second derivatives (xx, yy, zz, xy, xz, yz) of the basis
            %   functions at the reference point PTS_REF (1-by-3),
            %   mapped to the physical elements. PTS and PARTS are
            %   ignored.
            val = FE.mapDerivatives(obj.transformation,obj.refBasis(x,2),2);
        end

        function val = getD3Basis(obj,~,x,~)
            %getD3Basis - Evaluate third derivatives of the basis
            %   VAL = getD3Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   third derivatives (xxx, yyy, zzz, xxy, xxz, xyy, yyz,
            %   xzz, yzz, xyz) of the basis functions at the reference
            %   point PTS_REF (1-by-3), mapped to the physical elements.
            %   PTS and PARTS are ignored.
            val = FE.mapDerivatives(obj.transformation,obj.refBasis(x,3),3);
        end

        function val = getD4Basis(obj,~,x,~)
            %getD4Basis - Evaluate fourth derivatives of the basis
            %   VAL = getD4Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   fourth derivatives (xxxx, yyyy, zzzz, xxxy, xxxz, xyyy,
            %   yyyz, xzzz, yzzz, xxyy, xxzz, yyzz, xxyz, xyyz, xyzz) of
            %   the basis functions at the reference point PTS_REF
            %   (1-by-3), mapped to the physical elements. PTS and PARTS
            %   are ignored.
            val = FE.mapDerivatives(obj.transformation,obj.refBasis(x,4),4);
        end

    end

    methods (Access = private)

        function val = refBasis(obj,x,order)
            %refBasis - Evaluate the oriented reference basis derivatives
            %   VAL = refBasis(OBJ,X,ORDER) evaluates the derivatives of
            %   order ORDER = 0,...,4 of the nodal, edge, face and cell
            %   bubble basis functions on the reference tetrahedron at the
            %   point X (1-by-3), with the face polynomials of the
            %   permutations perm4e and the orientation signs
            %   basisCoefficients of the elements. VAL has one row per
            %   element, one column per function, and the components
            %   of Monomials.multiIndices(3,ORDER), one for ORDER = 0.
            op = {'','Grad','D2','D3','D4'};
            op = op{order+1};
            nrElems = size(obj.dof4e,1);
            val = zeros(nrElems,obj.nrBasisFunc,nchoosek(order+2,2));
            offset = 4*obj.nrBasis4n;
            if order == 0
                val(:,1:offset,:) = ones(nrElems,1).*permute(obj.p1RefBasis(x),[3 1 2]);
            elseif order == 1
                val(:,1:offset,:) = ones(nrElems,1).*permute(obj.p1GradRefBasis(),[3 1 2]);
            end
            val(:,offset+(1:6*obj.nrBasis4ed),:) = ones(nrElems,1).*permute( ...
                obj.(['getEdge',op,'RefBasis'])(obj.degree,x,obj.nrBasis4ed),[3 1 2]);
            offset = offset + 6*obj.nrBasis4ed;
            faceValues = obj.(['getFace',op,'RefBasis'])(obj.degree,x,obj.nrBasis4f);
            for k = 1 : obj.nrBasis4f
                for curFace = 1 : 4
                    for perm = 1 : 3
                        indx = obj.perm4e(:,curFace)==perm;
                        val(indx,offset+4*(k-1)+curFace,:) = ones(nnz(indx),1).*faceValues(curFace,k,:,perm);
                    end
                end
            end
            offset = offset + 4*obj.nrBasis4f;
            val(:,offset+1:end,:) = ones(nrElems,1).*permute( ...
                obj.(['getElem',op,'RefBasis'])(obj.degree,x,obj.nrBasis4e),[3 1 2]);
            val = obj.basisCoefficients.*val;
        end

        function obj = getBasisCoefficients(obj,mesh)
            %getBasisCoefficients - Set the orientation signs and face permutations
            obj.sign4e = ones(mesh.nrElems,10);
            obj.perm4e = zeros(mesh.nrElems,4);
            n4e = mesh.n4e;
            % edge orientation
            edgeNodes = [ 1 2; 2 3; 3 1; 1 4; 2 4; 3 4];
            for j = 1 : 6
                indx = n4e(:,edgeNodes(j,1)) > n4e(:,edgeNodes(j,2));
                obj.sign4e(indx,j) = -1;
            end
            % face orientation
            faceNodes = [1 2 3; 1 2 4; 4 2 3; 1 3 4];
            permutations = [1 2 3; 2 3 1; 3 1 2];
            negativepermutations = [2 1 3; 3 2 1; 1 3 2];
            for j = 1 : 4
                [~,sortIdx] = sort(n4e(:,faceNodes(j,:)),2);
                [~,perm] = ismember(sortIdx,permutations,'rows');
                [~,negperm] = ismember(sortIdx,negativepermutations,'rows');
                obj.sign4e(negperm>0,6+j) = -1;
                obj.perm4e(:,j) = perm + negperm;
            end
            % save indices i j for face polynomials
            obj.facePolyIndices = zeros(obj.nrBasis4f,2);
            indx = 1;
            p = obj.degree;
            for i = 2 : p - 1
                for j = 1 : p - i
                    obj.facePolyIndices(indx,:) = [i,j];
                    indx = indx + 1;
                end
            end
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, nrBasis4f, and nrBasis4e of nodal, edge,
            %   face, and interior degrees of freedom. The constructor
            %   calls it right before dofManager, so that subclasses
            %   may override it.
            obj.nrBasis4n = 1;
            obj.nrBasis4ed = obj.degree-1;
            obj.nrBasis4f = (obj.degree-1)*max(0,obj.degree-2)/2;
            obj.nrBasis4e = (obj.degree-1)*max(0,obj.degree-2)*max(0,obj.degree-3)/6;
        end

        function obj = dofManager(obj,mesh)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ,MESH) enumerates the nodal, edge,
            %   face, and interior degrees of freedom with FE.dofManager
            %   and fixes them on the Dirichlet boundary.
            obj.mesh = mesh;
            obj = dofManager@FE(obj);
            obj.fixedNodes = obj.dof4Db;
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
        end

    end

    methods(Static)

        % ordering of the reference tetrahedron
        % c4n: [ 0 0 0; 1 0 0; 0 1 0; 0 0 1];
        % edges: [1 2], [2 3], [3 1], [1 4], [2 4], [3 4]

        %% Edge Basis helper functions
        function val = getEdgeRefBasis(p,x,nrBasisFunc4ed)
            %getEdgeRefBasis - Evaluate the edge reference basis
            %   VAL = getEdgeRefBasis(P,X,nrBasisFunc4ed) evaluates the
            %   nrBasisFunc4ed = P-1 edge functions of degree i = 2,...,P
            %   of the six edges of the reference tetrahedron at the
            %   point X (1-by-3). VAL is 6*nrBasisFunc4ed-by-1, row
            %   6*(i-2)+e for edge e.
            val = HPk3d.edgeDerivatives(p,x,nrBasisFunc4ed,0);
        end

        function val = getEdgeGradRefBasis(p,x,nrBasisFunc4ed)
            %getEdgeGradRefBasis - Evaluate the edge gradient basis
            %   VAL = getEdgeGradRefBasis(P,X,nrBasisFunc4ed) evaluates
            %   the gradients of the edge functions of getEdgeRefBasis at
            %   the point X (1-by-3), one row per function and one column
            %   per derivative.
            val = HPk3d.edgeDerivatives(p,x,nrBasisFunc4ed,1);
        end

        function val = getEdgeD2RefBasis(p,x,nrBasisFunc4ed)
            %getEdgeD2RefBasis - Evaluate the edge D2 reference basis
            %   VAL = getEdgeD2RefBasis(P,X,nrBasisFunc4ed) evaluates the
            %   second derivatives (xx, yy, zz, xy, xz, yz) of the edge
            %   functions of getEdgeRefBasis at the point X (1-by-3).
            val = HPk3d.edgeDerivatives(p,x,nrBasisFunc4ed,2);
        end

        function val = getEdgeD3RefBasis(p,x,nrBasisFunc4ed)
            %getEdgeD3RefBasis - Evaluate the edge D3 reference basis
            %   VAL = getEdgeD3RefBasis(P,X,nrBasisFunc4ed) evaluates the
            %   third derivatives of the edge functions of
            %   getEdgeRefBasis at the point X (1-by-3), in the order of
            %   getD3Basis.
            val = HPk3d.edgeDerivatives(p,x,nrBasisFunc4ed,3);
        end

        function val = getEdgeD4RefBasis(p,x,nrBasisFunc4ed)
            %getEdgeD4RefBasis - Evaluate the edge D4 reference basis
            %   VAL = getEdgeD4RefBasis(P,X,nrBasisFunc4ed) evaluates the
            %   fourth derivatives of the edge functions of
            %   getEdgeRefBasis at the point X (1-by-3), in the order of
            %   getD4Basis.
            val = HPk3d.edgeDerivatives(p,x,nrBasisFunc4ed,4);
        end

        %% Face Basis helper functions
        function val = getFaceRefBasis(p,x,nrBasisFunc4f)
            %getFaceRefBasis - Evaluate the face reference basis
            %   VAL = getFaceRefBasis(P,X,nrBasisFunc4f) evaluates the
            %   nrBasisFunc4f = (P-1)*(P-2)/2 face bubbles of degree up
            %   to P of the four faces for the three cyclic permutations
            %   of their nodes at the point X (1-by-3). VAL is
            %   4-by-nrBasisFunc4f-by-1-by-3.
            val = HPk3d.faceDerivatives(p,x,nrBasisFunc4f,0);
        end

        function val = getFaceGradRefBasis(p,x,nrBasisFunc4f)
            %getFaceGradRefBasis - Evaluate the face gradient basis
            %   VAL = getFaceGradRefBasis(P,X,nrBasisFunc4f) evaluates the
            %   gradients of the face bubbles of getFaceRefBasis at the
            %   point X (1-by-3). VAL is 4-by-nrBasisFunc4f-by-3-by-3.
            val = HPk3d.faceDerivatives(p,x,nrBasisFunc4f,1);
        end

        function val = getFaceD2RefBasis(p,x,nrBasisFunc4f)
            %getFaceD2RefBasis - Evaluate the face D2 reference basis
            %   VAL = getFaceD2RefBasis(P,X,nrBasisFunc4f) evaluates the
            %   second derivatives (xx, yy, zz, xy, xz, yz) of the face
            %   bubbles of getFaceRefBasis at the point X (1-by-3).
            val = HPk3d.faceDerivatives(p,x,nrBasisFunc4f,2);
        end

        function val = getFaceD3RefBasis(p,x,nrBasisFunc4f)
            %getFaceD3RefBasis - Evaluate the face D3 reference basis
            %   VAL = getFaceD3RefBasis(P,X,nrBasisFunc4f) evaluates the
            %   third derivatives of the face bubbles of getFaceRefBasis
            %   at the point X (1-by-3), in the order of getD3Basis.
            val = HPk3d.faceDerivatives(p,x,nrBasisFunc4f,3);
        end

        function val = getFaceD4RefBasis(p,x,nrBasisFunc4f)
            %getFaceD4RefBasis - Evaluate the face D4 reference basis
            %   VAL = getFaceD4RefBasis(P,X,nrBasisFunc4f) evaluates the
            %   fourth derivatives of the face bubbles of getFaceRefBasis
            %   at the point X (1-by-3), in the order of getD4Basis.
            val = HPk3d.faceDerivatives(p,x,nrBasisFunc4f,4);
        end

        %% Elem Basis helper functions
        function val = getElemRefBasis(p,x,nrBasisFunc4e)
            %getElemRefBasis - Evaluate the element reference basis
            %   VAL = getElemRefBasis(P,X,nrBasisFunc4e) evaluates the
            %   nrBasisFunc4e = (P-1)*(P-2)*(P-3)/6 cell bubbles of
            %   degree up to P at the point X (1-by-3). VAL is
            %   nrBasisFunc4e-by-1 (zero rows for P < 4).
            val = HPk3d.elemDerivatives(p,x,nrBasisFunc4e,0);
        end

        function val = getElemGradRefBasis(p,x,nrBasisFunc4e)
            %getElemGradRefBasis - Evaluate the element gradient basis
            %   VAL = getElemGradRefBasis(P,X,nrBasisFunc4e) evaluates the
            %   gradients of the cell bubbles of getElemRefBasis at the
            %   point X (1-by-3). VAL is nrBasisFunc4e-by-3.
            val = HPk3d.elemDerivatives(p,x,nrBasisFunc4e,1);
        end

        function val = getElemD2RefBasis(p,x,nrBasisFunc4e)
            %getElemD2RefBasis - Evaluate the element D2 reference basis
            %   VAL = getElemD2RefBasis(P,X,nrBasisFunc4e) evaluates the
            %   second derivatives (xx, yy, zz, xy, xz, yz) of the cell
            %   bubbles of getElemRefBasis at the point X (1-by-3).
            val = HPk3d.elemDerivatives(p,x,nrBasisFunc4e,2);
        end

        function val = getElemD3RefBasis(p,x,nrBasisFunc4e)
            %getElemD3RefBasis - Evaluate the element D3 reference basis
            %   VAL = getElemD3RefBasis(P,X,nrBasisFunc4e) evaluates the
            %   third derivatives (xxx, yyy, zzz, xxy, xxz, xyy, yyz, xzz,
            %   yzz, xyz) of the cell bubbles of getElemRefBasis at the
            %   point X (1-by-3).
            val = HPk3d.elemDerivatives(p,x,nrBasisFunc4e,3);
        end

        function val = getElemD4RefBasis(p,x,nrBasisFunc4e)
            %getElemD4RefBasis - Evaluate the element D4 reference basis
            %   VAL = getElemD4RefBasis(P,X,nrBasisFunc4e) evaluates the
            %   fourth derivatives (xxxx, yyyy, zzzz, xxxy, xxxz, xyyy,
            %   yyyz, xzzz, yzzz, xxyy, xxzz, yyzz, xxyz, xyyz, xyzz) of
            %   the cell bubbles of getElemRefBasis at the point X
            %   (1-by-3).
            val = HPk3d.elemDerivatives(p,x,nrBasisFunc4e,4);
        end

        %% P1 basis functions
        function val = p1RefBasis(x)
            %p1RefBasis - Evaluate the P1 reference basis
            %   VAL = p1RefBasis(X) evaluates the four linear nodal basis
            %   functions [1-x-y-z; x; y; z] of the reference tetrahedron
            %   at the point X (1-by-3). VAL is 4-by-1.
            val = zeros(4,1);
            val(1) = 1 - x(1) - x(2) - x(3);
            val(2) = x(1);
            val(3) = x(2);
            val(4) = x(3);
        end

        function val = p1GradRefBasis()
            %p1GradRefBasis - Evaluate the gradients of the P1 basis
            %   VAL = p1GradRefBasis() returns the constant gradients of
            %   the four linear nodal basis functions of the reference
            %   tetrahedron, one row per function (4-by-3).
            val = zeros(4,3);
            val(1,:) = -1;
            val(2,1) = 1;
            val(3,2) = 1;
            val(4,3) = 1;
        end
    end

    methods (Static, Access = private)

        % The basis functions are products of the scaled integrated Jacobi
        % factors L_j(lambda_b,sum(lambda_nodes)) of BarycentricJacobi.

        function val = edgeDerivatives(p,x,nrBasisFunc4ed,order)
            %edgeDerivatives - Evaluate derivatives of the edge basis
            %   VAL = edgeDerivatives(P,X,nrBasisFunc4ed,ORDER) evaluates
            %   the derivatives of order ORDER of the edge
            %   functions L_i(lambda_b,lambda_a+lambda_b), i = 2,...,P,
            %   of the edges [a b] at the point X, one row per function
            %   and one column per component of
            %   Monomials.multiIndices(3,ORDER).
            edgeNodes = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4];
            lambda = BarycentricJacobi.coordinates(x);
            cols = BarycentricJacobi.columns(3,order);
            val = zeros(6*nrBasisFunc4ed,numel(cols));
            for e = 1 : 6
                F = BarycentricJacobi.factor(lambda, ...
                    BarycentricJacobi.chainRule(edgeNodes(e,2),edgeNodes(e,:),order,3),p,0);
                val(e:6:end,:) = F(3:end,cols);
            end
        end

        function val = faceDerivatives(p,x,nrBasisFunc4f,order)
            %faceDerivatives - Evaluate derivatives of the face basis
            %   VAL = faceDerivatives(P,X,nrBasisFunc4f,ORDER) evaluates
            %   the derivatives of order ORDER of the face bubble
            %   functions L_i(lambda_b,lambda_a+lambda_b)
            %   L_j^(2i-1)(lambda_c,lambda_a+lambda_b+lambda_c),
            %   i = 2,...,P-1, j = 1,...,P-i, of the four faces with the
            %   three cyclic permutations [a b c] of their nodes at the
            %   point X. VAL has the size 4 x nrBasisFunc4f x
            %   size(Monomials.multiIndices(3,ORDER),1) x 3.
            faceNodes = [1 2 3; 1 2 4; 4 2 3; 1 3 4];
            lambda = BarycentricJacobi.coordinates(x);
            D = BarycentricJacobi.multiIndices(3,order);
            cols = BarycentricJacobi.columns(3,order);
            terms = BarycentricJacobi.leibnizTerms(D(cols,:),order);
            val = zeros(4,nrBasisFunc4f,numel(cols),3);
            for f = 1 : 4
                for perm = 1 : 3
                    nodes = faceNodes(f,circshift(1:3,1-perm));
                    F = BarycentricJacobi.factor(lambda, ...
                        BarycentricJacobi.chainRule(nodes(2),nodes(1:2),order,3),p-1,0);
                    chain = BarycentricJacobi.chainRule(nodes(3),nodes,order,3);
                    indx = 1;
                    for i = 2 : p-1
                        G = BarycentricJacobi.factor(lambda,chain,p-i,2*i-1);
                        val(f,indx+(0:p-i-1),:,perm) = BarycentricJacobi.product(F(i+1,:),G(2:end,:),terms);
                        indx = indx + p-i;
                    end
                end
            end
        end

        function val = elemDerivatives(p,x,nrBasisFunc4e,order)
            %elemDerivatives - Evaluate derivatives of the cell basis
            %   VAL = elemDerivatives(P,X,nrBasisFunc4e,ORDER) evaluates
            %   the derivatives of order ORDER of the cell bubble
            %   functions L_i(lambda_2,lambda_1+lambda_2)
            %   L_j^(2i-1)(lambda_3,lambda_1+lambda_2+lambda_3)
            %   L_k^(2i+2j-2)(lambda_4,1), i = 2,...,P-2, j = 1,...,P-i-1,
            %   k = 1,...,P-i-j, at the point X, one row per function and
            %   one column per component of Monomials.multiIndices(3,ORDER).
            D = BarycentricJacobi.multiIndices(3,order);
            cols = BarycentricJacobi.columns(3,order);
            val = zeros(nrBasisFunc4e,numel(cols));
            if p < 4
                return;
            end
            lambda = BarycentricJacobi.coordinates(x);
            termsAll = BarycentricJacobi.leibnizTerms(D,order);
            terms = BarycentricJacobi.leibnizTerms(D(cols,:),order);
            F = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(2,[1 2],order,3),p-2,0);
            chainG = BarycentricJacobi.chainRule(3,[1 2 3],order,3);
            chainH = BarycentricJacobi.chainRule(4,1:4,order,3);
            indx = 1;
            for i = 2 : p-2
                FG = BarycentricJacobi.product(F(i+1,:),BarycentricJacobi.factor(lambda,chainG,p-i-1,2*i-1),termsAll);
                for j = 1 : p-i-1
                    H = BarycentricJacobi.factor(lambda,chainH,p-i-j,2*i+2*j-2);
                    val(indx+(0:p-i-j-1),:) = BarycentricJacobi.product(FG(j+1,:),H(2:end,:),terms);
                    indx = indx + p-i-j;
                end
            end
        end
    end
end