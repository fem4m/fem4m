classdef HPk < FE
    %HPk - hp-version Pk element with elementwise degrees
    %   OBJ = HPk(DEGREE,MESH) creates the H1-conforming hp-version Pk
    %   element on the triangular mesh MESH. DEGREE >= 1 is a scalar or
    %   an nrElems-by-1 vector with one polynomial degree per element;
    %   each edge gets the minimum degree of its two elements. The
    %   hierarchical basis consists of the P1 nodal functions, the
    %   integrated Legendre edge functions, and the element bubbles of
    %   integrated Jacobi polynomials of Beuchler and Schoeberl. Every
    %   element carries the basis functions of the maximum degree; those
    %   above the degree of the element or of the edge are unused
    %   (nonDoFNodes) and excluded from freeNodes. Coefficient vectors
    %   and system matrices of the space have allnrDoF entries or rows.
    %
    %   HPk functions:
    %       getBasis                 - Evaluate the basis functions
    %       getGradBasis             - Evaluate the gradients of the basis
    %       getD2Basis               - Evaluate the Hessians of the basis
    %       getD3Basis               - Evaluate third derivatives of the basis
    %       getD4Basis               - Evaluate fourth derivatives of the basis
    %       getBasisCoefficients     - Build the orientation factors
    %       projectDbBoundary        - Project data onto the Dirichlet boundary
    %       intBoundaryMatrix        - Integrand of the boundary mass matrix
    %       intBoundaryRHS           - Integrand of the boundary load vector
    %       getEdgeRefBasis          - Evaluate the edge reference basis
    %       getEdgeGradRefBasis      - Evaluate the edge gradient basis
    %       getEdgeD2RefBasis        - Evaluate the edge Hessian basis
    %       getEdgeD3RefBasis        - Evaluate the edge third-order basis
    %       getEdgeD4RefBasis        - Evaluate the edge fourth-order basis
    %       getElementRefBasis       - Evaluate the element reference basis
    %       getElementGradRefBasis   - Evaluate the element gradient basis
    %       getElementD2RefBasis     - Evaluate the element Hessian basis
    %       getElementD3RefBasis     - Evaluate the element third-order basis
    %       getElementD4RefBasis     - Evaluate the element fourth-order basis
    %       p1RefBasis               - Evaluate the P1 reference basis
    %       p1GradRefBasis           - Evaluate gradients of the P1 basis
    %
    %   HPk properties:
    %       mesh              - Mesh on which the element is defined
    %       freeNodes         - Indices of the free degrees of freedom
    %       fixedNodes        - Indices of the fixed degrees of freedom
    %       degree4e          - Polynomial degree of each element
    %       degree4ed         - Minimum degree of the two elements at each edge
    %       degree            - Maximum of the elementwise degrees
    %       basisCoefficients - Orientation factors per element
    %       nrBasisFunc       - Number of basis functions per element
    %       dof4e             - Local-to-global dofs per element
    %       dof4ed            - Local-to-global dofs per edge
    %       nrDoF             - Number of used degrees of freedom
    %       allnrDoF          - Total number of degrees of freedom
    %       nonDoFNodes       - Unused dofs above the element or edge degree
    %       degree4DoF        - Polynomial degree of each local dof
    %       dof4Db            - Dofs on the Dirichlet boundary
    %       nrBasis4n         - Number of nodal basis functions
    %       nrBasis4ed        - Number of edge basis functions
    %       nrBasis4e         - Number of interior basis functions
    %       isRefElemBased    - Basis lives on the reference element
    %
    %   See also Pk, HPk3d, DGHPk, DGDk, MarkHP, FE

    % Copyright (c) 2026 by the fem4m authors

    properties
        mesh              % mesh on which the element is defined
        freeNodes         % indices of the free degrees of freedom
        fixedNodes        % indices of the fixed degrees of freedom
        degree4e          % polynomial degree of each element
        degree4ed         % minimum degree of the two elements at each edge
        degree            % maximum of the elementwise degrees
        basisCoefficients % orientation factors per element
        nrBasisFunc       % number of basis functions per element
        dof4e             % local-to-global dofs per element
        dof4ed            % local-to-global dofs per edge
        nrDoF             % number of used degrees of freedom
        allnrDoF          % total number of degrees of freedom
        nonDoFNodes       % unused dofs above the element or edge degree
        degree4DoF        % polynomial degree of each local dof
        dof4Db            % dofs on the Dirichlet boundary
        nrBasis4n         % number of nodal basis functions
        nrBasis4ed        % number of edge basis functions
        nrBasis4e         % number of interior basis functions
        isRefElemBased = true % basis lives on the reference element
    end

    methods
        function obj = HPk(degree,mesh)
            %HPk - Create an hp-version Pk finite element object
            if any(degree(:) < 1)
                error('fem4m:HPk:invalidDegree','HPk degree must be at least 1');
            end
            obj.degree4e = degree;
            obj.degree = max(degree);
            if isscalar(obj.degree4e)
                obj.degree4e = obj.degree4e*ones(mesh.nrElems,1);
            end
            obj.degree4ed = obj.degree4e(mesh.e4ed(:,1))*ones(1,2);
            obj.degree4ed(mesh.e4ed(:,2)~=0,2) = obj.degree4e(mesh.e4ed(mesh.e4ed(:,2)~=0,2));
            obj.degree4ed = min(obj.degree4ed,[],2);
            obj.nrBasisFunc = 3 + 3*(obj.degree-1) + (obj.degree-2)*(obj.degree-1)/2;
            obj.mesh = mesh;
            obj = obj.dofCounts();
            obj = obj.dofManager;
            obj = obj.getBasisCoefficients;
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the basis
            %   functions, oriented by basisCoefficients, in all elements
            %   at the reference point PTS_REF (1-by-2). VAL is
            %   nrElems-by-nrBasisFunc. PTS and PARTS are ignored.
            %
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) with six entries
            %   [x1 y1 x2 y2 x3 y3] in PTS_REF (or a 3-by-2 PTS_REF), one
            %   point per local edge of the reference triangle, evaluates
            %   the basis functions of the first element e4ed(:,1) of
            %   every edge at the point of its local edge edNr4ed(:,1).
            %   VAL is nrEdges-by-nrBasisFunc.
            if numel(x) == 6
                if size(x,1) == 1
                    x = reshape(x,2,3)';
                end
                edNr4ed = obj.mesh.edNr4ed;
                if obj.degree == 0
                    val = ones(obj.mesh.nrEdges,obj.nrBasisFunc);
                else
                    val = zeros(obj.mesh.nrEdges,obj.nrBasisFunc);
                end
                if obj.degree > 0
                    p1RefBasis = obj.p1RefBasis(x);
                end
                if obj.degree > 1
                    basisE = obj.getEdgeRefBasis(obj.degree,x);
                end
                for ed = 1:3
                    ind = find(edNr4ed(:,1) == ed);
                    if ~isempty(ind)
                        if obj.degree > 0
                            val(ind,1:3) = ones(numel(ind),1)*p1RefBasis(ed,:);
                        end
                        if obj.degree > 1
                            val(ind,4:3+3*(obj.degree-1)) = ones(numel(ind),1)*basisE(ed,:);
                        end
                        % the element bubbles vanish on the edges
                    end
                end
                val = obj.basisCoefficients(obj.mesh.e4ed(:,1),:).*val;
            else
                val = obj.refBasis(x,0);
            end
        end

        function val = getGradBasis(obj,~,x,~)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   gradients of the basis functions at the reference point
            %   PTS_REF (1-by-2) and maps them to the physical elements.
            %   VAL is nrElems-by-nrBasisFunc-by-d with the dimension d
            %   of the mesh coordinates. PTS and PARTS are ignored.
            %
            %   VAL = getGradBasis(OBJ,PTS,PTS_REF,PARTS) with six
            %   entries in PTS_REF evaluates, as getBasis does, the
            %   gradients in the first element of every edge at the
            %   point of its local edge (nrEdges rows).
            if numel(x) == 6
                if size(x,1) == 1
                    x = reshape(x,2,3)';
                end
                edNr4ed = obj.mesh.edNr4ed;
                val = zeros(obj.mesh.nrEdges,obj.nrBasisFunc,2);
                if obj.degree > 0
                    p1RefBasis = obj.p1GradRefBasis(x);
                end
                if obj.degree > 1
                    basisE = obj.getEdgeGradRefBasis(obj.degree,x);
                end
                if obj.degree > 2
                    basisT = obj.getElementGradRefBasis(obj.degree,x);
                end
                for ed = 1 : 3
                    ind = find(edNr4ed(:,1) == ed);
                    if ~isempty(ind)
                        if obj.degree > 0
                            val(ind,1:3,1:2) = bsxfun(@times,ones(numel(ind),1),p1RefBasis(ed,:,:));
                        end
                        if obj.degree > 1
                            val(ind,4:3+3*(obj.degree-1),1:2) = bsxfun(@times,ones(numel(ind),1),basisE(ed,:,:));
                        end
                        if obj.degree > 2
                            val(ind,3+3*(obj.degree-1)+1:obj.nrBasisFunc,1:2) = bsxfun(@times,ones(numel(ind),1),basisT(ed,:,:));
                        end
                    end
                end
                val = bsxfun(@times,obj.basisCoefficients(obj.mesh.e4ed(:,1),:),val);
                c = obj.mesh.transformation;
                c = c(obj.mesh.e4ed(:,1),:,:);
            else
                val = obj.refBasis(x,1);
                c = obj.mesh.transformation(x);
            end
            val = matMul(val,permute(c,[1 3 2]));
        end

        function val = getD2Basis(obj,~,x,~)
            %getD2Basis - Evaluate the Hessians of the basis
            %   VAL = getD2Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   Hessians [u_xx,u_yy,u_xy] of the basis functions at the
            %   reference point PTS_REF (1-by-2), mapped to the physical
            %   elements (on a MeshSurface [xx yy zz xy xz yz]). On a
            %   curved MeshSurface (degree > 1) FE.curvedSurfaceD2 adds
            %   the derivative of the transformation applied to the
            %   gradients and projects the result onto the tangent
            %   planes. PTS and PARTS are ignored.
            c = obj.mesh.transformation(x);
            val = FE.mapDerivatives(c,obj.refBasis(x,2),2);
            if isa(obj.mesh,'MeshSurface') && obj.mesh.degree > 1
                val = FE.curvedSurfaceD2(obj.mesh,x,c,val,obj.refBasis(x,1));
            end
        end

        function val = getD3Basis(obj,~,x,~)
            %getD3Basis - Evaluate third derivatives of the basis
            %   VAL = getD3Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   third derivatives [u_xxx,u_yyy,u_xxy,u_xyy] of the basis
            %   functions at the reference point PTS_REF (1-by-2),
            %   mapped to the physical elements. Surface meshes get no
            %   curvature correction. PTS and PARTS are ignored.
            val = FE.mapDerivatives(obj.mesh.transformation(x),obj.refBasis(x,3),3);
        end

        function val = getD4Basis(obj,~,x,~)
            %getD4Basis - Evaluate fourth derivatives of the basis
            %   VAL = getD4Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   fourth derivatives [u_xxxx,u_yyyy,u_xxyy,u_xxxy,u_xyyy]
            %   of the basis functions at the reference point PTS_REF
            %   (1-by-2), mapped to the physical elements. Surface meshes
            %   get no curvature correction. PTS and PARTS are ignored.
            val = FE.mapDerivatives(obj.mesh.transformation(x),obj.refBasis(x,4),4);
        end

        function obj = getBasisCoefficients(obj)
            %getBasisCoefficients - Build the orientation factors
            %   OBJ = getBasisCoefficients(OBJ) sets basisCoefficients
            %   (nrElems-by-nrBasisFunc) to 1, except for the factor -1
            %   of the edge functions of odd degree i = 3,5,... in the
            %   first element e4ed(:,1) of every edge, boundary edges
            %   included. The two elements of an interior edge traverse
            %   it in opposite directions, so the odd functions then
            %   agree on the edge. Nothing is flipped for degree <= 2
            %   (the maximum degree).
            obj.basisCoefficients = ones(obj.mesh.nrElems,obj.nrBasisFunc);
            % odd edge basis functions: -1 in the first element of each edge
            if obj.degree > 2
                obj.basisCoefficients(obj.mesh.e4ed(obj.mesh.edNr4ed(:,1)==1,1),...
                    3+2:2:3+1*(obj.degree-1)) = -1;
                obj.basisCoefficients(obj.mesh.e4ed(obj.mesh.edNr4ed(:,1)==2,1),...
                    3+2+(obj.degree-1):2:3+2*(obj.degree-1)) = -1;
                obj.basisCoefficients(obj.mesh.e4ed(obj.mesh.edNr4ed(:,1)==3,1),...
                    3+2+2*(obj.degree-1):2:3+3*(obj.degree-1)) = -1;
            end
        end

        function val = projectDbBoundary(obj,u4Db,uDegree)
            %projectDbBoundary - Project data onto the Dirichlet boundary
            %   VAL = projectDbBoundary(OBJ,U4DB) returns the
            %   allnrDoF-by-1 vector of the L2 projection of the
            %   Dirichlet data onto the traces of the basis functions on
            %   the Dirichlet edges ed4Db. It solves the boundary mass
            %   matrix system for the nodal and edge dofs of these edges
            %   up to the edge degree degree4ed; VAL is zero at the unused
            %   dofs nonDoFNodes and elsewhere. U4DB is a function handle;
            %   U4DB(X) returns the N-by-1 values at the points X
            %   (N-by-2).
            %
            %   VAL = projectDbBoundary(OBJ,U4DB,UDEGREE) also sets the
            %   polynomial degree UDEGREE of the data (default: the
            %   maximum degree), which only enters the quadrature
            %   degree (maximum degree + UDEGREE) of the right-hand side.
            if nargin < 3 || isempty(uDegree)
                uDegree = obj.degree;
            end
            M = integrate(@(pts,pts_ref,indices)HPk.intBoundaryMatrix(pts,pts_ref,indices,obj,obj.mesh),...
                obj.mesh.n4ed(obj.mesh.ed4Db,:),obj.mesh.ed4Db,obj.mesh,2*obj.degree);
            e4ed = obj.mesh.e4ed(obj.mesh.ed4Db,1);
            dofs = obj.dof4e(e4ed,:);
            [I,J] = localDoFtoGlobalDoF(dofs);
            M = sparse(I(:),J(:),M(:),obj.allnrDoF,obj.allnrDoF);
            b = integrate(@(pts,pts_ref,indices)HPk.intBoundaryRHS(pts,pts_ref,indices,obj,obj.mesh,u4Db),...
                obj.mesh.n4ed(obj.mesh.ed4Db,:),obj.mesh.ed4Db,obj.mesh,obj.degree+uDegree);
            b = accumarray(dofs(:),b(:),[obj.allnrDoF 1]);
            dofs = [obj.dof4ed,obj.mesh.n4ed];
            val = zeros(obj.allnrDoF,1);
            freeNodesDb = setdiff(dofs(obj.mesh.ed4Db,:),obj.nonDoFNodes);
            perm = symamd(M(freeNodesDb,freeNodesDb));
            freeNodesDb = freeNodesDb(perm);
            val(freeNodesDb) = M(freeNodesDb,freeNodesDb)\b(freeNodesDb);
        end

    end

    methods (Access = private)

        function val = refBasis(obj,x,order)
            %refBasis - Evaluate the oriented reference basis derivatives
            %   VAL = refBasis(OBJ,X,ORDER) evaluates the derivatives of
            %   order ORDER = 0,...,4 of the basis functions on the
            %   reference triangle at the reference point X (1-by-2),
            %   oriented per element by basisCoefficients. VAL has one
            %   row per element, one column per function, and the
            %   components of Monomials.multiIndices(2,ORDER), one for
            %   ORDER = 0.
            op = {'','Grad','D2','D3','D4'};
            op = op{order+1};
            nrComponents = order+1;
            val = zeros(obj.nrBasisFunc,nrComponents);
            if order == 0
                val(1:3,:) = reshape(obj.p1RefBasis(x),[],1);
            elseif order == 1
                val(1:3,:) = reshape(obj.p1GradRefBasis(x),[],2);
            end
            if obj.degree > 1
                val(4:3+3*(obj.degree-1),:) = reshape( ...
                    obj.(['getEdge',op,'RefBasis'])(obj.degree,x),[],nrComponents);
            end
            if obj.degree > 2
                val(3+3*(obj.degree-1)+1:end,:) = reshape( ...
                    obj.(['getElement',op,'RefBasis'])(obj.degree,x),[],nrComponents);
            end
            val = obj.basisCoefficients.*reshape(val,1,obj.nrBasisFunc,nrComponents);
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of nodal, edge, and interior
            %   degrees of freedom. The constructor calls it right
            %   before dofManager, so that subclasses may override it.
            obj.nrBasis4n = 1;
            obj.nrBasis4ed = obj.degree-1;
            obj.nrBasis4e = (obj.degree-1)*(obj.degree-2)/2;
        end

        function obj = dofManager(obj)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ) enumerates the nodal, edge, and
            %   interior degrees of freedom with FE.dofManager (the edge
            %   dofs of each element listed edge by edge), sets
            %   degree4DoF, fixes the dofs dof4Db of the Dirichlet
            %   boundary, and collects in nonDoFNodes the unused dofs
            %   above degree4e (bubbles) and degree4ed (edges). allnrDoF
            %   counts all dofs, nrDoF the used ones, and freeNodes are
            %   the used dofs off the Dirichlet boundary.
            obj = dofManager@FE(obj,false,false,[false true false]);
            obj.degree4DoF(1:3) = 1;
            % 3*(p-1) for edges E1;E2;E3
            if obj.degree > 1
                obj.degree4DoF(4:3+3*(obj.degree-1)) = [2:obj.degree,2:obj.degree,2:obj.degree];
            end
            % (p-2)*(p-1)/2 interior, sorted by degree (elementDerivatives)
            if obj.degree > 2
                for p = 3 : obj.degree
                    obj.degree4DoF(:,3+3*(obj.degree-1)+(1+(p-3)*(p-2)/2:(p-2)*(p-1)/2)) = p;
                end
            end
            obj.allnrDoF = obj.nrDoF;
            obj.fixedNodes = obj.dof4Db;
            degr4ed = reshape(obj.degree4ed(obj.mesh.ed4e),obj.mesh.nrElems,3);
            nonDoFNodes = [];
            for curDegree = 1 : obj.degree
                %inner dofs
                nonDoFNodes = [nonDoFNodes;reshape(obj.dof4e(obj.degree4e==curDegree,...
                    3+3*(obj.degree-1)+(curDegree-2)*(curDegree-1)/2+1:end),[],1)]; %#ok<AGROW>
                % E1
                nonDoFNodes = [nonDoFNodes;reshape(obj.dof4e(degr4ed(:,1)==curDegree,...
                    3+(curDegree-1)+1:3+(obj.degree-1)),[],1)]; %#ok<AGROW>
                % E2
                nonDoFNodes = [nonDoFNodes;reshape(obj.dof4e(degr4ed(:,2)==curDegree,...
                    3+(obj.degree-1)+(curDegree-1)+1:3+2*(obj.degree-1)),[],1)]; %#ok<AGROW>
                % E3
                nonDoFNodes = [nonDoFNodes;reshape(obj.dof4e(degr4ed(:,3)==curDegree,...
                    3+2*(obj.degree-1)+(curDegree-1)+1:3+3*(obj.degree-1)),[],1)]; %#ok<AGROW>
            end
            nonDoFNodes = unique(nonDoFNodes);
            obj.nonDoFNodes = nonDoFNodes;
            obj.nrDoF = obj.allnrDoF - numel(nonDoFNodes);
            obj.freeNodes = setdiff(1:obj.allnrDoF,[obj.fixedNodes;nonDoFNodes]);
        end

    end

    methods(Static)

        function val = intBoundaryMatrix(~,pts_ref,indices,pk,mesh)
            %intBoundaryMatrix - Integrand of the boundary mass matrix
            %   VAL = intBoundaryMatrix(PTS,PTS_REF,INDICES,PK,MESH)
            %   evaluates the basis functions of PK in the first element
            %   e4ed(:,1) of the edges INDICES at the edge-local
            %   quadrature point PTS_REF and returns their products, the
            %   integrand of the boundary mass matrix (one
            %   nrBasisFunc-by-nrBasisFunc matrix per edge). PTS is
            %   ignored.
            basis = getEdgeValues4RefElemFunktion(pts_ref,@(x)pk.getBasis([],x,[]),mesh);
            basis = basis(indices,:,1);
            val = matMul(basis,permute(basis,[1 3 2]));
        end

        function val = intBoundaryRHS(pts,pts_ref,indices,pk,mesh,u4Db)
            %intBoundaryRHS - Integrand of the boundary load vector
            %   VAL = intBoundaryRHS(PTS,PTS_REF,INDICES,PK,MESH,U4DB)
            %   evaluates the basis functions of PK in the first element
            %   e4ed(:,1) of the edges INDICES at the edge-local
            %   quadrature point PTS_REF and multiplies them by the
            %   Dirichlet data U4DB(PTS) at the points PTS, the integrand
            %   of the right-hand side (one row per edge). A U4DB with
            %   three input arguments is called as U4DB(PTS,X,INDICES)
            %   with an undefined X and errors.
            basis = getEdgeValues4RefElemFunktion(pts_ref,@(x)pk.getBasis([],x,[]),mesh);
            basis = basis(indices,:,1);
            if nargin(u4Db) == 3
                curU4Db = u4Db(pts,x,indices);
            else
                curU4Db = u4Db(pts);
            end
            val = curU4Db.*basis;
        end

        % S. Beuchler, J. Schoeberl,
        % New shape functions for triangular p-FEM using integrated Jacobi
        % polynomials. Numer. Math. 103 (2006), no. 3, 339-366.

        %% Edge Basis helper functions
        function val = getEdgeRefBasis(p,x)
            %getEdgeRefBasis - Evaluate the edge reference basis
            %   VAL = getEdgeRefBasis(P,X) evaluates the edge functions
            %   of degree i = 2,...,P of the three edges of the
            %   reference triangle, edge by edge, at the reference
            %   points X (N-by-2). VAL is N-by-3*(P-1).
            val = HPk.edgeDerivatives(p,x,0);
        end

        function val = getEdgeGradRefBasis(p,x)
            %getEdgeGradRefBasis - Evaluate the edge gradient basis
            %   VAL = getEdgeGradRefBasis(P,X) evaluates the
            %   gradients [dx,dy] of the edge functions of
            %   getEdgeRefBasis at the reference points X (N-by-2). VAL
            %   is N-by-3*(P-1)-by-2.
            %[dx,dy]
            val = HPk.edgeDerivatives(p,x,1);
        end

        function val = getEdgeD2RefBasis(p,x)
            %getEdgeD2RefBasis - Evaluate the edge Hessian basis
            %   VAL = getEdgeD2RefBasis(P,X) evaluates the Hessians
            %   [dxx,dyy,dxy] of the edge functions of getEdgeRefBasis
            %   at the reference points X (N-by-2). VAL is
            %   N-by-3*(P-1)-by-3.
            % [dxx,dyy,dxy]
            val = HPk.edgeDerivatives(p,x,2);
            val = val(:,:,[1 3 2]);
        end

        function val = getEdgeD3RefBasis(p,x)
            %getEdgeD3RefBasis - Evaluate the edge third-order basis
            %   VAL = getEdgeD3RefBasis(P,X) evaluates the third
            %   derivatives [dxxx,dyyy,dxxy,dxyy] of the edge functions
            %   of getEdgeRefBasis at the reference points X (N-by-2).
            %   VAL is N-by-3*(P-1)-by-4.
            % [dxxx,dyyy,dxxy,dxyy]
            val = HPk.edgeDerivatives(p,x,3);
            val = val(:,:,[1 4 2 3]);
        end

        function val = getEdgeD4RefBasis(p,x)
            %getEdgeD4RefBasis - Evaluate the edge fourth-order basis
            %   VAL = getEdgeD4RefBasis(P,X) evaluates the fourth
            %   derivatives [d_xxxx,d_yyyy,d_xxyy,d_xxxy,d_xyyy] of
            %   the edge functions of getEdgeRefBasis at the reference
            %   points X (N-by-2). VAL is N-by-3*(P-1)-by-5.
            %  [d_xxxx,d_yyyy,d_xxyy,d_xxxy,d_xyyy]
            val = HPk.edgeDerivatives(p,x,4);
            val = val(:,:,[1 5 3 2 4]);
        end

        %% Elem basis helper functions
        function val = getElementRefBasis(p,x)
            %getElementRefBasis - Evaluate the element reference basis
            %   VAL = getElementRefBasis(P,X) evaluates the element
            %   bubbles of degree 3,...,P, products of integrated
            %   Legendre and integrated Jacobi polynomials sorted by
            %   degree, at the reference points X (N-by-2). VAL is
            %   N-by-(P-1)*(P-2)/2.
            val = HPk.elementDerivatives(p,x,0);
        end

        function val = getElementGradRefBasis(p,x)
            %getElementGradRefBasis - Evaluate the element gradient basis
            %   VAL = getElementGradRefBasis(P,X) evaluates the
            %   gradients [dx,dy] of the element bubbles of
            %   getElementRefBasis at the reference points X (N-by-2).
            %   VAL is N-by-(P-1)*(P-2)/2-by-2.
            val = HPk.elementDerivatives(p,x,1);
        end

        function val = getElementD2RefBasis(p,x)
            %getElementD2RefBasis - Evaluate the element Hessian basis
            %   VAL = getElementD2RefBasis(P,X) evaluates the
            %   Hessians [dxx,dyy,dxy] of the element bubbles of
            %   getElementRefBasis at the reference points X (N-by-2).
            %   VAL is (P-1)*(P-2)/2-by-3-by-N, one row per function.
            val = HPk.elementDerivatives(p,x,2);
            val = permute(val(:,:,[1 3 2]),[2 3 1]);
        end

        function val = getElementD3RefBasis(p,x)
            %getElementD3RefBasis - Evaluate the element third-order basis
            %   VAL = getElementD3RefBasis(P,X) evaluates the third
            %   derivatives [u_xxx,u_yyy,u_xxy,u_xyy] of the element
            %   bubbles of getElementRefBasis at the reference points X
            %   (N-by-2). VAL is (P-1)*(P-2)/2-by-4-by-N.
            %  order: [u_xxx,u_yyy,u_xxy, u_xyy]
            val = HPk.elementDerivatives(p,x,3);
            val = permute(val(:,:,[1 4 2 3]),[2 3 1]);
        end

        function val = getElementD4RefBasis(p,x)
            %getElementD4RefBasis - Evaluate the element fourth-order basis
            %   VAL = getElementD4RefBasis(P,X) evaluates the fourth
            %   derivatives [u_xxxx,u_yyyy,u_xxyy,u_xxxy,u_xyyy] of
            %   the element bubbles of getElementRefBasis at the
            %   reference points X (N-by-2). VAL is
            %   (P-1)*(P-2)/2-by-5-by-N.
            %  order: [u_xxxx,u_yyyy,u_xxyy, u_xxxy, u_xyyy]
            val = HPk.elementDerivatives(p,x,4);
            val = permute(val(:,:,[1 5 3 2 4]),[2 3 1]);
        end

        %% P1 basis functions
        function val = p1RefBasis(x)
            %p1RefBasis - Evaluate the P1 reference basis
            %   VAL = p1RefBasis(X) evaluates the linear nodal basis
            %   functions [1-x-y, x, y] at the reference points X
            %   (N-by-2). VAL is N-by-3.
            val = zeros(size(x,1),3);
            val(:,2) = x(:,1);
            val(:,3) = x(:,2);
            val(:,1) = 1 - x(:,1) - x(:,2);
        end

        function val = p1GradRefBasis(x)
            %p1GradRefBasis - Evaluate gradients of the P1 basis
            %   VAL = p1GradRefBasis(X) returns the constant gradients
            %   of the linear nodal basis functions as an N-by-3-by-2
            %   array for the N rows of X.
            val = zeros(size(x,1),3,2);
            val(:,2,1) = 1;
            val(:,3,2) = 1;
            val(:,1,:) = -1;
        end

    end

    methods (Static, Access = private)

        function val = edgeDerivatives(p,x,m)
            %edgeDerivatives - Evaluate derivatives of the edge basis
            %   VAL = edgeDerivatives(P,X,M) evaluates the derivatives
            %   d^M/(dx^(M-b) dy^b), b = 0,...,M, of the edge
            %   functions L_i(x,1-y), x g_i(y), and (-1)^i (1-x-y) g_i(y),
            %   i = 2,...,P, at the reference points X, where L_i(x,s) =
            %   s^i L_i(2x/s-1) is the scaled integrated Legendre
            %   polynomial and g_i(y) = L_i(2y-1)/(1-y) =
            %   -2y/(i-1) P_{i-2}^(1,1)(2y-1). VAL has the size
            %   size(X,1) x 3*(P-1) x (M+1).
            y = x(:,2);
            x = x(:,1);
            i = 2 : p;
            % G(:,i-1,b+1) = d^b/dy^b g_i(y)
            G = zeros(numel(x),p-1,m+1);
            for b = 0 : m
                G(:,:,b+1) = y.*2^b.*JacobiPolynomials.poly(2*y-1,p-2,1,1,b)';
                if b > 0
                    G(:,:,b+1) = G(:,:,b+1) + b*2^(b-1)*JacobiPolynomials.poly(2*y-1,p-2,1,1,b-1)';
                end
                G(:,:,b+1) = -2./(i-1).*G(:,:,b+1);
            end
            val = zeros(numel(x),3*(p-1),m+1);
            for b = 0 : m
                a = m-b;
                % e1: L_i(x,1-y)
                L = JacobiPolynomials.scaledIntegratedPoly(x,1-y,p,0,a,b)';
                val(:,1:p-1,b+1) = (-1)^b*L(:,i+1);
                % e2: x g_i(y), e3: (-1)^i (1-x-y) g_i(y)
                if a == 0
                    val(:,p:2*(p-1),b+1) = x.*G(:,:,b+1);
                    e3 = (1-x-y).*G(:,:,b+1);
                    if b > 0
                        e3 = e3 - b*G(:,:,b);
                    end
                elseif a == 1
                    val(:,p:2*(p-1),b+1) = G(:,:,b+1);
                    e3 = -G(:,:,b+1);
                else
                    e3 = zeros(numel(x),p-1);
                end
                val(:,2*(p-1)+1:3*(p-1),b+1) = (-1).^i.*e3;
            end
        end

        function val = elementDerivatives(p,x,m)
            %elementDerivatives - Evaluate derivatives of the element basis
            %   VAL = elementDerivatives(P,X,M) evaluates the derivatives
            %   d^M/(dx^(M-b) dy^b), b = 0,...,M, of the element bubble
            %   functions phi_ij = L_i(x,1-y) L_j^(2i-1)(2y-1), i =
            %   2,...,P-1, j = 1,...,P-i, at the reference points X,
            %   where L_i(x,s) = s^i L_i(2x/s-1) is the scaled integrated
            %   Legendre polynomial. VAL has the size size(X,1) x
            %   (P-1)*(P-2)/2 x (M+1). The functions are sorted by their
            %   degree i+j and then by i, so that the first (q-1)*(q-2)/2
            %   of them are the bubbles of degree at most q, as
            %   degree4DoF and nonDoFNodes assume.
            y = x(:,2);
            x = x(:,1);
            % L{a+1,c+1} = d^a/dx^a d^c/ds^c L_i(x,s) with s = 1-y
            L = cell(m+1);
            for a = 0 : m
                for c = 0 : m-a
                    L{a+1,c+1} = JacobiPolynomials.scaledIntegratedPoly(x,1-y,p-1,0,a,c)';
                end
            end
            val = zeros(numel(x),(p-1)*(p-2)/2,m+1);
            for i = 2 : p-1
                % K{d+1} = d^d/dy^d L_j^(2i-1)(2y-1) for j = 0,...,p-i
                K = cell(m+1,1);
                for d = 0 : m
                    K{d+1} = 2^d*JacobiPolynomials.integratedPoly(2*y-1,p-i,2*i-1,d)';
                end
                % columns of phi_ij behind the (k-3)*(k-2)/2 bubbles of
                % degree < k = i+j and the i-2 of degree k with smaller i
                k = i + (1:p-i);
                col = (k-3).*(k-2)/2 + i-1;
                % Leibniz rule in y with ds/dy = -1
                for b = 0 : m
                    for c = 0 : b
                        val(:,col,b+1) = val(:,col,b+1) + (-1)^c*nchoosek(b,c)*L{m-b+1,c+1}(:,i+1).*K{b-c+1}(:,2:end);
                    end
                end
            end
        end
    end
end