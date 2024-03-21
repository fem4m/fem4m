classdef Nd2k < FE
    %Nd2k - Nedelec second family element of degree k
    %   OBJ = Nd2k(DEGREE,MESH) creates the H(curl)-conforming Nedelec
    %   second family finite element of polynomial degree DEGREE >= 1 on
    %   the triangular mesh MESH. The hierarchical basis of Beuchler,
    %   Pillwein and Zaglmayr consists of DEGREE+1 edge functions per
    %   edge, whose tangential components are continuous across the
    %   edges, and of interior bubble functions with vanishing
    %   tangential trace, (DEGREE+1)*(DEGREE+2) basis functions per
    %   element in total. The degrees of freedom are the coefficients
    %   in this basis; those on the Dirichlet edges are fixed.
    %
    %   Nd2k functions:
    %       getBasis            - Evaluate the basis functions
    %       getRotBasis         - Evaluate the rot of the basis
    %       projectBoundary     - Project tangential data onto boundary
    %       intBoundaryMatrix   - Integrand of the boundary mass matrix
    %       intBoundaryRHS      - Integrand of the boundary load vector
    %       getEdgeRefBasis     - Evaluate the edge reference basis
    %       getEdgeRotRefBasis  - Evaluate the edge rot reference basis
    %       getElemRefBasis     - Evaluate the element reference basis
    %       getElemRotRefBasis  - Evaluate the element rot ref basis
    %       p1RefBasis          - Evaluate the P1 reference basis
    %       p1GradRefBasis      - Evaluate gradients of the P1 basis
    %
    %   Nd2k properties:
    %       dof4e          - Local-to-global dofs per element
    %       dof4ed         - Local-to-global degrees of freedom per edge
    %       nrDoF          - Total number of degrees of freedom
    %       mesh           - Mesh on which the element is defined
    %       freeNodes      - Indices of the free degrees of freedom
    %       fixedNodes     - Indices of the fixed degrees of freedom
    %       degree         - Polynomial degree of the element
    %       nrBasisFunc    - Number of basis functions per element
    %       jacobian       - Cofactor det(J)*J^-T of the mesh Jacobian
    %       negativeEdges  - Orientation sign of each edge per element
    %       degree4DoF     - Polynomial degree of each local dof
    %       dof4Db         - Dofs on the Dirichlet boundary
    %       nrBasis4n      - Number of nodal basis functions
    %       nrBasis4ed     - Number of basis functions per edge
    %       nrBasis4e      - Number of interior basis functions per element
    %       isRefElemBased - Basis lives on the reference element
    %
    %   See also Nd1k, Nd2k3d, RTk, FE

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
        jacobian       % cofactor det(J)*J^-T of the mesh Jacobian
        negativeEdges  % orientation sign of each edge per element
        degree4DoF     % polynomial degree of each local dof
        dof4Db         % dofs on the Dirichlet boundary
        nrBasis4n      % number of nodal basis functions
        nrBasis4ed     % number of basis functions per edge
        nrBasis4e      % number of interior basis functions per element
        isRefElemBased = true % basis lives on the reference element
    end

    methods

        %% Constructor
        function obj = Nd2k(degree,mesh)
            %Nd2k - Create a Nedelec second family finite element object
            if degree < 1
                error('fem4m:Nd2k:invalidDegree','Nd2 degree must be >= 1');
            end
            obj.degree = degree;
            obj.mesh = mesh;
            obj.nrBasisFunc = 2*(3+3*(degree-1)+(degree-1)*max(0,degree-2)/2);
            obj = obj.dofCounts();
            obj = dofManager(obj);
            %covariant transformation (cofactor of J, i.e. det(J)*J^-T)
            obj.jacobian = obj.mesh.jacobian;
            tmp = obj.jacobian(:,1,1);
            obj.jacobian(:,1,1) = obj.jacobian(:,2,2);
            obj.jacobian(:,2,2) = tmp;
            tmp = obj.jacobian(:,1,2);
            obj.jacobian(:,1,2) = -obj.jacobian(:,2,1);
            obj.jacobian(:,2,1) = -tmp;
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
            %   covariant Piola-mapped Nedelec second family basis
            %   functions in all elements at the reference point PTS_REF
            %   (1-by-2). VAL has the size nrElems x nrBasisFunc x 2. PTS
            %   and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,2);
            val(:,1:3*(obj.degree+1),:) = ones(obj.mesh.nrElems,1).*obj.getEdgeRefBasis(obj.degree,x);
            for k = 1 : 2 : obj.degree+1
                val(:,3*(k-1)+[1,2,3],1) = val(:,3*(k-1)+[1,2,3],1).*obj.negativeEdges;
                val(:,3*(k-1)+[1,2,3],2) = val(:,3*(k-1)+[1,2,3],2).*obj.negativeEdges;
            end
            val(:,3*(obj.degree+1)+1 : obj.nrBasisFunc,:) = ones(obj.mesh.nrElems,1).*obj.getElemRefBasis(obj.degree,x,obj.nrBasisFunc);

            % covariant Piola transformation
            val = matMul(val,obj.jacobian);
            val = (0.5./obj.mesh.area4e).*val;
        end

        function val = getRotBasis(obj,~,x,~)
            %getRotBasis - Evaluate the rot of the basis
            %   VAL = getRotBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   rots of the Nedelec second family basis functions in all
            %   elements at the reference point PTS_REF (1-by-2). VAL has
            %   the size nrElems x nrBasisFunc. PTS and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc);
            val(:,1:3*(obj.degree+1)) = ones(obj.mesh.nrElems,1).*obj.getEdgeRotRefBasis(obj.degree,x);
            for k = 1 : 2 : obj.degree+1
                val(:,3*(k-1)+[1,2,3]) = val(:,3*(k-1)+[1,2,3]).*obj.negativeEdges;
            end
            val(:,3*(obj.degree+1)+1 : obj.nrBasisFunc) = ones(obj.mesh.nrElems,1).*obj.getElemRotRefBasis(obj.degree,x,obj.nrBasisFunc);

            % Piola transformation
            val = (0.5./obj.mesh.area4e).*val;
        end

        function val = projectBoundary(obj,ed4Bd,g,gDegree)
            %projectBoundary - Project tangential data onto boundary edges
            %   VAL = projectBoundary(OBJ,ED4BD,G) computes the L2
            %   projection of the tangential trace G onto the edge degrees
            %   of freedom of the boundary edges ED4BD by solving the
            %   boundary mass matrix system restricted to these degrees of
            %   freedom. The function handle G is called as G(PTS,TANGENTS)
            %   or, if it takes four arguments, as
            %   G(PTS,PTS_REF,INDICES,TANGENTS) with the quadrature points
            %   PTS (N-by-2) on the edges INDICES = ED4BD, the edge-local
            %   reference point PTS_REF and the unit tangents TANGENTS
            %   (N-by-2) of these edges, N = numel(ED4BD), and must return
            %   the tangential component g*t (N-by-1), e.g.
            %   @(pts,tangents)sum(u4Db(pts).*tangents,2). VAL is the
            %   nrDoF-by-1 coefficient vector, zero outside the edge
            %   degrees of freedom of ED4BD.
            %
            %   VAL = projectBoundary(OBJ,ED4BD,G,GDEGREE) also specifies
            %   the polynomial degree GDEGREE of G, which sets the
            %   quadrature degree OBJ.degree+GDEGREE of the right-hand
            %   side (default GDEGREE = OBJ.degree).
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
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of nodal, edge, and interior
            %   degrees of freedom. The constructor calls it right
            %   before dofManager, so that subclasses may override it.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = obj.degree+1;
            obj.nrBasis4e = obj.nrBasisFunc-3*obj.nrBasis4ed;
        end

        %% Organisation of the degrees of freedom
        function obj = dofManager(obj)
            %dofManager - Assemble the combined degrees of freedom
            %   OBJ = dofManager(OBJ) enumerates the edge and element
            %   degrees of freedom of the Nedelec second family
            %   element with FE.dofManager, fixes them on the Dirichlet
            %   boundary, and records the polynomial degree of each
            %   local degree of freedom.
            obj = dofManager@FE(obj);
            obj.fixedNodes = obj.dof4Db;
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
            % determine degree4DoF (useful for hierarchical FEM)
            obj.degree4DoF = ones(1,size(obj.dof4e,2));
            for p = 1 : obj.degree
                obj.degree4DoF(3*p+(1:3)) = p;
            end
            indx = 3*(obj.degree+1)+1;
            for i = 2 : obj.degree
                for j = 1 : obj.degree - i + 1
                    obj.degree4DoF(indx) = i+j-1;
                    indx = indx + 1;
                end
            end
            for j = 1 : obj.degree - 1
                obj.degree4DoF(indx) = j+1;
                indx = indx + 1;
            end
            for i = 2 : obj.degree
                for j = 1 : obj.degree + 1 - i
                    obj.degree4DoF(indx) = i+j-1;
                    indx = indx + 1;
                end
            end
        end

    end

    methods(Static)

        function val = intBoundaryMatrix(~,pts_ref,indices,fe,mesh,~,tangents)
            %intBoundaryMatrix - Integrand of the boundary mass matrix
            %   VAL = intBoundaryMatrix(PTS,PTS_REF,INDICES,FE,MESH,NORMALS,TANGENTS)
            %   returns for the edges INDICES the products (b_i*t)(b_j*t)
            %   of the tangential components of all basis functions b_i of
            %   FE in the first neighbouring element of each edge
            %   (MESH.e4ed(INDICES,1)) at the edge-local reference point
            %   PTS_REF, with the unit tangents TANGENTS (N-by-2,
            %   N = numel(INDICES)). VAL has the size
            %   N x nrBasisFunc x nrBasisFunc. PTS and NORMALS are
            %   ignored.
            basis = getEdgeValues4RefElemFunktion(pts_ref,@(x)fe.getBasis([],x,[]),mesh);
            basis = basis(indices,:,:,1);
            basis = matMul(basis,tangents);
            val = matMul(basis,permute(basis,[1 3 2]));
        end

        function val = intBoundaryRHS(pts,pts_ref,indices,fe,mesh,~,tangents,g)
            %intBoundaryRHS - Integrand of the boundary load vector
            %   VAL = intBoundaryRHS(PTS,PTS_REF,INDICES,FE,MESH,NORMALS,TANGENTS,G)
            %   returns for the edges INDICES the products g*(b_i*t) of
            %   the tangential trace G with the tangential components of
            %   all basis functions b_i of FE in the first neighbouring
            %   element of each edge (MESH.e4ed(INDICES,1)) at the
            %   edge-local reference point PTS_REF and the points PTS
            %   (N-by-2), with the unit tangents TANGENTS (N-by-2,
            %   N = numel(INDICES)). G is called as G(PTS,TANGENTS) or,
            %   if it takes four arguments, as
            %   G(PTS,PTS_REF,INDICES,TANGENTS) and returns an N-by-1
            %   vector. VAL has the size N x nrBasisFunc. NORMALS is
            %   ignored.
            basis = getEdgeValues4RefElemFunktion(pts_ref,@(x)fe.getBasis([],x,[]),mesh);
            basis = basis(indices,:,:,1);
            basis = matMul(basis,tangents);
            if nargin(g) == 4
                curG = g(pts,pts_ref,indices,tangents);
            else
                curG = g(pts,tangents);
            end
            val = curG.*basis;
        end

        % S. Beuchler, V. Pillwein, S. Zaglmayr, Sparsity optimized high
        % order finite element functions for H(curl) on tetrahedra,
        % Advances in Applied Mathematics 50 (2013) 749-769

        %% Edge Basis helper functions
        function val = getEdgeRefBasis(p,x)
            %getEdgeRefBasis - Evaluate the edge reference basis
            %   VAL = getEdgeRefBasis(P,X) evaluates the integrated
            %   Jacobi-based Nedelec second family edge basis functions of
            %   degree P on the reference triangle at the point X
            %   (1-by-2), the Nd1k edge functions of degree P+1. VAL has
            %   the size 1 x 3*(P+1) x 2.
            val = Nd1k.getEdgeRefBasis(p+1,x);
        end

        function val = getEdgeRotRefBasis(p,x)
            %getEdgeRotRefBasis - Evaluate the edge rot reference basis
            %   VAL = getEdgeRotRefBasis(P,X) evaluates the rots of the
            %   Nedelec second family edge basis functions of degree P at
            %   the point X (1-by-2), which are -2 for the lowest order and
            %   zero for all higher orders. VAL has the size 1 x 3*(P+1).
            val = Nd1k.getEdgeRotRefBasis(p+1,x);
        end

        %% Elem Basis helper functions
        function val = getElemRefBasis(p,x,nrBasisFunc)
            %getElemRefBasis - Evaluate the element reference basis
            %   VAL = getElemRefBasis(P,X,NRBASISFUNC) evaluates the
            %   interior bubble basis functions of the Nedelec second
            %   family element of degree P with NRBASISFUNC basis
            %   functions on the reference triangle at the point X
            %   (1-by-2). VAL has the size 1 x NRBASISFUNC-3*(P+1) x 2.
            val = Nd2k.elemFunctions(p,x,nrBasisFunc,'');
        end

        function val = getElemRotRefBasis(p,x,nrBasisFunc)
            %getElemRotRefBasis - Evaluate the element rot ref basis
            %   VAL = getElemRotRefBasis(P,X,NRBASISFUNC) evaluates the
            %   rots of the interior bubble basis functions of the
            %   Nedelec second family element of degree P with NRBASISFUNC
            %   basis functions at the point X (1-by-2). VAL has the size
            %   1 x NRBASISFUNC-3*(P+1).
            val = Nd2k.elemFunctions(p,x,nrBasisFunc,'curl');
        end

        %% P1 basis functions
        function val = p1RefBasis(x)
            %p1RefBasis - Evaluate the P1 reference basis
            %   VAL = p1RefBasis(X) evaluates the three linear nodal basis
            %   functions of the reference triangle at the points X
            %   (N-by-2), one row per point. VAL has the size N x 3.
            val = zeros(size(x,1),3);
            val(:,2) = x(:,1);
            val(:,3) = x(:,2);
            val(:,1) = 1 - x(:,1) - x(:,2);
        end

        function val = p1GradRefBasis(x)
            %p1GradRefBasis - Evaluate gradients of the P1 basis
            %   VAL = p1GradRefBasis(X) returns the constant gradients of
            %   the three linear nodal basis functions of the reference
            %   triangle at the points X (N-by-2). VAL has the size
            %   N x 3 x 2.
            val = zeros(size(x,1),3,2);
            val(:,2,1) = 1;
            val(:,3,2) = 1;
            val(:,1,:) = -1;
        end

    end

    methods (Static, Access = private)

        function val = elemFunctions(p,x,nrBasisFunc,op)
            %elemFunctions - Evaluate the cell basis or its rot
            %   VAL = elemFunctions(P,X,NRBASISFUNC,OP) evaluates the cell
            %   basis functions at the point X, or their rots for
            %   OP = 'curl': the rot-free functions of Nd1k of degree
            %   P+1, the type I functions of Nd1k of degree P, and
            %   v_j grad(u_i) - u_i grad(v_j) with the factors
            %   u_i = L_i(lambda_2,lambda_1+lambda_2) and
            %   v_j = L_j^(2i-1)(lambda_3,1), i = 2,...,P, j = 1,...,P+1-i
            %   (type II). VAL has the size 1 x NRBASISFUNC-3*(P+1) x 2,
            %   or 1 x NRBASISFUNC-3*(P+1) for the rots.
            lambda = BarycentricJacobi.coordinates(x);
            free = Nd1k.rotFreeElemFunctions(p+1,x);
            if ~isempty(op)
                free = zeros(size(free,1),1);
            end
            U = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(2,[1 2],1,2),p,0);
            chainV = BarycentricJacobi.chainRule(3,1:3,1,2);
            typeII = zeros(0,2-~isempty(op));
            for i = 2 : p
                V = BarycentricJacobi.factor(lambda,chainV,p+1-i,2*i-1);
                if isempty(op)
                    typeII = [typeII; V(2:end,1)*U(i+1,2:3) - U(i+1,1)*V(2:end,2:3)]; %#ok<AGROW>
                else
                    typeII = [typeII; 2*(V(2:end,2)*U(i+1,3) - V(2:end,3)*U(i+1,2))]; %#ok<AGROW>
                end
            end
            val = [free; Nd1k.typeIElemFunctions(p,x,op,1); typeII];
            val = reshape(val,1,nrBasisFunc-3*(p+1),size(val,2));
        end
    end
end