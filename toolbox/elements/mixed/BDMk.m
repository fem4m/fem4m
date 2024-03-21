classdef BDMk < FE
    %BDMk - Brezzi-Douglas-Marini element of degree k
    %   OBJ = BDMk(DEGREE,MESH) creates the Brezzi-Douglas-Marini
    %   H(div)-conforming vector-valued finite element of full
    %   polynomial degree DEGREE >= 1 on the triangular mesh MESH with
    %   the hierarchical edge and cell basis functions of Beuchler,
    %   Pillwein, and Zaglmayr. Each element has (DEGREE+1)*(DEGREE+2)
    %   basis functions, DEGREE+1 per edge. A DEGREE < 1 raises the
    %   error fem4m:BDMk:invalidDegree.
    %
    %   BDMk functions:
    %       getBasis           - Evaluate the Brezzi-Douglas-Marini basis
    %       getDivBasis        - Evaluate the divergence of the basis functions
    %       projectNbBoundary  - Project data onto the Neumann boundary
    %       intBoundaryMatrix  - Integrand of the boundary mass matrix
    %       intBoundaryRHS     - Integrand of the boundary load vector
    %       getEdgeRefBasis    - Evaluate the edge reference basis
    %       getEdgeDivRefBasis - Evaluate the edge div reference basis
    %       getElemRefBasis    - Evaluate the element reference basis
    %       getElemDivRefBasis - Evaluate the element div ref basis
    %       p1RefBasis         - Evaluate the P1 reference basis
    %       p1GradRefBasis     - Evaluate the gradients of the P1 basis
    %       p1PerpGradRefBasis - Evaluate perpendicular gradients of P1
    %
    %   BDMk properties:
    %       dof4e          - Local-to-global dofs per element
    %       dof4ed         - Local-to-global degrees of freedom per edge
    %       nrDoF          - Total number of degrees of freedom
    %       mesh           - Mesh on which the element is defined
    %       freeNodes      - Indices of the free degrees of freedom
    %       fixedNodes     - Indices of the fixed degrees of freedom
    %       degree         - Polynomial degree of the element
    %       nrBasisFunc    - Number of basis functions per element
    %       negativeEdges  - Orientation sign of each edge per element
    %       degree4DoF     - Polynomial degree of each local dof
    %       dof4Nb         - Dofs on the Neumann boundary
    %       nrBasis4n      - Number of nodal basis functions
    %       nrBasis4ed     - Number of edge basis functions
    %       nrBasis4e      - Number of interior basis functions
    %       divDegree4DoF  - Degree of the divergence of each local dof
    %       isRefElemBased - True if the basis lives on the reference element
    %
    %   See also RTk, HBDMk, BDMk3d, FEMatrix

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
        negativeEdges % orientation sign of each edge per element
        degree4DoF     % polynomial degree of each local dof
        dof4Nb         % dofs on the Neumann boundary
        nrBasis4n      % number of nodal basis functions
        nrBasis4ed     % number of edge basis functions
        nrBasis4e      % number of interior basis functions
        divDegree4DoF  % degree of the divergence of each local dof
        isRefElemBased = true % true if the basis lives on the reference element
    end

    methods

        %% Constructor
        function obj = BDMk(degree,mesh)
            %BDMk - Create a Brezzi-Douglas-Marini element object
            if degree < 1
                error('fem4m:BDMk:invalidDegree','BDM degree must be >= 1');
            end
            obj.degree = degree;
            obj.mesh = mesh;
            obj.nrBasisFunc = 2*(3+3*(degree-1)+(degree-1)*max(0,degree-2)/2);
            obj = obj.dofCounts();
            obj = dofManager(obj);
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
            %getBasis - Evaluate the Brezzi-Douglas-Marini basis
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   Piola-mapped Brezzi-Douglas-Marini basis functions of all
            %   elements at the point PTS_REF (1-by-2) of the reference
            %   triangle. VAL is nrElems-by-nrBasisFunc-by-2 (by-3 on a
            %   MeshSurface). PTS and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,2);
            val(:,1:3*(obj.degree+1),:) = ones(obj.mesh.nrElems,1).*obj.getEdgeRefBasis(obj.degree,x);
            for k = 1 : 2 : obj.degree + 1
                val(:,3*(k-1)+[1,2,3],1) = val(:,3*(k-1)+[1,2,3],1).*obj.negativeEdges;
                val(:,3*(k-1)+[1,2,3],2) = val(:,3*(k-1)+[1,2,3],2).*obj.negativeEdges;
            end
            val(:,3*(obj.degree+1)+1 : obj.nrBasisFunc,:) = ones(obj.mesh.nrElems,1).*obj.getElemRefBasis(obj.degree,x,obj.nrBasisFunc);

            % Piola transformation
            if isa(obj.mesh,'MeshSurface')
                area = obj.mesh.area4e(x);
            else
                area = obj.mesh.area4e;
            end
            val = matMul(val,obj.mesh.jacobian(x));
            val = (0.5./area).*val;
        end

        function val = getDivBasis(obj,~,x,~)
            %getDivBasis - Evaluate the divergence of the basis functions
            %   VAL = getDivBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   divergences of the Brezzi-Douglas-Marini basis functions
            %   of all elements at the point PTS_REF (1-by-2) of the
            %   reference triangle. VAL is nrElems-by-nrBasisFunc. PTS
            %   and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc);
            val(:,1:3*(obj.degree+1)) = ones(obj.mesh.nrElems,1).*obj.getEdgeDivRefBasis(obj.degree,x);
            for k = 1 : 2 : obj.degree + 1
                val(:,3*(k-1)+[1,2,3]) = val(:,3*(k-1)+[1,2,3]).*obj.negativeEdges;
            end
            val(:,3*(obj.degree+1)+1 : obj.nrBasisFunc) = ones(obj.mesh.nrElems,1).*obj.getElemDivRefBasis(obj.degree,x,obj.nrBasisFunc);

            % Piola transformation
            if isa(obj.mesh,'MeshSurface')
                area = obj.mesh.area4e(x);
            else
                area = obj.mesh.area4e;
            end
            val = (0.5./area).*val;
        end

        function val = projectNbBoundary(obj,g,gDegree)
            %projectNbBoundary - Project data onto the Neumann boundary
            %   VAL = projectNbBoundary(OBJ,G) computes the L2 projection
            %   of the normal flux G onto the normal traces of the basis
            %   functions on the Neumann edges mesh.ed4Nb by solving the
            %   boundary mass matrix system for their degrees of freedom
            %   dof4ed(mesh.ed4Nb,:). G is a function handle
            %   G(PTS,NORMALS), or G(PTS,PTS_REF,INDICES,NORMALS) if it
            %   takes four arguments, that returns the normal flux at the
            %   points PTS (N-by-2) of the Neumann edges INDICES as an
            %   N-by-1 vector; NORMALS (N-by-2) are the outer unit
            %   normals and PTS_REF is the reference coordinate on the
            %   edges. VAL is an nrDoF-by-1 vector that is zero outside
            %   the Neumann edge dofs.
            %
            %   VAL = projectNbBoundary(OBJ,G,GDEGREE) also specifies the
            %   polynomial degree GDEGREE of G, which sets the quadrature
            %   degree OBJ.degree+GDEGREE of the right-hand side. It
            %   defaults to OBJ.degree, also for an empty GDEGREE.
            if nargin < 3 || isempty(gDegree)
                gDegree = obj.degree;
            end
            normals = obj.mesh.normal4ed;
            normals = normals(obj.mesh.ed4Nb,:);
            M = integrate(@(pts,pts_ref,indices)BDMk.intBoundaryMatrix(pts,pts_ref,indices,obj,obj.mesh,normals),...
                obj.mesh.n4ed(obj.mesh.ed4Nb,:),obj.mesh.ed4Nb,obj.mesh,2*obj.degree);
            e4ed = obj.mesh.e4ed(obj.mesh.ed4Nb,1);
            dofs = obj.dof4e(e4ed,:);
            [I,J] = localDoFtoGlobalDoF(dofs);
            M = sparse(I(:),J(:),M(:),obj.nrDoF,obj.nrDoF);
            b = integrate(@(pts,pts_ref,indices)BDMk.intBoundaryRHS(pts,pts_ref,indices,obj,obj.mesh,normals,g),...
                obj.mesh.n4ed(obj.mesh.ed4Nb,:),obj.mesh.ed4Nb,obj.mesh,obj.degree+gDegree);
            b = accumarray(dofs(:),b(:),[obj.nrDoF 1]);
            freeNodesNb = unique(obj.dof4ed(obj.mesh.ed4Nb,:));
            val = zeros(obj.nrDoF,1);
            perm = symamd(M(freeNodesNb,freeNodesNb));
            freeNodesNb = freeNodesNb(perm);
            val(freeNodesNb) = M(freeNodesNb,freeNodesNb)\b(freeNodesNb);
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
            %   degrees of freedom of the Brezzi-Douglas-Marini element
            %   with FE.dofManager for the counts nrBasis4n, nrBasis4ed,
            %   and nrBasis4e set by the constructor, fixes the dofs
            %   dof4Nb of the Neumann edges (fixedNodes, all others are
            %   freeNodes), and records the polynomial degree of each
            %   local degree of freedom in degree4DoF and the polynomial
            %   degree of its divergence in divDegree4DoF (-1 for
            %   divergence-free basis functions).
            obj = dofManager@FE(obj);
            obj.fixedNodes = obj.dof4Nb;
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
            % determine divDegree4DoF (useful for divergence constraints)
            % degree == -1 -> div-free
            obj.divDegree4DoF = zeros(1,size(obj.dof4e,2));
            obj.divDegree4DoF(1:3) = 0;
            obj.divDegree4DoF(4:3*(obj.degree+1)) = -1;
            indx = 3*(obj.degree+1)+1;
            for i = 2 : obj.degree
                for j = 1 : obj.degree - i + 1
                    obj.divDegree4DoF(indx) = -1; % div-free
                    indx = indx + 1;
                end
            end
            for j = 1 : obj.degree - 1
                obj.divDegree4DoF(indx) = j; % non-div-free
                indx = indx + 1;
            end
            for i = 2 : obj.degree
                for j = 1 : obj.degree + 1 - i
                    obj.divDegree4DoF(indx) = i+j-2; % non-div-free
                    indx = indx + 1;
                end
            end
        end

    end

    methods(Static)

        function val = intBoundaryMatrix(~,pts_ref,indices,rt,mesh,normals)
            %intBoundaryMatrix - Integrand of the boundary mass matrix
            %   VAL = intBoundaryMatrix(PTS,PTS_REF,INDICES,BDM,MESH,NORMALS)
            %   returns the integrand of the local boundary mass matrices
            %   for integrate: for each edge in INDICES the
            %   nrBasisFunc-by-nrBasisFunc matrix of the products of the
            %   normal components of all basis functions of the BDMk
            %   object BDM on the first adjacent element mesh.e4ed(:,1),
            %   evaluated at the reference coordinate PTS_REF of the
            %   edges. NORMALS (N-by-2) holds the unit normals of the N
            %   edges INDICES and VAL is N-by-nrBasisFunc-by-nrBasisFunc.
            %   PTS is ignored.
            basis = getEdgeValues4RefElemFunktion(pts_ref,@(x)rt.getBasis([],x,[]),mesh);
            basis = basis(indices,:,:,1);
            basis = matMul(basis,normals);
            val = matMul(basis,permute(basis,[1 3 2]));
        end

        function val = intBoundaryRHS(pts,pts_ref,indices,rt,mesh,normals,g)
            %intBoundaryRHS - Integrand of the boundary load vector
            %   VAL = intBoundaryRHS(PTS,PTS_REF,INDICES,BDM,MESH,NORMALS,G)
            %   returns the integrand of the boundary load vector for
            %   integrate: for each edge in INDICES the values of G at
            %   the points PTS (N-by-2) times the normal components of
            %   all basis functions of the BDMk object BDM on the first
            %   adjacent element, evaluated at the reference coordinate
            %   PTS_REF of the edges, an N-by-nrBasisFunc array. G is
            %   called as G(PTS,PTS_REF,INDICES,NORMALS) if it takes four
            %   arguments and as G(PTS,NORMALS) otherwise.
            basis = getEdgeValues4RefElemFunktion(pts_ref,@(x)rt.getBasis([],x,[]),mesh);
            basis = basis(indices,:,:,1);
            basis = matMul(basis,normals);
            if nargin(g) == 4
                curG = g(pts,pts_ref,indices,normals);
            else
                curG = g(pts,normals);
            end
            val = curG.*basis;
        end

        % S. Beuchler, V. Pillwein, S. Zaglmayr, Sparsity optimized high
        % order finite element functions for H(div) on simplices, Numer.
        % Math., Volume 122, Issue 2, pp. 197-225

        %% Edge Basis helper functions
        function val = getEdgeRefBasis(p,x)
            %getEdgeRefBasis - Evaluate the edge reference basis
            %   VAL = getEdgeRefBasis(P,X) evaluates the 3*(P+1)
            %   integrated Jacobi-based Brezzi-Douglas-Marini edge basis
            %   functions of degree P at the point X (1-by-2) of the
            %   reference triangle, the RTk edge functions of degree P+1.
            %   VAL is 1-by-3*(P+1)-by-2.
            val = RTk.getEdgeRefBasis(p+1,x);
        end

        function val = getEdgeDivRefBasis(p,x)
            %getEdgeDivRefBasis - Evaluate the edge div reference basis
            %   VAL = getEdgeDivRefBasis(P,X) evaluates the divergences
            %   of the Brezzi-Douglas-Marini edge basis functions of
            %   degree P at the point X (1-by-2) of the reference
            %   triangle, which are 2 for the three lowest-order
            %   functions and zero otherwise. VAL is 1-by-3*(P+1).
            val = RTk.getEdgeDivRefBasis(p+1,x);
        end

        %% Elem Basis helper functions
        function val = getElemRefBasis(p,x,nrBasisFunc)
            %getElemRefBasis - Evaluate the element reference basis
            %   VAL = getElemRefBasis(P,X,NRBASISFUNC) evaluates the
            %   interior bubble basis functions of the
            %   Brezzi-Douglas-Marini element of degree P with
            %   NRBASISFUNC basis functions at the point X (1-by-2) of
            %   the reference triangle. VAL is
            %   1-by-(NRBASISFUNC-3*(P+1))-by-2.
            val = BDMk.elemFunctions(p,x,nrBasisFunc,'');
        end

        function val = getElemDivRefBasis(p,x,nrBasisFunc)
            %getElemDivRefBasis - Evaluate the element div ref basis
            %   VAL = getElemDivRefBasis(P,X,NRBASISFUNC) evaluates the
            %   divergences of the interior bubble basis functions of the
            %   Brezzi-Douglas-Marini element of degree P with
            %   NRBASISFUNC basis functions at the point X (1-by-2) of
            %   the reference triangle; the first P*(P-1)/2 are zero.
            %   VAL is 1-by-(NRBASISFUNC-3*(P+1)).
            val = BDMk.elemFunctions(p,x,nrBasisFunc,'div');
        end

        %% P1 basis functions
        function val = p1RefBasis(x)
            %p1RefBasis - Evaluate the P1 reference basis
            %   VAL = p1RefBasis(X) evaluates the three linear nodal
            %   basis functions of the reference triangle at the points X
            %   (N-by-2, one row per point). VAL is N-by-3.
            val = zeros(size(x,1),3);
            val(:,2) = x(:,1);
            val(:,3) = x(:,2);
            val(:,1) = 1 - x(:,1) - x(:,2);
        end

        function val = p1GradRefBasis(x)
            %p1GradRefBasis - Evaluate the gradients of the P1 basis
            %   VAL = p1GradRefBasis(X) returns the constant gradients of
            %   the three linear nodal basis functions of the reference
            %   triangle at the points X (N-by-2). VAL is N-by-3-by-2.
            val = zeros(size(x,1),3,2);
            val(:,2,1) = 1;
            val(:,3,2) = 1;
            val(:,1,:) = -1;
        end

        function val = p1PerpGradRefBasis(x)
            %p1PerpGradRefBasis - Evaluate perpendicular gradients of P1
            %   VAL = p1PerpGradRefBasis(X) returns the constant
            %   perpendicular gradients (-d_y phi, d_x phi) = -curl(phi)
            %   of the three linear nodal basis functions phi of the
            %   reference triangle at the points X (N-by-2). VAL is
            %   N-by-3-by-2.
            val = zeros(size(x,1),3,2);
            val(:,2,2) = 1;
            val(:,3,1) = -1;
            val(:,1,1) = 1;
            val(:,1,2) = -1;
        end


    end

    methods (Static, Access = private)

        function val = elemFunctions(p,x,nrBasisFunc,op)
            %elemFunctions - Evaluate the cell basis or its divergence
            %   VAL = elemFunctions(P,X,NRBASISFUNC,OP) evaluates the cell
            %   basis functions at the point X (1-by-2), or their
            %   divergences for OP = 'div': the div-free cell functions
            %   of RTk of degree P+1, the type I cell functions of RTk of
            %   degree P, and v_j perp(grad(u_i)) with perp(v) = (-v_2,v_1), i.e.
            %   v_j (-d_y u_i, d_x u_i), and the factors
            %   u_i = L_i(lambda_2,lambda_1+lambda_2) and
            %   v_j = L_j^(2i-1)(lambda_3,1), i = 2,...,P, j = 1,...,P+1-i
            %   (type II). VAL is 1-by-(NRBASISFUNC-3*(P+1))-by-2, or
            %   1-by-(NRBASISFUNC-3*(P+1)) for the divergences.
            lambda = BarycentricJacobi.coordinates(x);
            div = strcmp(op,'div');
            U = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(2,[1 2],1,2),p,0);
            chainV = BarycentricJacobi.chainRule(3,1:3,1,2);
            typeII = zeros(0,2-div);
            for i = 2 : p
                V = BarycentricJacobi.factor(lambda,chainV,p+1-i,2*i-1);
                perpGradU = [-U(i+1,3), U(i+1,2)];
                if div
                    typeII = [typeII; V(2:end,2:3)*perpGradU']; %#ok<AGROW>
                else
                    typeII = [typeII; V(2:end,1)*perpGradU]; %#ok<AGROW>
                end
            end
            if div
                free = zeros(size(Nd1k.rotFreeElemFunctions(p+1,x),1),1);
                val = [free; -Nd1k.typeIElemFunctions(p,x,'curl',3); typeII];
            else
                val = [Nd1k.rotFreeElemFunctions(p+1,x); Nd1k.typeIElemFunctions(p,x,'',3)];
                val = [-val(:,2), val(:,1); typeII];
            end
            val = reshape(val,1,nrBasisFunc-3*(p+1),size(val,2));
        end
    end
end