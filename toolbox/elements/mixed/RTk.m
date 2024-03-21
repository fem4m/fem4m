classdef RTk < FE
    %RTk - Raviart-Thomas element of degree k
    %   OBJ = RTk(DEGREE,MESH) creates the Raviart-Thomas H(div)-
    %   conforming vector-valued finite element of polynomial degree
    %   DEGREE >= 1 on the triangular mesh MESH with the hierarchical
    %   edge and cell basis functions of Beuchler, Pillwein, and
    %   Zaglmayr. DEGREE = 1 is the lowest-order element; each element
    %   has DEGREE*(DEGREE+2) basis functions, DEGREE per edge. A
    %   DEGREE < 1 raises the error fem4m:RTk:invalidDegree.
    %
    %   RTk functions:
    %       getBasis           - Evaluate the Raviart-Thomas basis functions
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
    %   RTk properties:
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
    %       isRefElemBased - True if the basis lives on the reference element
    %
    %   See also BDMk, RTk3d, HRTk, FEMatrix

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
        isRefElemBased = true % true if the basis lives on the reference element
    end

    methods

        %% Constructor
        function obj = RTk(degree,mesh)
            %RTk - Create a Raviart-Thomas finite element object
            if degree < 1
                error('fem4m:RTk:invalidDegree','RT degree must be >= 1');
            end
            obj.degree = degree;
            obj.mesh = mesh;
            obj.nrBasisFunc = 2*(max(1,3*(degree-1))+max(0,degree-2)*max(0,degree-3)/2)+factorial(degree)/(factorial(degree-1));
            %obj = prepareCoefficients(obj);
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
            %getBasis - Evaluate the Raviart-Thomas basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   Piola-mapped Raviart-Thomas basis functions of all
            %   elements at the point PTS_REF (1-by-2) of the reference
            %   triangle. VAL is nrElems-by-nrBasisFunc-by-2 (by-3 on a
            %   MeshSurface). PTS and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,2);
            val(:,1:3*obj.degree,:) = ones(obj.mesh.nrElems,1).*obj.getEdgeRefBasis(obj.degree,x);
            for k = 1 : 2 : obj.degree
                val(:,3*(k-1)+[1,2,3],1) = val(:,3*(k-1)+[1,2,3],1).*obj.negativeEdges;
                val(:,3*(k-1)+[1,2,3],2) = val(:,3*(k-1)+[1,2,3],2).*obj.negativeEdges;
            end
            val(:,3*obj.degree+1 : obj.nrBasisFunc,:) = ones(obj.mesh.nrElems,1).*obj.getElemRefBasis(obj.degree,x,obj.nrBasisFunc);

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
            %   divergences of the Raviart-Thomas basis functions of all
            %   elements at the point PTS_REF (1-by-2) of the reference
            %   triangle. VAL is nrElems-by-nrBasisFunc. PTS and PARTS
            %   are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc);
            val(:,1:3*obj.degree) = ones(obj.mesh.nrElems,1).*obj.getEdgeDivRefBasis(obj.degree,x);
            for k = 1 : 2 : obj.degree
                val(:,3*(k-1)+[1,2,3]) = val(:,3*(k-1)+[1,2,3]).*obj.negativeEdges;
            end
            val(:,3*obj.degree+1 : obj.nrBasisFunc) = ones(obj.mesh.nrElems,1).*obj.getElemDivRefBasis(obj.degree,x,obj.nrBasisFunc);

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
            M = integrate(@(pts,pts_ref,indices)RTk.intBoundaryMatrix(pts,pts_ref,indices,obj,obj.mesh,normals),...
                obj.mesh.n4ed(obj.mesh.ed4Nb,:),obj.mesh.ed4Nb,obj.mesh,2*obj.degree);
            e4ed = obj.mesh.e4ed(obj.mesh.ed4Nb,1);
            dofs = obj.dof4e(e4ed,:);
            [I,J] = localDoFtoGlobalDoF(dofs);
            M = sparse(I(:),J(:),M(:),obj.nrDoF,obj.nrDoF);
            b = integrate(@(pts,pts_ref,indices)RTk.intBoundaryRHS(pts,pts_ref,indices,obj,obj.mesh,normals,g),...
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
            obj.nrBasis4ed = obj.degree;
            obj.nrBasis4e = obj.nrBasisFunc-3*obj.nrBasis4ed;
        end

        %% Organisation of the degrees of freedom
        function obj = dofManager(obj)
            %dofManager - Assemble the combined degrees of freedom
            %   OBJ = dofManager(OBJ) enumerates the edge and element
            %   degrees of freedom of the Raviart-Thomas element with
            %   FE.dofManager for the counts nrBasis4n, nrBasis4ed, and
            %   nrBasis4e set by the constructor, fixes the dofs dof4Nb
            %   of the Neumann edges (fixedNodes, all others are
            %   freeNodes), and records the polynomial degree of each
            %   local degree of freedom in degree4DoF.
            obj = dofManager@FE(obj);
            obj.fixedNodes = obj.dof4Nb;
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
            % determine degree4DoF (useful for hierarchical FEM)
            obj.degree4DoF = zeros(1,size(obj.dof4e,2));
            for p = 1 : obj.degree
                obj.degree4DoF(3*(p-1)+(1:3)) = p;
            end
            indx = 3*obj.degree+1;
            for i = 2 : obj.degree - 1
                for j = 1 : obj.degree - i
                    obj.degree4DoF(indx) = i+j;
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

        function val = intBoundaryMatrix(~,pts_ref,indices,rt,mesh,normals)
            %intBoundaryMatrix - Integrand of the boundary mass matrix
            %   VAL = intBoundaryMatrix(PTS,PTS_REF,INDICES,RT,MESH,NORMALS)
            %   returns the integrand of the local boundary mass matrices
            %   for integrate: for each edge in INDICES the
            %   nrBasisFunc-by-nrBasisFunc matrix of the products of the
            %   normal components of all basis functions of the RTk
            %   object RT on the first adjacent element mesh.e4ed(:,1),
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
            %   VAL = intBoundaryRHS(PTS,PTS_REF,INDICES,RT,MESH,NORMALS,G)
            %   returns the integrand of the boundary load vector for
            %   integrate: for each edge in INDICES the values of G at
            %   the points PTS (N-by-2) times the normal components of
            %   all basis functions of the RTk object RT on the first
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

        % The Raviart-Thomas functions are the Nd1k functions rotated by
        % perp(v) = (-v_2,v_1), except for the parameter alpha = 3 of the
        % type I cell functions, and div(perp(v)) = -rot(v).

        %% Edge Basis helper functions
        function val = getEdgeRefBasis(p,x)
            %getEdgeRefBasis - Evaluate the edge reference basis
            %   VAL = getEdgeRefBasis(P,X) evaluates the 3*P integrated
            %   Jacobi-based Raviart-Thomas edge basis functions of
            %   degree P at the point X (1-by-2) of the reference
            %   triangle, the Nd1k edge functions v rotated to
            %   (-v_2,v_1). VAL is 1-by-3*P-by-2.
            val = Nd1k.getEdgeRefBasis(p,x);
            val = cat(3,-val(:,:,2),val(:,:,1));
        end

        function val = getEdgeDivRefBasis(p,x)
            %getEdgeDivRefBasis - Evaluate the edge div reference basis
            %   VAL = getEdgeDivRefBasis(P,X) evaluates the divergences
            %   of the Raviart-Thomas edge basis functions of degree P
            %   at the point X (1-by-2) of the reference triangle, which
            %   are 2 for the three lowest-order functions and zero for
            %   all higher orders. VAL is 1-by-3*P.
            val = -Nd1k.getEdgeRotRefBasis(p,x);
        end

        %% Elem Basis helper functions
        function val = getElemRefBasis(p,x,nrBasisFunc)
            %getElemRefBasis - Evaluate the element reference basis
            %   VAL = getElemRefBasis(P,X,NRBASISFUNC) evaluates the
            %   interior bubble basis functions of the Raviart-Thomas
            %   element of degree P with NRBASISFUNC basis functions at
            %   the point X (1-by-2) of the reference triangle. VAL is
            %   1-by-(NRBASISFUNC-3*P)-by-2.
            val = [Nd1k.rotFreeElemFunctions(p,x); Nd1k.typeIElemFunctions(p,x,'',3); ...
                Nd1k.typeIIElemFunctions(p,x,'')];
            val = reshape([-val(:,2), val(:,1)],1,nrBasisFunc-3*p,2);
        end

        function val = getElemDivRefBasis(p,x,nrBasisFunc)
            %getElemDivRefBasis - Evaluate the element div ref basis
            %   VAL = getElemDivRefBasis(P,X,NRBASISFUNC) evaluates the
            %   divergences of the interior bubble basis functions of the
            %   Raviart-Thomas element of degree P with NRBASISFUNC basis
            %   functions at the point X (1-by-2) of the reference
            %   triangle; the first (P-1)*(P-2)/2 are zero. VAL is
            %   1-by-(NRBASISFUNC-3*P).
            nrDivFree = (p-1)*(p-2)/2;
            val = -[zeros(nrDivFree,1); Nd1k.typeIElemFunctions(p,x,'curl',3); ...
                Nd1k.typeIIElemFunctions(p,x,'curl')];
            val = reshape(val,1,nrBasisFunc-3*p);
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

end