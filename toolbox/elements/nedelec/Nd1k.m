classdef Nd1k < FE
    %Nd1k - Nedelec first family element of degree k
    %   OBJ = Nd1k(DEGREE,MESH) creates the H(curl)-conforming Nedelec
    %   first family finite element of polynomial degree DEGREE >= 1 on
    %   the triangular mesh MESH. The hierarchical basis of Beuchler,
    %   Pillwein and Zaglmayr consists of DEGREE edge functions per
    %   edge, whose tangential components are continuous across the
    %   edges, and of interior bubble functions with vanishing
    %   tangential trace, DEGREE*(DEGREE+2) basis functions per element
    %   in total. The degrees of freedom are the coefficients in this
    %   basis; those on the Dirichlet edges are fixed.
    %
    %   Nd1k functions:
    %       getBasis                - Evaluate the basis functions
    %       getRotBasis             - Evaluate the rot of the basis
    %       getCurlRotBasis         - Evaluate the curl of the rot of the basis
    %       projectBoundary         - Project tangential data onto boundary
    %       intBoundaryMatrix       - Integrand of the boundary mass matrix
    %       intBoundaryRHS          - Integrand of the boundary load vector
    %       getEdgeRefBasis         - Evaluate the edge reference basis
    %       getEdgeRotRefBasis      - Evaluate the edge rot reference basis
    %       getEdgeCurlRotRefBasis  - Evaluate the edge curl-rot ref basis
    %       getElemRefBasis         - Evaluate the element reference basis
    %       getElemRotRefBasis      - Evaluate the element rot ref basis
    %       getElemCurlRotRefBasis  - Evaluate the element curl-rot ref basis
    %       p1RefBasis              - Evaluate the P1 reference basis
    %       p1GradRefBasis          - Evaluate gradients of the P1 basis
    %
    %   Nd1k properties:
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
    %   See also Nd2k, Nd1k3d, RTk, FE

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
        function obj = Nd1k(degree,mesh)
            %Nd1k - Create a Nedelec first family finite element object
            if degree < 1
                error('fem4m:Nd1k:invalidDegree','Nd1 degree must be >= 1');
            end
            obj.degree = degree;
            obj.mesh = mesh;
            obj.nrBasisFunc = 2*(max(1,3*(degree-1))+max(0,degree-2)*max(0,degree-3)/2)+factorial(degree)/(factorial(degree-1));
            %obj = prepareCoefficients(obj);
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
            %   covariant Piola-mapped Nedelec first family basis
            %   functions in all elements at the reference point PTS_REF
            %   (1-by-2). VAL has the size nrElems x nrBasisFunc x 2. PTS
            %   and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,2);
            val(:,1:3*obj.degree,:) = ones(obj.mesh.nrElems,1).*obj.getEdgeRefBasis(obj.degree,x);
            for k = 1 : 2 : obj.degree
                val(:,3*(k-1)+[1,2,3],1) = val(:,3*(k-1)+[1,2,3],1).*obj.negativeEdges;
                val(:,3*(k-1)+[1,2,3],2) = val(:,3*(k-1)+[1,2,3],2).*obj.negativeEdges;
            end
            val(:,3*obj.degree+1 : obj.nrBasisFunc,:) = ones(obj.mesh.nrElems,1).*obj.getElemRefBasis(obj.degree,x,obj.nrBasisFunc);

            % covariant Piola transformation
            val = matMul(val,obj.jacobian);
            val = (0.5./obj.mesh.area4e).*val;
        end

        function val = getRotBasis(obj,~,x,~)
            %getRotBasis - Evaluate the rot of the basis
            %   VAL = getRotBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   rots of the Nedelec first family basis functions in all
            %   elements at the reference point PTS_REF (1-by-2). VAL has
            %   the size nrElems x nrBasisFunc. PTS and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc);
            val(:,1:3*obj.degree) = ones(obj.mesh.nrElems,1).*obj.getEdgeRotRefBasis(obj.degree,x);
            for k = 1 : 2 : obj.degree
                val(:,3*(k-1)+[1,2,3]) = val(:,3*(k-1)+[1,2,3]).*obj.negativeEdges;
            end
            val(:,3*obj.degree+1 : obj.nrBasisFunc) = ones(obj.mesh.nrElems,1).*obj.getElemRotRefBasis(obj.degree,x,obj.nrBasisFunc);

            % Piola transformation
            val = (0.5./obj.mesh.area4e).*val;
        end

        function val = getCurlRotBasis(obj,~,x,~)
            %getCurlRotBasis - Evaluate the curl of the rot of the basis
            %   VAL = getCurlRotBasis(OBJ,PTS,PTS_REF,PARTS) evaluates
            %   curl(rot(v)) = (d_y rot(v), -d_x rot(v)) of the Nedelec
            %   first family basis functions v in all elements at the
            %   reference point PTS_REF (1-by-2). Only the interior bubble
            %   functions contribute. VAL has the size
            %   nrElems x nrBasisFunc x 2. PTS and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,2);
            %             val(:,1:3*obj.degree,:) = ones(obj.mesh.nrElems,1).*obj.getEdgeCurlRotRefBasis(obj.degree,x);
            %             for k = 1 : 2 : obj.degree
            %                 val(:,3*(k-1)+[1,2,3],1) = val(:,3*(k-1)+[1,2,3],1).*obj.negativeEdges;
            %                 val(:,3*(k-1)+[1,2,3],2) = val(:,3*(k-1)+[1,2,3],2).*obj.negativeEdges;
            %             end
            val(:,3*obj.degree+1 : obj.nrBasisFunc,:) = ones(obj.mesh.nrElems,1).*obj.getElemCurlRotRefBasis(obj.degree,x,obj.nrBasisFunc);

            % transform back to gradient
            val = matMul(val,permute([0 1 ; -1 0 ],[3 1 2]));
            % transform gradient
            c = obj.mesh.transformation;
            val = matMul(val,permute(c,[1 3 2]));
            % transform back to curl
            val = matMul(val,permute([0 -1 ; 1 0 ],[3 1 2]));
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
            obj.nrBasis4ed = obj.degree;
            obj.nrBasis4e = obj.nrBasisFunc-3*obj.nrBasis4ed;
        end

        %% Organisation of the degrees of freedom
        function obj = dofManager(obj)
            %dofManager - Assemble the combined degrees of freedom
            %   OBJ = dofManager(OBJ) enumerates the edge and element
            %   degrees of freedom of the Nedelec element with
            %   FE.dofManager, fixes them on the Dirichlet boundary, and
            %   records the polynomial degree of each local degree of
            %   freedom.
            obj = dofManager@FE(obj);
            obj.fixedNodes = obj.dof4Db;
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
            %   Jacobi-based Nedelec first family edge basis functions of
            %   degree P on the reference triangle at the point X
            %   (1-by-2). VAL has the size 1 x 3*P x 2.
            val = Nd1k.edgeFunctions(p,x,'');
        end

        function val = getEdgeRotRefBasis(p,x)
            %getEdgeRotRefBasis - Evaluate the edge rot reference basis
            %   VAL = getEdgeRotRefBasis(P,X) evaluates the rots of the
            %   Nedelec first family edge basis functions of degree P at the
            %   point X (1-by-2), which are -2 for the lowest order and zero
            %   for all higher orders. VAL has the size 1 x 3*P.
            val = Nd1k.edgeFunctions(p,x,'curl');
        end

        function val = getEdgeCurlRotRefBasis(p,~)
            %getEdgeCurlRotRefBasis - Evaluate the edge curl-rot ref basis
            %   VAL = getEdgeCurlRotRefBasis(P,X) returns curl(rot(v)) of
            %   the Nedelec first family edge basis functions v of degree
            %   P, which vanish identically, as a 1 x 3*P x 2 array of
            %   zeros. X is ignored.
            val = zeros(1,3*p,2);
        end

        %% Elem Basis helper functions
        function val = getElemRefBasis(p,x,nrBasisFunc)
            %getElemRefBasis - Evaluate the element reference basis
            %   VAL = getElemRefBasis(P,X,NRBASISFUNC) evaluates the
            %   interior bubble basis functions of the Nedelec first
            %   family element of degree P with NRBASISFUNC basis
            %   functions on the reference triangle at the point X
            %   (1-by-2). VAL has the size 1 x NRBASISFUNC-3*P x 2.
            val = Nd1k.elemFunctions(p,x,nrBasisFunc,'');
        end

        function val = getElemRotRefBasis(p,x,nrBasisFunc)
            %getElemRotRefBasis - Evaluate the element rot ref basis
            %   VAL = getElemRotRefBasis(P,X,NRBASISFUNC) evaluates the
            %   rots of the interior bubble basis functions of the
            %   Nedelec first family element of degree P with NRBASISFUNC
            %   basis functions at the point X (1-by-2). VAL has the size
            %   1 x NRBASISFUNC-3*P.
            val = Nd1k.elemFunctions(p,x,nrBasisFunc,'curl');
        end

        function val = getElemCurlRotRefBasis(p,x,nrBasisFunc)
            %getElemCurlRotRefBasis - Evaluate the element curl-rot ref basis
            %   VAL = getElemCurlRotRefBasis(P,X,NRBASISFUNC) evaluates
            %   curl(rot(v)) of the interior bubble basis functions v of
            %   the Nedelec first family element of degree P with
            %   NRBASISFUNC basis functions at the point X (1-by-2). VAL
            %   has the size 1 x NRBASISFUNC-3*P x 2.
            val = Nd1k.elemFunctions(p,x,nrBasisFunc,'curlrot');
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

    methods (Static, Access = {?Nd1k,?Nd2k,?RTk,?BDMk})

        % The basis functions are built from the scaled integrated Jacobi
        % factors u_i = L_i^alpha(lambda_b,lambda_a+lambda_b) and
        % v_j = L_j^alpha(lambda_3,1) of BarycentricJacobi and the
        % Whitney functions w_ab = lambda_a grad(lambda_b) -
        % lambda_b grad(lambda_a). Rotated by perp = (-d_y,d_x), they
        % also give the Raviart-Thomas and Brezzi-Douglas-Marini bases.

        function val = rotFreeElemFunctions(p,x)
            %rotFreeElemFunctions - Evaluate the rot-free cell basis
            %   VAL = rotFreeElemFunctions(P,X) evaluates the rot-free
            %   cell basis functions grad(u_i v_j) of degree P with
            %   alpha = 0 and 2i-1, i = 2,...,P-1, j = 1,...,P-i, at the
            %   point X, one row per function.
            lambda = BarycentricJacobi.coordinates(x);
            terms = BarycentricJacobi.leibnizTerms(BarycentricJacobi.multiIndices(2,1),1);
            U = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(2,[1 2],1,2),p-1,0);
            chainV = BarycentricJacobi.chainRule(3,1:3,1,2);
            val = zeros(0,2);
            for i = 2 : p-1
                UV = BarycentricJacobi.product(U(i+1,:),BarycentricJacobi.factor(lambda,chainV,p-i,2*i-1),terms);
                val = [val; UV(2:end,2:3)]; %#ok<AGROW>
            end
        end

        function val = typeIElemFunctions(p,x,op,alpha)
            %typeIElemFunctions - Evaluate the type I cell basis
            %   VAL = typeIElemFunctions(P,X,OP,ALPHA) evaluates the
            %   cell basis functions 2 v_j w_21, j = 1,...,P-1, with the
            %   factor parameter ALPHA at the point X, or their rots for
            %   OP = 'curl' or curl(rot) for OP = 'curlrot', one row per
            %   function.
            [lambda,gradLambda] = BarycentricJacobi.coordinates(x);
            m = 1 + strcmp(op,'curlrot');
            V = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(3,1:3,m,2),p-1,alpha);
            val = 2*Nd1k.timesWhitney(V(2:end,:),lambda,gradLambda,2,1,op);
        end

        function val = typeIIElemFunctions(p,x,op)
            %typeIIElemFunctions - Evaluate the type II cell basis
            %   VAL = typeIIElemFunctions(P,X,OP) evaluates the cell
            %   basis functions u_i v_j w_13 with alpha = 1 and 2i-1,
            %   i = 1,...,P-1, j = 0,...,P-i-1, at the point X, or their
            %   rots for OP = 'curl' or curl(rot) for OP = 'curlrot', one
            %   row per function.
            [lambda,gradLambda] = BarycentricJacobi.coordinates(x);
            m = 1 + strcmp(op,'curlrot');
            terms = BarycentricJacobi.leibnizTerms(BarycentricJacobi.multiIndices(2,m),m);
            U = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(2,[1 2],m,2),p-1,1);
            chainV = BarycentricJacobi.chainRule(3,1:3,m,2);
            val = zeros(0,2-strcmp(op,'curl'));
            for i = 1 : p-1
                UV = BarycentricJacobi.product(U(i+1,:),BarycentricJacobi.factor(lambda,chainV,p-i-1,2*i-1),terms);
                val = [val; Nd1k.timesWhitney(UV,lambda,gradLambda,1,3,op)]; %#ok<AGROW>
            end
        end
    end

    methods (Static, Access = private)

        function val = edgeFunctions(p,x,op)
            %edgeFunctions - Evaluate the edge basis or its rot
            %   VAL = edgeFunctions(P,X,OP) evaluates the edge basis
            %   functions w_ba and grad(u_i), i = 2,...,P, with alpha = 0
            %   of the edges [a b] at the point X, or their rots for
            %   OP = 'curl'. VAL has the size 1 x 3*P x 2, or 1 x 3*P for
            %   the rots.
            edgeNodes = [1 2; 2 3; 3 1];
            [lambda,gradLambda] = BarycentricJacobi.coordinates(x);
            if isempty(op)
                val = zeros(1,3*p,2);
            else
                val = zeros(1,3*p);
            end
            for e = 1 : 3
                a = edgeNodes(e,1);
                b = edgeNodes(e,2);
                [w,rotw] = BarycentricJacobi.whitney(lambda,gradLambda,b,a);
                if isempty(op)
                    val(1,e,:) = w;
                    U = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(b,[a b],1,2),p,0);
                    val(1,e+3:3:end,:) = U(3:end,2:3);
                else
                    val(1,e) = rotw;
                end
            end
        end

        function val = elemFunctions(p,x,nrBasisFunc,op)
            %elemFunctions - Evaluate the cell basis or its derivatives
            %   VAL = elemFunctions(P,X,NRBASISFUNC,OP) evaluates the cell
            %   basis functions of rotFreeElemFunctions, of
            %   typeIElemFunctions with ALPHA = 1, and of
            %   typeIIElemFunctions at the point X, or their rots for
            %   OP = 'curl' or curl(rot) for OP = 'curlrot'. VAL has the
            %   size 1 x NRBASISFUNC-3*P x 2, or 1 x NRBASISFUNC-3*P for
            %   the rots (OP = 'curl').
            free = Nd1k.rotFreeElemFunctions(p,x);
            if ~isempty(op)
                free = zeros(size(free,1),2-strcmp(op,'curl'));
            end
            val = [free; Nd1k.typeIElemFunctions(p,x,op,1); Nd1k.typeIIElemFunctions(p,x,op)];
            val = reshape(val,1,nrBasisFunc-3*p,size(val,2));
        end

        function val = timesWhitney(F,lambda,gradLambda,a,b,op)
            %timesWhitney - Products with a Whitney function
            %   VAL = timesWhitney(F,LAMBDA,GRADLAMBDA,A,B,OP) returns the
            %   products f_k w of the functions f_k with the Whitney
            %   function w = w_AB for OP = '', their rots for OP = 'curl',
            %   or curl(rot(f_k w)) = (d_y r_k,-d_x r_k) with
            %   r_k = rot(f_k w) for OP = 'curlrot', where the rows of F
            %   hold the derivatives of f_k as in
            %   BarycentricJacobi.multiIndices(2,2), up to order one
            %   except for 'curlrot'.
            [w,rotw] = BarycentricJacobi.whitney(lambda,gradLambda,a,b);
            if ~strcmp(op,'curlrot')
                val = BarycentricJacobi.timesField(F(:,1:3),w,rotw,op);
                return;
            end
            % r = f_x w_2 - f_y w_1 + f rot(w), Dw(m,k) = d_k w_m
            Dw = gradLambda(b,:)'*gradLambda(a,:) - gradLambda(a,:)'*gradLambda(b,:);
            rx = F(:,4)*w(2) + F(:,2)*(Dw(2,1)+rotw) - F(:,6)*w(1) - F(:,3)*Dw(1,1);
            ry = F(:,6)*w(2) + F(:,2)*Dw(2,2) - F(:,5)*w(1) - F(:,3)*(Dw(1,2)-rotw);
            val = [ry, -rx];
        end
    end
end