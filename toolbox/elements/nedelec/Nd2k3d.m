classdef Nd2k3d < FE
    %Nd2k3d - Nedelec second family element of degree k on tetrahedra
    %   OBJ = Nd2k3d(DEGREE,MESH) creates the tetrahedral Nedelec second
    %   family H(curl)-conforming vector-valued finite element of
    %   polynomial degree DEGREE >= 1 on the tetrahedral mesh MESH. The
    %   space consists of all vector polynomials of degree DEGREE. The
    %   degrees of freedom are the coefficients of the hierarchical
    %   edge, face and cell basis functions of Beuchler, Pillwein and
    %   Zaglmayr; the shared edge and face degrees of freedom make the
    %   space tangentially continuous. The degrees of freedom on the
    %   Dirichlet boundary are fixed. DEGREE < 1 raises the error
    %   fem4m:Nd2k3d:invalidDegree.
    %
    %   Nd2k3d functions:
    %       getBasis            - Evaluate the Nedelec basis functions
    %       getCurlBasis        - Evaluate the curl of the basis functions
    %       getEdgeRefBasis     - Evaluate the edge reference basis
    %       getEdgeCurlRefBasis - Evaluate the edge curl reference basis
    %       getFaceRefBasis     - Evaluate the face reference basis
    %       getFaceCurlRefBasis - Evaluate the face curl reference basis
    %       getElemRefBasis     - Evaluate the element reference basis
    %       getElemCurlRefBasis - Evaluate the element curl ref basis
    %       p1RefBasis          - Evaluate the P1 reference basis
    %       p1GradRefBasis      - Evaluate the gradients of the P1 basis
    %
    %   Nd2k3d properties:
    %       nrDoF          - Total number of degrees of freedom
    %       dof4e          - Local-to-global dofs per element
    %       dof4f          - Local-to-global dof numbers per face
    %       dof4ed         - Local-to-global dof numbers per edge
    %       pkFE           - Placeholder for the associated Pk element
    %       nrBasisFunc    - Number of basis functions per element
    %       degree         - Polynomial degree of the element
    %       fixedNodes     - Indices of the fixed degrees of freedom
    %       dof4Db         - Dofs on the Dirichlet boundary
    %       freeNodes      - Indices of the free degrees of freedom
    %       nrBasis4n      - Number of nodal basis functions
    %       nrBasis4ed     - Number of edge basis functions
    %       nrBasis4f      - Number of face basis functions
    %       nrBasis4e      - Number of element bubble functions
    %       mesh           - Mesh on which the element is defined
    %       isRefElemBased - Basis lives on the reference element
    %       sign4e         - Edge and face orientation signs per element
    %       perm4e         - Cyclic node permutation index of each face
    %
    %   See also Nd2k, Nd1k3d, RTk3d, FE

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = public)
        nrDoF          % total number of degrees of freedom
        dof4e          % local-to-global dofs per element
        dof4f          % local-to-global dof numbers per face
        dof4ed         % local-to-global dof numbers per edge
        pkFE           % placeholder for the associated Pk element
        nrBasisFunc    % number of basis functions per element
        degree         % polynomial degree of the element
        fixedNodes     % indices of the fixed degrees of freedom
        dof4Db         % dofs on the Dirichlet boundary
        freeNodes      % indices of the free degrees of freedom
        nrBasis4n      % number of nodal basis functions
        nrBasis4ed     % number of edge basis functions
        nrBasis4f      % number of face basis functions
        nrBasis4e      % number of element bubble functions
        mesh           % mesh on which the element is defined
        isRefElemBased = true % basis lives on the reference element
    end

    properties(SetAccess = public)
        sign4e            % edge and face orientation signs per element
        perm4e            % cyclic node permutation index of each face
    end

    properties(Access = private)
        jacobian
        transformation
        facePolyIndices
        area4e
    end

    methods
        function obj = Nd2k3d(degree,mesh)
            %Nd2k3d - Create a tetrahedral Nedelec element object
            if degree < 1
                error('fem4m:Nd2k3d:invalidDegree','Nedelec does only work for degree >= 1');
            end
            obj.degree = degree;
            obj = obj.dofCounts();
            obj.nrBasisFunc = 6*obj.nrBasis4ed + 4*obj.nrBasis4f + obj.nrBasis4e;
            obj = obj.getBasisCoefficients(mesh);
            obj = obj.dofManager(mesh);
            obj.transformation = reshape(mesh.transformation(),mesh.nrElems,3,3);
            obj.jacobian = reshape(mesh.jacobian(),mesh.nrElems,3,3);
            obj.area4e = mesh.area4e;
            obj.mesh = mesh;
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the Nedelec basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   covariant Piola-mapped tetrahedral Nedelec second family
            %   basis functions in all elements at the reference point
            %   PTS_REF (1-by-3). VAL is an nrElems-by-nrBasisFunc-by-3
            %   array. PTS and PARTS are ignored.
            nrElems = size(obj.area4e,1);
            val = zeros(nrElems,obj.nrBasisFunc,3);
            val(:,1:6*obj.nrBasis4ed,:) = ones(nrElems,1).*obj.getEdgeRefBasis(obj.degree,x);
            for k = 1 : 2 : obj.nrBasis4ed
                val(:,6*(k-1)+(1:6),1) = val(:,6*(k-1)+(1:6),1).*obj.sign4e(:,1:6);
                val(:,6*(k-1)+(1:6),2) = val(:,6*(k-1)+(1:6),2).*obj.sign4e(:,1:6);
                val(:,6*(k-1)+(1:6),3) = val(:,6*(k-1)+(1:6),3).*obj.sign4e(:,1:6);
            end
            faceValues = ones(nrElems,1).*obj.getFaceRefBasis(obj.degree,x,obj.nrBasis4f);
            for k = 1 : obj.nrBasis4f
                for curFace = 1 : 4
                    for perm = 1 : 3
                        indx = obj.perm4e(:,curFace)==perm;
                        for coordinate = 1 : 3
                            if mod(obj.facePolyIndices(k,1),2)==1
                                val(indx,6*obj.nrBasis4ed + 4*(k-1)+curFace,coordinate) = faceValues(indx,curFace,k,coordinate,perm).*obj.sign4e(indx,6+curFace);
                            else
                                val(indx,6*obj.nrBasis4ed + 4*(k-1)+curFace,coordinate) = faceValues(indx,curFace,k,coordinate,perm);
                            end
                        end
                    end
                end
            end
            val(:,6*obj.nrBasis4ed+4*obj.nrBasis4f+1:end,:) = ones(nrElems,1).*obj.getElemRefBasis(obj.degree,x,obj.nrBasis4e);

            % covariant Piola transformation
            val = matMul(val,obj.transformation);
        end

        function val = getCurlBasis(obj,~,x,~)
            %getCurlBasis - Evaluate the curl of the basis functions
            %   VAL = getCurlBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   curls of the tetrahedral Nedelec second family basis
            %   functions, mapped by the contravariant Piola
            %   transformation, in all elements at the reference point
            %   PTS_REF (1-by-3). VAL is an nrElems-by-nrBasisFunc-by-3
            %   array. PTS and PARTS are ignored.
            nrElems = size(obj.area4e,1);
            val = zeros(nrElems,obj.nrBasisFunc,3);
            val(:,1:6*obj.nrBasis4ed,:) = ones(nrElems,1).*obj.getEdgeCurlRefBasis(obj.degree,x);
            for k = 1 : 2 : obj.nrBasis4ed
                val(:,6*(k-1)+(1:6),1) = val(:,6*(k-1)+(1:6),1).*obj.sign4e(:,1:6);
                val(:,6*(k-1)+(1:6),2) = val(:,6*(k-1)+(1:6),2).*obj.sign4e(:,1:6);
                val(:,6*(k-1)+(1:6),3) = val(:,6*(k-1)+(1:6),3).*obj.sign4e(:,1:6);
            end
            faceValues = ones(nrElems,1).*obj.getFaceCurlRefBasis(obj.degree,x,obj.nrBasis4f);
            for k = 1 : obj.nrBasis4f
                for curFace = 1 : 4
                    for perm = 1 : 3
                        indx = obj.perm4e(:,curFace)==perm;
                        for coordinate = 1 : 3
                            if mod(obj.facePolyIndices(k,1),2)==1
                                val(indx,6*obj.nrBasis4ed + 4*(k-1)+curFace,coordinate) = faceValues(indx,curFace,k,coordinate,perm).*obj.sign4e(indx,6+curFace);
                            else
                                val(indx,6*obj.nrBasis4ed + 4*(k-1)+curFace,coordinate) = faceValues(indx,curFace,k,coordinate,perm);
                            end
                        end
                    end
                end
            end
            val(:,6*obj.nrBasis4ed+4*obj.nrBasis4f+1:end,:) = ones(nrElems,1).*obj.getElemCurlRefBasis(obj.degree,x,obj.nrBasis4e);

            % contravariant Piola transformation of the curl
            val = matMul(val,obj.jacobian);
            val = (1./(6*obj.area4e)).*val;
        end

    end

    methods (Access = private)

        function obj = getBasisCoefficients(obj,mesh)
            %getBasisCoefficients - Set the signs, face permutations and face polynomial indices
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
            % curl free
            for i = 2 : p
                for j = 1 : p + 1 - i
                    obj.facePolyIndices(indx,:) = [i,j];
                    indx = indx + 1;
                end
            end
            % non-curl free
            % Type I
            for j = 1 : p - 1
                obj.facePolyIndices(indx,:) = [1,j];
                indx = indx + 1;
            end
            % Type II
            for i = 2 : p
                for j = 1 : p + 1 - i
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
            %   calls it before dofManager, so that subclasses may
            %   override it.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = obj.degree+1;
            obj.nrBasis4f = (obj.degree-1) + (obj.degree-1)*obj.degree/2 + (obj.degree-1)*obj.degree/2;
            obj.nrBasis4e = (obj.degree-1)*(obj.degree-2)/2 + 3*obj.degree*(obj.degree-1)*max(0,obj.degree-2)/6;
        end

        function obj = dofManager(obj,mesh)
            %dofManager - Assemble the degrees of freedom of the element
            %   OBJ = dofManager(OBJ,MESH) enumerates the edge, face,
            %   and interior degrees of freedom with FE.dofManager and
            %   fixes those on the Dirichlet boundary (dof4Db).
            obj.mesh = mesh;
            obj = dofManager@FE(obj);
            obj.fixedNodes = obj.dof4Db;
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
        end

    end

    methods(Static)

        % S. Beuchler, V. Pillwein, S. Zaglmayr, Sparsity optimized high
        % order finite element functions for H(curl) on tetrahedra,
        % Advances in Applied Mathematics 50 (2013) 749-769

        % ordering of the reference tetrahedron
        % c4n: [ 0 0 0; 1 0 0; 0 1 0; 0 0 1];
        % edges: [1 2], [2 3], [3 1], [1 4], [2 4], [3 4]

        %% Edge Basis helper functions
        function val = getEdgeRefBasis(p,x)
            %getEdgeRefBasis - Evaluate the edge reference basis
            %   VAL = getEdgeRefBasis(P,X) evaluates the integrated
            %   Jacobi-based tetrahedral Nedelec second family edge basis
            %   functions of degree P on the reference tetrahedron at the
            %   point X (1-by-3), the Nd1k3d edge functions of degree P+1.
            %   VAL is a 1-by-6*(P+1)-by-3 array.
            val = Nd1k3d.getEdgeRefBasis(p+1,x);
        end

        function val = getEdgeCurlRefBasis(p,x)
            %getEdgeCurlRefBasis - Evaluate the edge curl reference basis
            %   VAL = getEdgeCurlRefBasis(P,X) evaluates the curls of the
            %   tetrahedral Nedelec second family edge basis functions of
            %   degree P at the point X (1-by-3), which are constant for
            %   the Whitney functions and zero for the gradient functions.
            %   VAL is a 1-by-6*(P+1)-by-3 array.
            val = Nd1k3d.getEdgeCurlRefBasis(p+1,x);
        end

        %% Face Basis helper functions
        % Order: [1 2 3], [1 2 4], [4 2 3], [1 3 4]
        % Permutations: [a b c] = [1 2 3], [2 3 1], [3 1 2] of the face
        % nodes
        function val = getFaceRefBasis(p,x,nrBasisFunc4f)
            %getFaceRefBasis - Evaluate the face reference basis
            %   VAL = getFaceRefBasis(P,X,NRBASISFUNC4F) evaluates the
            %   NRBASISFUNC4F integrated Jacobi-based face basis functions
            %   per face of degree P for all four faces and the three
            %   cyclic node permutations of the tetrahedral Nedelec second
            %   family element at the point X (1-by-3); getBasis obtains
            %   the odd permutations by changing the sign of the
            %   functions with odd index i. VAL is a
            %   1-by-4-by-NRBASISFUNC4F-by-3-by-3 array.
            val = Nd2k3d.faceFunctions(p,x,nrBasisFunc4f,'');
        end

        function val = getFaceCurlRefBasis(p,x,nrBasisFunc4f)
            %getFaceCurlRefBasis - Evaluate the face curl reference basis
            %   VAL = getFaceCurlRefBasis(P,X,NRBASISFUNC4F) evaluates the
            %   curls of the NRBASISFUNC4F face basis functions per face of
            %   degree P for all four faces and the three cyclic node
            %   permutations of the tetrahedral Nedelec second family
            %   element at the point X (1-by-3); getCurlBasis obtains the
            %   odd permutations by changing the sign of the functions
            %   with odd index i. VAL is a
            %   1-by-4-by-NRBASISFUNC4F-by-3-by-3 array.
            val = Nd2k3d.faceFunctions(p,x,nrBasisFunc4f,'curl');
        end

        %% Elem Basis helper functions
        function val = getElemRefBasis(p,x,nrBasisFunc4e)
            %getElemRefBasis - Evaluate the element reference basis
            %   VAL = getElemRefBasis(P,X,NRBASISFUNC4E) evaluates the
            %   NRBASISFUNC4E interior bubble basis functions of degree P
            %   of the tetrahedral Nedelec second family element at the
            %   point X (1-by-3). VAL is a 1-by-NRBASISFUNC4E-by-3 array.
            val = Nd2k3d.elemFunctions(p,x,nrBasisFunc4e,'');
        end

        function val = getElemCurlRefBasis(p,x,nrBasisFunc4e)
            %getElemCurlRefBasis - Evaluate the element curl ref basis
            %   VAL = getElemCurlRefBasis(P,X,NRBASISFUNC4E) evaluates the
            %   curls of the NRBASISFUNC4E interior bubble basis functions
            %   of degree P of the tetrahedral Nedelec second family
            %   element at the point X (1-by-3). VAL is a
            %   1-by-NRBASISFUNC4E-by-3 array.
            val = Nd2k3d.elemFunctions(p,x,nrBasisFunc4e,'curl');
        end

        %% P1 basis functions
        function val = p1RefBasis(x)
            %p1RefBasis - Evaluate the P1 reference basis
            %   VAL = p1RefBasis(X) evaluates the four linear nodal basis
            %   functions of the reference tetrahedron at the points X
            %   (N-by-3, one row per point). VAL is an N-by-4 array.
            val = zeros(size(x,1),4);
            val(:,1) = 1 - x(:,1) - x(:,2) - x(:,3);
            val(:,2) = x(:,1);
            val(:,3) = x(:,2);
            val(:,4) = x(:,3);
        end

        function val = p1GradRefBasis(x)
            %p1GradRefBasis - Evaluate the gradients of the P1 basis
            %   VAL = p1GradRefBasis(X) returns the constant gradients of
            %   the four linear nodal basis functions of the reference
            %   tetrahedron at the points X (N-by-3, one row per point).
            %   VAL is an N-by-4-by-3 array.
            val = zeros(size(x,1),4,3);
            val(:,1,:) = -1;
            val(:,2,1) = 1;
            val(:,3,2) = 1;
            val(:,4,3) = 1;
        end

    end

    methods (Static, Access = private)

        % The basis functions are built from the scaled integrated Jacobi
        % factors u_i = L_i(lambda_b,lambda_a+lambda_b),
        % v_j = L_j^alpha(lambda_c,lambda_a+lambda_b+lambda_c), and
        % z_k = L_k^alpha(lambda_4,1) of BarycentricJacobi and the Whitney
        % functions w_ab = lambda_a grad(lambda_b) - lambda_b grad(lambda_a).

        function val = faceFunctions(p,x,nrBasisFunc4f,op)
            %faceFunctions - Evaluate the face basis or its curl
            %   VAL = faceFunctions(P,X,NRBASISFUNC4F,OP) evaluates the
            %   face basis functions grad(u_i v_j), i = 2,...,P,
            %   j = 1,...,P+1-i (curl-free), -v_j w_ab with alpha = 1,
            %   j = 1,...,P-1 (type I), and v_j grad(u_i) - u_i grad(v_j),
            %   i = 2,...,P, j = 1,...,P+1-i (type II), with alpha = 2i-1
            %   otherwise, of the four faces in the three cyclic
            %   permutations [a b c] of their nodes at the point X, or
            %   their curls for OP = 'curl'. VAL is a
            %   1-by-4-by-NRBASISFUNC4F-by-3-by-3 array.
            faceNodes = [1 2 3; 1 2 4; 4 2 3; 1 3 4];
            [lambda,gradLambda] = BarycentricJacobi.coordinates(x);
            terms = BarycentricJacobi.leibnizTerms(BarycentricJacobi.multiIndices(3,1),1);
            val = zeros(1,4,nrBasisFunc4f,3,3);
            for f = 1 : 4
                for perm = 1 : 3
                    nodes = faceNodes(f,circshift(1:3,1-perm));
                    U = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(nodes(2),nodes(1:2),1,3),p,0);
                    chain = BarycentricJacobi.chainRule(nodes(3),nodes,1,3);
                    [wab,curlwab] = BarycentricJacobi.whitney(lambda,gradLambda,nodes(1),nodes(2));
                    face = zeros(nrBasisFunc4f,3);
                    indx = 0;
                    % curl-free
                    for i = 2 : p
                        if isempty(op)
                            UV = BarycentricJacobi.product(U(i+1,:),BarycentricJacobi.factor(lambda,chain,p+1-i,2*i-1),terms);
                            face(indx+(1:p+1-i),:) = UV(2:end,2:4);
                        end
                        indx = indx + p+1-i;
                    end
                    % type I
                    V = BarycentricJacobi.factor(lambda,chain,p-1,1);
                    face(indx+(1:p-1),:) = -BarycentricJacobi.timesField(V(2:end,:),wab,curlwab,op);
                    indx = indx + p-1;
                    % type II
                    for i = 2 : p
                        V = BarycentricJacobi.factor(lambda,chain,p+1-i,2*i-1);
                        if ~isempty(op)
                            face(indx+(1:p+1-i),:) = 2*BarycentricJacobi.cross(V(2:end,2:4),U(i+1,2:4));
                        else
                            face(indx+(1:p+1-i),:) = V(2:end,1)*U(i+1,2:4) - U(i+1,1)*V(2:end,2:4);
                        end
                        indx = indx + p+1-i;
                    end
                    val(1,f,:,:,perm) = face;
                end
            end
        end

        function val = elemFunctions(p,x,nrBasisFunc4e,op)
            %elemFunctions - Evaluate the cell basis or its curl
            %   VAL = elemFunctions(P,X,NRBASISFUNC4E,OP) evaluates the
            %   cell basis functions of the nodes [a b c] = [1 2 3]:
            %   grad(u_i v_j z_k), i = 2,...,P-1, j = 1,...,P-i,
            %   k = 1,...,P+1-i-j (curl-free), -v_j z_k w_12 with
            %   alpha = 1 and 2j, j = 1,...,P-2, k = 1,...,P-j-1 (type I),
            %   v_j z_k grad(u_i) (type II), and u_i v_j grad(z_k)
            %   (type III), i = 2,...,P-1, j = 1,...,P-i,
            %   k = 1,...,P+1-i-j, with alpha = 2i-1 and 2i+2j-2
            %   otherwise, at the point X, or their curls for
            %   OP = 'curl'. VAL is a 1-by-NRBASISFUNC4E-by-3 array.
            [lambda,gradLambda] = BarycentricJacobi.coordinates(x);
            terms = BarycentricJacobi.leibnizTerms(BarycentricJacobi.multiIndices(3,1),1);
            U = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(2,[1 2],1,3),p-1,0);
            chainV = BarycentricJacobi.chainRule(3,[1 2 3],1,3);
            chainW = BarycentricJacobi.chainRule(4,1:4,1,3);
            [w12,curlw12] = BarycentricJacobi.whitney(lambda,gradLambda,1,2);
            val = zeros(nrBasisFunc4e,3);
            indx = 0;
            % curl-free
            for i = 2 : p-1
                V = BarycentricJacobi.factor(lambda,chainV,p-i,2*i-1);
                for j = 1 : p-i
                    if isempty(op)
                        UVW = BarycentricJacobi.product(BarycentricJacobi.product(U(i+1,:),V(j+1,:),terms), ...
                            BarycentricJacobi.factor(lambda,chainW,p+1-i-j,2*i+2*j-2),terms);
                        val(indx+(1:p+1-i-j),:) = UVW(2:end,2:4);
                    end
                    indx = indx + p+1-i-j;
                end
            end
            % type I
            V = BarycentricJacobi.factor(lambda,chainV,p-2,1);
            for j = 1 : p-2
                VW = BarycentricJacobi.product(V(j+1,:),BarycentricJacobi.factor(lambda,chainW,p-1-j,2*j),terms);
                val(indx+(1:p-1-j),:) = -BarycentricJacobi.timesField(VW(2:end,:),w12,curlw12,op);
                indx = indx + p-1-j;
            end
            % type II
            for i = 2 : p-1
                V = BarycentricJacobi.factor(lambda,chainV,p-i,2*i-1);
                for j = 1 : p-i
                    VW = BarycentricJacobi.product(V(j+1,:),BarycentricJacobi.factor(lambda,chainW,p+1-i-j,2*i+2*j-2),terms);
                    if ~isempty(op)
                        val(indx+(1:p+1-i-j),:) = BarycentricJacobi.cross(VW(2:end,2:4),U(i+1,2:4));
                    else
                        val(indx+(1:p+1-i-j),:) = VW(2:end,1)*U(i+1,2:4);
                    end
                    indx = indx + p+1-i-j;
                end
            end
            % type III
            for i = 2 : p-1
                V = BarycentricJacobi.factor(lambda,chainV,p-i,2*i-1);
                for j = 1 : p-i
                    UV = BarycentricJacobi.product(U(i+1,:),V(j+1,:),terms);
                    W = BarycentricJacobi.factor(lambda,chainW,p+1-i-j,2*i+2*j-2);
                    if ~isempty(op)
                        val(indx+(1:p+1-i-j),:) = BarycentricJacobi.cross(UV(2:4),W(2:end,2:4));
                    else
                        val(indx+(1:p+1-i-j),:) = UV(1)*W(2:end,2:4);
                    end
                    indx = indx + p+1-i-j;
                end
            end
            val = reshape(val,1,nrBasisFunc4e,3);
        end
    end
end