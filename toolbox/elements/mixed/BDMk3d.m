classdef BDMk3d < FE
    %BDMk3d - Brezzi-Douglas-Marini element of degree k on tetrahedra
    %   OBJ = BDMk3d(DEGREE,MESH) creates the tetrahedral
    %   Brezzi-Douglas-Marini H(div)-conforming vector-valued finite
    %   element of polynomial degree DEGREE >= 1 on the tetrahedral mesh
    %   MESH. The local space consists of all vector-valued polynomials
    %   of degree at most DEGREE, and DEGREE = 1 gives the lowest-order
    %   element with three degrees of freedom per face. The degrees of
    %   freedom are the coefficients of a hierarchical face and cell
    %   basis after Beuchler, Pillwein, and Zaglmayr; those on the
    %   Neumann boundary are fixed. DEGREE < 1 raises the error
    %   fem4m:BDMk3d:invalidDegree.
    %
    %   BDMk3d functions:
    %       getBasis           - Evaluate the basis functions
    %       getDivBasis        - Evaluate the divergence of the basis functions
    %       getFaceRefBasis    - Evaluate the face reference basis
    %       getFaceDivRefBasis - Evaluate the face div reference basis
    %       getElemRefBasis    - Evaluate the element reference basis
    %       getElemDivRefBasis - Evaluate the element div ref basis
    %       p1RefBasis         - Evaluate the P1 reference basis
    %       p1GradRefBasis     - Evaluate the gradients of the P1 basis
    %
    %   BDMk3d properties:
    %       nrDoF          - Total number of degrees of freedom
    %       dof4e          - Local-to-global dofs per element
    %       pkFE           - Placeholder for the associated Pk element
    %       nrBasisFunc    - Number of basis functions per element
    %       degree         - Polynomial degree of the element
    %       fixedNodes     - Indices of the fixed degrees of freedom
    %       dof4Nb         - Dofs on the Neumann boundary
    %       freeNodes      - Indices of the free degrees of freedom
    %       nrBasis4n      - Number of nodal basis functions
    %       nrBasis4ed     - Number of edge basis functions
    %       nrBasis4f      - Number of face basis functions
    %       nrBasis4e      - Number of element bubble functions
    %       mesh           - Mesh on which the element is defined
    %       isRefElemBased - True if the basis lives on the reference element
    %       sign4e         - Face orientation signs per element
    %       perm4e         - Node permutation index of each face per element
    %
    %   See also BDMk, RTk3d, HBDMk, FEFace

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = public)
        nrDoF          % total number of degrees of freedom
        dof4e          % local-to-global dofs per element
        pkFE           % placeholder for the associated Pk element
        nrBasisFunc    % number of basis functions per element
        degree         % polynomial degree of the element
        fixedNodes     % indices of the fixed degrees of freedom
        dof4Nb         % dofs on the Neumann boundary
        freeNodes      % indices of the free degrees of freedom
        nrBasis4n      % number of nodal basis functions
        nrBasis4ed     % number of edge basis functions
        nrBasis4f      % number of face basis functions
        nrBasis4e      % number of element bubble functions
        mesh           % mesh on which the element is defined
        isRefElemBased = true % true if the basis lives on the reference element
    end

    properties(SetAccess = public)
        sign4e            % face orientation signs per element
        perm4e            % node permutation index of each face per element
    end

    properties(Access = private)
        jacobian         % transposed Jacobians of the element maps
        facePolyIndices  % type indices [i j] of the face functions
        area4e           % volumes of the elements
    end

    methods
        function obj = BDMk3d(degree,mesh)
            %BDMk3d - Create a tetrahedral BDM finite element object
            if degree < 1
                error('fem4m:BDMk3d:invalidDegree','BDM does only work for degree >= 1');
            end
            obj.degree = degree;
            obj = obj.dofCounts();
            obj.nrBasisFunc = 4*obj.nrBasis4f + obj.nrBasis4e;
            obj = obj.getBasisCoefficients(mesh);
            obj = obj.dofManager(mesh);
            obj.jacobian = reshape(mesh.jacobian,mesh.nrElems,3,3);
            obj.area4e = mesh.area4e;
            obj.mesh = mesh;
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   contravariant Piola-mapped Brezzi-Douglas-Marini basis
            %   functions at the reference point PTS_REF (1-by-3) in all
            %   elements of the mesh. VAL has the size
            %   nrElems-by-nrBasisFunc-by-3. PTS and PARTS are ignored; the
            %   callers keep the rows of the elements PARTS with
            %   FormBasis.elementRows.
            nrElems = size(obj.area4e,1);
            val = zeros(nrElems,obj.nrBasisFunc,3);
            faceValues = ones(nrElems,1).*obj.getFaceRefBasis(obj.degree,x,obj.nrBasis4f);
            for k = 1 : obj.nrBasis4f
                for curFace = 1 : 4
                    for perm = 1 : 3
                        indx = obj.perm4e(:,curFace)==perm;
                        for coordinate = 1 : 3
                            if mod(obj.facePolyIndices(k,1),2)==1
                                val(indx,4*(k-1)+curFace,coordinate) = faceValues(indx,curFace,k,coordinate,perm).*obj.sign4e(indx,curFace);
                            else
                                val(indx,4*(k-1)+curFace,coordinate) = faceValues(indx,curFace,k,coordinate,perm);
                            end
                        end
                    end
                end
            end
            val(:,4*obj.nrBasis4f+1:end,:) = ones(nrElems,1).*obj.getElemRefBasis(obj.degree,x,obj.nrBasis4e);

            % contravariant Piola transformation
            val = matMul(val,obj.jacobian);
            val = (1./(6*obj.area4e)).*val;
        end

        function val = getDivBasis(obj,~,x,~)
            %getDivBasis - Evaluate the divergence of the basis functions
            %   VAL = getDivBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   divergences of the tetrahedral Brezzi-Douglas-Marini
            %   basis functions at the reference point PTS_REF (1-by-3)
            %   in all elements of the mesh. VAL has the size
            %   nrElems-by-nrBasisFunc. PTS and PARTS are ignored; the
            %   callers keep the rows of the elements PARTS with
            %   FormBasis.elementRows.
            nrElems = size(obj.area4e,1);
            val = zeros(nrElems,obj.nrBasisFunc);
            faceValues = ones(nrElems,1).*obj.getFaceDivRefBasis(obj.degree,x,obj.nrBasis4f);
            for k = 1 : 1 % obj.nrBasis4f; (only lowest order has div)
                for curFace = 1 : 4
                    for perm = 1 : 3
                        indx = obj.perm4e(:,curFace)==perm;
                        if mod(obj.facePolyIndices(k,1),2)==1
                            val(indx,4*(k-1)+curFace) = faceValues(indx,curFace,k,1,perm).*obj.sign4e(indx,curFace);
                        else
                            val(indx,4*(k-1)+curFace) = faceValues(indx,curFace,k,1,perm);
                        end
                    end
                end
            end
            val(:,4*obj.nrBasis4f+1:end,:) = ones(nrElems,1).*obj.getElemDivRefBasis(obj.degree,x,obj.nrBasis4e);

            % contravariant Piola transformation (div scales by 1/det J)
            val = (1./(6*obj.area4e)).*val;
        end

    end

    methods (Access = private)

        function obj = getBasisCoefficients(obj,mesh)
            %getBasisCoefficients - Set face orientations and index types
            obj.sign4e = ones(mesh.nrElems,4);
            obj.perm4e = zeros(mesh.nrElems,4);
            n4e = mesh.n4e;
            % face orientation
            faceNodes = [1 2 3; 1 2 4; 4 2 3; 1 3 4];
            permutations = [1 2 3; 2 3 1; 3 1 2];
            negativepermutations = [2 1 3; 3 2 1; 1 3 2];
            for j = 1 : 4
                [~,sortIdx] = sort(n4e(:,faceNodes(j,:)),2);
                [~,perm] = ismember(sortIdx,permutations,'rows');
                [~,negperm] = ismember(sortIdx,negativepermutations,'rows');
                obj.sign4e(negperm>0,j) = -1;
                obj.perm4e(:,j) = perm + negperm;
            end
            % save indices i j for face polynomials
            obj.facePolyIndices = zeros(obj.nrBasis4f,2);
            indx = 1;
            p = obj.degree;
            % non-div free
            % Type I
            obj.facePolyIndices(indx,:) = [1,0];
            indx = indx + 1;
            % div-free
            % Type II
            for j = 1 : p
                obj.facePolyIndices(indx,:) = [1,j];
                indx = indx + 1;
            end
            % Type III
            for i = 2 : p + 1
                for j = 1 : p + 2 - i
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
            %   calls it before getBasisCoefficients and dofManager.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 0;
            obj.nrBasis4f = 1 + obj.degree + (obj.degree+1)*obj.degree/2;
            obj.nrBasis4e = obj.degree*(obj.degree-1)/2 + 2*(obj.degree+1)*obj.degree*(obj.degree-1)/6 ...
                + (obj.degree-1) + obj.degree*(obj.degree-1)/2 + (obj.degree+1)*obj.degree*(obj.degree-1)/6;
        end

        function obj = dofManager(obj,mesh)
            %dofManager - Assemble the degrees of freedom of the element
            %   OBJ = dofManager(OBJ,MESH) enumerates the face and
            %   interior degrees of freedom with FE.dofManager, fixes
            %   those on the Neumann boundary (dof4Nb), and frees all
            %   others.
            obj.mesh = mesh;
            obj = dofManager@FE(obj);
            obj.fixedNodes = obj.dof4Nb;
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
        end

    end

    methods(Static)

        % S. Beuchler, V. Pillwein, S. Zaglmayr, Sparsity optimized high order
        % finite element functions for H(div) on simplices,
        % Numer. Math. 122 (2012), no. 2, 197-225.

        %% Face Basis helper functions
        function val = getFaceRefBasis(p,x,nrBasisFunc4f)
            %getFaceRefBasis - Evaluate the face reference basis
            %   VAL = getFaceRefBasis(P,X,NRBASISFUNC4F) evaluates the
            %   NRBASISFUNC4F face-based Brezzi-Douglas-Marini basis
            %   functions of degree P on the reference tetrahedron at the
            %   point X (1-by-3) for all four faces and the three cyclic
            %   permutations of their nodes, the RTk3d face functions of
            %   degree P+1. VAL has the size 1-by-4-by-NRBASISFUNC4F-by-3-by-3
            %   (face, function, component, permutation).
            val = RTk3d.getFaceRefBasis(p+1,x,nrBasisFunc4f);
        end

        function val = getFaceDivRefBasis(~,x,nrBasisFunc4f)
            %getFaceDivRefBasis - Evaluate the face div reference basis
            %   VAL = getFaceDivRefBasis(P,X,NRBASISFUNC4F) evaluates the
            %   divergences of the NRBASISFUNC4F face-based
            %   Brezzi-Douglas-Marini basis functions at the point X
            %   (1-by-3) on the reference tetrahedron. Only the
            %   lowest-order (Whitney) face functions have a nonzero
            %   divergence, so the degree P is ignored. VAL has the size
            %   1-by-4-by-NRBASISFUNC4F-by-1-by-3 (face, function, 1,
            %   permutation).
            val = RTk3d.getFaceDivRefBasis([],x,nrBasisFunc4f);
        end

        %% Elem Basis helper functions
        function val = getElemRefBasis(p,x,nrBasisFunc4e)
            %getElemRefBasis - Evaluate the element reference basis
            %   VAL = getElemRefBasis(P,X,NRBASISFUNC4E) evaluates the
            %   NRBASISFUNC4E interior bubble basis functions of the
            %   tetrahedral Brezzi-Douglas-Marini element of degree P at
            %   the point X (1-by-3) on the reference tetrahedron. VAL has
            %   the size 1-by-NRBASISFUNC4E-by-3.
            val = BDMk3d.elemFunctions(p,x,nrBasisFunc4e,'');
        end

        function val = getElemDivRefBasis(p,x,nrBasisFunc4e)
            %getElemDivRefBasis - Evaluate the element div ref basis
            %   VAL = getElemDivRefBasis(P,X,NRBASISFUNC4E) evaluates the
            %   divergences of the NRBASISFUNC4E interior bubble basis
            %   functions of the tetrahedral Brezzi-Douglas-Marini element
            %   of degree P at the point X (1-by-3) on the reference
            %   tetrahedron. VAL has the size 1-by-NRBASISFUNC4E.
            val = BDMk3d.elemFunctions(p,x,nrBasisFunc4e,'div');
        end

        %% P1 basis functions
        function val = p1RefBasis(x)
            %p1RefBasis - Evaluate the P1 reference basis
            %   VAL = p1RefBasis(X) evaluates the four linear nodal basis
            %   functions of the reference tetrahedron at the points X
            %   (N-by-3). VAL has the size N-by-4.
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
            %   tetrahedron at the points X (N-by-3). VAL has the size
            %   N-by-4-by-3.
            val = zeros(size(x,1),4,3);
            val(:,1,:) = -1;
            val(:,2,1) = 1;
            val(:,3,2) = 1;
            val(:,4,3) = 1;
        end

    end

    methods (Static, Access = private)

        function val = elemFunctions(p,x,nrBasisFunc4e,op)
            %elemFunctions - Evaluate the cell basis or its divergence
            %   VAL = elemFunctions(P,X,NRBASISFUNC4E,OP) evaluates the
            %   cell basis functions at the point X (1-by-3), or their
            %   divergences for OP = 'div': the div-free cell functions of
            %   RTk3d of degree P+1, 4 w_k phi_123 with alpha = 4,
            %   k = 1,...,P-1 (type I), 2 v_j grad(w_k) x w_12 with alpha = 3
            %   and 2j+2, j = 1,...,P-1, k = 1,...,P-j (type II), and
            %   w_k grad(u_i) x grad(v_j), i = 2,...,P, j = 1,...,P+1-i,
            %   k = 1,...,P+2-i-j, with alpha = 2i-1 and 2i+2j-2 (type
            %   III), where u_i, v_j, w_k, w_12, and phi_123 are as in
            %   RTk3d. VAL has the size 1-by-NRBASISFUNC4E-by-3, or
            %   1-by-NRBASISFUNC4E for the divergences.
            [lambda,gradLambda] = BarycentricJacobi.coordinates(x);
            U = BarycentricJacobi.factor(lambda,BarycentricJacobi.chainRule(2,[1 2],1,3),p,0);
            chainV = BarycentricJacobi.chainRule(3,[1 2 3],1,3);
            chainW = BarycentricJacobi.chainRule(4,1:4,1,3);
            [w12,curlw12] = BarycentricJacobi.whitney(lambda,gradLambda,1,2);
            [phi123,divphi123] = BarycentricJacobi.whitneyFace(lambda,gradLambda,1,2,3);
            free = RTk3d.divFreeElemFunctions(p+1,x);
            indx = size(free,1);
            div = strcmp(op,'div');
            if div
                val = zeros(nrBasisFunc4e,1);
            else
                val = zeros(nrBasisFunc4e,3);
                val(1:indx,:) = free;
            end
            % type I
            W = BarycentricJacobi.factor(lambda,chainW,p-1,4);
            val(indx+(1:p-1),:) = 4*BarycentricJacobi.timesField(W(2:end,:),phi123,divphi123,op);
            indx = indx + p-1;
            % type II, div(v grad(w) x w_12)
            % = grad(v).(grad(w) x w_12) - v grad(w).curl(w_12)
            V = BarycentricJacobi.factor(lambda,chainV,p-1,3);
            for j = 1 : p-1
                W = BarycentricJacobi.factor(lambda,chainW,p-j,2*j+2);
                gradWxw12 = BarycentricJacobi.cross(W(2:end,2:4),w12);
                if div
                    val(indx+(1:p-j),:) = 2*(gradWxw12*V(j+1,2:4)' - V(j+1,1)*W(2:end,2:4)*curlw12');
                else
                    val(indx+(1:p-j),:) = 2*V(j+1,1)*gradWxw12;
                end
                indx = indx + p-j;
            end
            % type III
            for i = 2 : p
                V = BarycentricJacobi.factor(lambda,chainV,p+1-i,2*i-1);
                for j = 1 : p+1-i
                    W = BarycentricJacobi.factor(lambda,chainW,p+2-i-j,2*i+2*j-2);
                    gradUxgradV = BarycentricJacobi.cross(U(i+1,2:4),V(j+1,2:4));
                    if div
                        val(indx+(1:p+2-i-j),:) = W(2:end,2:4)*gradUxgradV';
                    else
                        val(indx+(1:p+2-i-j),:) = W(2:end,1)*gradUxgradV;
                    end
                    indx = indx + p+2-i-j;
                end
            end
            val = reshape(val,1,nrBasisFunc4e,size(val,2));
        end
    end
end