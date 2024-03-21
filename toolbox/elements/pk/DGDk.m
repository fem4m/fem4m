classdef DGDk < FE
    %DGDk - Discontinuous Dubiner element of total degree k
    %   OBJ = DGDk(DEGREE,MESH) creates the discontinuous element of
    %   total degree DEGREE >= 0 on the triangular mesh MESH, a scalar or
    %   one degree per element. Its element basis is the hierarchical
    %   Dubiner basis of Jacobi polynomials on the reference triangle.
    %   With one degree per element, the basis functions of higher
    %   degree than the element degree are unused (nonDoFNodes) and not
    %   free.
    %
    %   DGDk functions:
    %       getBasis             - Evaluate the basis functions
    %       getGradBasis         - Evaluate the gradients of the basis
    %       getD2Basis           - Evaluate the Hessians of the basis
    %       getD3Basis           - Evaluate third derivatives of the basis
    %       getD4Basis           - Evaluate fourth derivatives of the basis
    %       getElementBasis      - Evaluate the Dubiner element basis
    %       getGradElementBasis  - Evaluate the gradient of the Dubiner basis
    %       getD2ElementBasis    - Evaluate the Hessian of the Dubiner basis
    %       getD3ElementBasis    - Evaluate third Dubiner basis derivatives
    %       getD4ElementBasis    - Evaluate fourth Dubiner basis derivatives
    %
    %   DGDk properties:
    %       mesh           - Mesh on which the element is defined
    %       freeNodes      - Indices of the free degrees of freedom
    %       fixedNodes     - Indices of the fixed degrees of freedom
    %       degree4e       - Polynomial degree of each element
    %       degree         - Maximum polynomial degree of the mesh
    %       nrBasisFunc    - Number of basis functions per element
    %       dof4e          - Local-to-global dofs per element
    %       nrDoF          - Number of used degrees of freedom
    %       allnrDoF       - Total number of degrees of freedom
    %       nonDoFNodes    - Unused dofs of elements of lower degree
    %       degree4DoF     - Polynomial degree of each local dof
    %       nrBasis4n      - Number of nodal basis functions
    %       nrBasis4ed     - Number of edge basis functions
    %       nrBasis4e      - Number of element basis functions (all)
    %       isRefElemBased - Basis lives on the reference element
    %
    %   See also HPk, Pk, DGDk3d, DGDkVector, DGDkMatrix, Bk

    % Copyright (c) 2026 by the fem4m authors

    properties
        mesh              % mesh on which the element is defined
        freeNodes         % indices of the free degrees of freedom
        fixedNodes        % indices of the fixed degrees of freedom
        degree4e          % polynomial degree of each element
        degree            % maximum polynomial degree of the mesh
        nrBasisFunc       % number of basis functions per element
        dof4e             % local-to-global dofs per element
        nrDoF             % number of used degrees of freedom
        allnrDoF          % total number of degrees of freedom
        nonDoFNodes       % unused dofs of elements of lower degree
        degree4DoF        % polynomial degree of each local dof
        nrBasis4n         % number of nodal basis functions
        nrBasis4ed        % number of edge basis functions
        nrBasis4e         % number of element basis functions (all)
        isRefElemBased = true % basis lives on the reference element
    end

    methods

        function obj = DGDk(degree,mesh)
            %DGDk - Create a DGDk finite element object
            if any(degree(:) < 0)
                error('fem4m:DGDk:invalidDegree','DGDk degree must be at least 0');
            end
            obj.degree4e = degree;
            obj.degree = max(degree);
            if isscalar(obj.degree4e)
                obj.degree4e = obj.degree4e*ones(mesh.nrElems,1);
            end
            obj.nrBasisFunc = 3 + 3*(obj.degree-1) + (obj.degree-2)*(obj.degree-1)/2;
            obj.mesh = mesh;
            obj = obj.dofCounts();
            obj = obj.dofManager;
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,X,PARTS) evaluates the DGDk
            %   basis functions at the reference point X (1-by-2)
            %   and replicates the result for all elements, so VAL
            %   is nrElems-by-nrBasisFunc. PTS and PARTS are
            %   ignored.
            val = ones(1,obj.nrBasisFunc);
            val(1,:) = obj.getElementBasis(obj.degree,x);
            val = ones(obj.mesh.nrElems,1).*val;
        end

        function val = getGradBasis(obj,~,x,~)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,X,PARTS) evaluates the
            %   gradients of the DGDk basis functions at the
            %   reference point X (1-by-2), mapped to the physical
            %   elements by the mesh transformation. VAL is
            %   nrElems-by-nrBasisFunc-by-2 (-by-3 on a surface
            %   mesh). PTS and PARTS are ignored.
            val = zeros(1,obj.nrBasisFunc,2);
            val(1,:,:) = obj.getGradElementBasis(obj.degree,x);
            c = obj.mesh.transformation(x);
            val = matMul(val,permute(c,[1 3 2]));
        end

        function val = getD2Basis(obj,~,x,~)
            %getD2Basis - Evaluate the Hessians of the basis
            %   VAL = getD2Basis(OBJ,PTS,X,PARTS) evaluates the
            %   Hessian [dxx,dyy,dxy] of the DGDk basis functions at
            %   the reference point X (1-by-2), mapped to the
            %   physical elements. VAL is nrElems-by-nrBasisFunc-by-3
            %   (-by-6 with [dxx,dyy,dzz,dxy,dxz,dyz] on a surface mesh,
            %   without a curvature correction for curved surfaces).
            %   PTS and PARTS are ignored.
            d2 = zeros(1,obj.nrBasisFunc,3);
            d2(1,:,:) = obj.getD2ElementBasis(obj.degree,x);
            val = FE.mapDerivatives(obj.mesh.transformation(x),d2,2);
        end

        function val = getD3Basis(obj,~,x,~)
            %getD3Basis - Evaluate third derivatives of the basis
            %   VAL = getD3Basis(OBJ,PTS,X,PARTS) evaluates the
            %   third derivatives [u_xxx,u_yyy,u_xxy,u_xyy] of the
            %   DGDk basis functions at the reference point X
            %   (1-by-2), mapped to the physical elements. VAL is
            %   nrElems-by-nrBasisFunc-by-4 (-by-10 on a surface
            %   mesh, without a curvature correction). PTS and PARTS
            %   are ignored.
            d3 = zeros(1,obj.nrBasisFunc,4);
            d3(1,:,:) = obj.getD3ElementBasis(obj.degree,x);
            val = FE.mapDerivatives(obj.mesh.transformation(x),d3,3);
        end

        function val = getD4Basis(obj,~,x,~)
            %getD4Basis - Evaluate fourth derivatives of the basis
            %   VAL = getD4Basis(OBJ,PTS,X,PARTS) evaluates the
            %   fourth derivatives [u_xxxx,u_yyyy,u_xxyy,u_xxxy,
            %   u_xyyy] of the DGDk basis functions at the reference
            %   point X (1-by-2), mapped to the physical elements.
            %   VAL is nrElems-by-nrBasisFunc-by-5 (-by-15 on a
            %   surface mesh, without a curvature correction). PTS
            %   and PARTS are ignored.
            d4 = zeros(1,obj.nrBasisFunc,5);
            d4(1,:,:) = obj.getD4ElementBasis(obj.degree,x);
            val = FE.mapDerivatives(obj.mesh.transformation(x),d4,4);
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
            obj.nrBasis4ed = 0;
            obj.nrBasis4e = obj.nrBasisFunc;
        end

        function obj = dofManager(obj)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ) numbers the element-local
            %   degrees of freedom with FE.dofManager, records the
            %   polynomial degree of each local dof, and determines
            %   the free degrees of freedom of the elementwise
            %   degrees degree4e.
            obj = dofManager@FE(obj);
            % degrees
            obj.degree4DoF = zeros(1,obj.nrBasisFunc);
            indx = 1;
            for i = 0 : obj.degree
                for j = 0 : obj.degree - i
                    obj.degree4DoF(indx) = i+j;
                    indx = indx + 1;
                end
            end
            % free dofs
            obj.allnrDoF = obj.nrDoF;
            obj.fixedNodes = [];
            nonDoFNodes = [];
            for curDegree = 0 : obj.degree
                nonDoFNodes = [nonDoFNodes;reshape(obj.dof4e(obj.degree4e==curDegree,obj.degree4DoF>curDegree),[],1)]; %#ok<AGROW>
            end
            nonDoFNodes = unique(nonDoFNodes);
            obj.nonDoFNodes = nonDoFNodes;
            obj.nrDoF = obj.allnrDoF - numel(nonDoFNodes);
            obj.freeNodes = setdiff(1:obj.allnrDoF,nonDoFNodes);
        end

    end

    methods (Static)
        %% Elem basis helper functions
        function val = getElementBasis(p,x)
            %getElementBasis - Evaluate the Dubiner element basis
            %   VAL = getElementBasis(P,X) evaluates the Dubiner
            %   basis of total degree P at the reference points X
            %   (N-by-2). VAL is N-by-(P+1)(P+2)/2.
            val = permute(DGDk.elementDerivatives(p,x,0),[3 1 2]);
        end

        function val = getGradElementBasis(p,x)
            %getGradElementBasis - Evaluate the gradient of the Dubiner basis
            %   VAL = getGradElementBasis(P,X) evaluates the
            %   gradients of the Dubiner basis of total degree P
            %   at the reference points X (N-by-2). VAL is
            %   (P+1)(P+2)/2-by-2-by-N.
            val = DGDk.elementDerivatives(p,x,1);
        end

        function val = getD2ElementBasis(p,x)
            %getD2ElementBasis - Evaluate the Hessian of the Dubiner basis
            %   VAL = getD2ElementBasis(P,X) evaluates the
            %   Hessian [dxx,dyy,dxy] of the Dubiner basis of
            %   total degree P at the reference points X (N-by-2).
            %   VAL is (P+1)(P+2)/2-by-3-by-N.
            val = DGDk.elementDerivatives(p,x,2);
            val = val(:,[1 3 2],:);
        end

        function val = getD3ElementBasis(p,x)
            %getD3ElementBasis - Evaluate third Dubiner basis derivatives
            %   VAL = getD3ElementBasis(P,X) evaluates the third
            %   derivatives [u_xxx,u_yyy,u_xxy,u_xyy] of the
            %   Dubiner basis of total degree P at the reference
            %   points X (N-by-2). VAL is (P+1)(P+2)/2-by-4-by-N.
            val = DGDk.elementDerivatives(p,x,3);
            val = val(:,[1 4 2 3],:);
        end

        function val = getD4ElementBasis(p,x)
            %getD4ElementBasis - Evaluate fourth Dubiner basis derivatives
            %   VAL = getD4ElementBasis(P,X) evaluates the fourth
            %   derivatives [u_xxxx,u_yyyy,u_xxyy,u_xxxy,u_xyyy]
            %   of the Dubiner basis of total degree P at the
            %   reference points X (N-by-2). VAL is
            %   (P+1)(P+2)/2-by-5-by-N.
            val = DGDk.elementDerivatives(p,x,4);
            val = val(:,[1 5 3 2 4],:);
        end
    end

    methods (Static, Access = private)

        function val = elementDerivatives(p,x,m)
            %elementDerivatives - Evaluate derivatives of order M
            %   VAL = elementDerivatives(P,X,M) evaluates the
            %   derivatives d^M/(dx^(M-b) dy^b), b = 0,...,M, of the
            %   Dubiner basis phi_ij = Q_i(x,1-y) P_j^(2i+1,0)(2y-1)
            %   of total degree P at the reference points X, where
            %   Q_i(x,s) = s^i P_i(2x/s-1) is the scaled Legendre
            %   polynomial. VAL is (P+1)(P+2)/2-by-(M+1)-by-N for the
            %   N-by-2 points X.
            y = x(:,2);
            x = x(:,1);
            % Q{a+1,c+1} = d^a/dx^a d^c/ds^c Q with s = 1-y
            Q = cell(m+1);
            for a = 0 : m
                for c = 0 : m-a
                    Q{a+1,c+1} = JacobiPolynomials.scaledPoly(x,1-y,p,0,0,a,c);
                end
            end
            val = zeros((p^2+3*p+2)/2,m+1,numel(x));
            indx = 1;
            for i = 0 : p
                % K{d+1} = d^d/dy^d P_j^(2i+1,0)(2y-1) for j = 0,...,p-i
                K = cell(m+1,1);
                for d = 0 : m
                    K{d+1} = 2^d*JacobiPolynomials.poly(2*y-1,p-i,2*i+1,0,d);
                end
                j = indx + (0:p-i);
                % Leibniz rule in y with ds/dy = -1
                for b = 0 : m
                    for c = 0 : b
                        val(j,b+1,:) = val(j,b+1,:) + reshape((-1)^c*nchoosek(b,c)*Q{m-b+1,c+1}(i+1,:).*K{b-c+1},[],1,numel(x));
                    end
                end
                indx = indx + p-i+1;
            end
        end
    end
end