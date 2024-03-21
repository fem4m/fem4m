classdef DGDk3d < FE
    %DGDk3d - Discontinuous Dubiner element of degree k on tetrahedra
    %   OBJ = DGDk3d(DEGREE,MESH) creates the discontinuous
    %   element of total polynomial degree DEGREE >= 0 (a scalar) on
    %   the tetrahedral mesh MESH. Its element basis is the
    %   hierarchical Dubiner basis of Jacobi polynomials on the
    %   reference tetrahedron.
    %
    %   DGDk3d functions:
    %       getBasis             - Evaluate the basis functions
    %       getGradBasis         - Evaluate the gradients of the basis
    %       getD2Basis           - Evaluate the Hessians of the basis
    %       getD3Basis           - Evaluate third derivatives of the basis
    %       getD4Basis           - Evaluate fourth derivatives of the basis
    %       getElemBasis         - Evaluate the Dubiner basis on a tetrahedron
    %       getElemBasisGrad     - Evaluate the gradient of the Dubiner basis
    %       getElemBasisD2       - Evaluate the Hessians of the Dubiner basis
    %       getElemBasisD3       - Evaluate third derivatives of the basis
    %       getElemBasisD4       - Evaluate fourth derivatives of the basis
    %       basisN               - Evaluate the P1 basis functions
    %
    %   DGDk3d properties:
    %       nrDoF          - Total number of degrees of freedom
    %       dof4e          - Local-to-global dofs per element
    %       freeNodes      - Indices of the free degrees of freedom
    %       fixedNodes     - Indices of the fixed degrees of freedom
    %       degree         - Polynomial degree of the element
    %       nrBasis4n      - Number of nodal basis functions
    %       nrBasis4ed     - Number of edge basis functions
    %       nrBasis4f      - Number of face basis functions
    %       nrBasis4e      - Number of element basis functions (all)
    %       nrBasisFunc    - Number of basis functions per element
    %       transformation - Inverse Jacobians of the element maps
    %       mesh           - Mesh on which the element is defined
    %       isRefElemBased - Basis lives on the reference element
    %
    %   See also Pk3d, DGDk, DGDkVector, FE

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = public)
        nrDoF            % total number of degrees of freedom
        dof4e            % local-to-global dofs per element
        freeNodes        % indices of the free degrees of freedom
        fixedNodes       % indices of the fixed degrees of freedom
        degree           % polynomial degree of the element
        nrBasis4n        % number of nodal basis functions
        nrBasis4ed       % number of edge basis functions
        nrBasis4f        % number of face basis functions
        nrBasis4e        % number of element basis functions (all)
        nrBasisFunc      % number of basis functions per element
        transformation   % inverse Jacobians of the element maps
        mesh             % mesh on which the element is defined
        isRefElemBased = true % basis lives on the reference element
    end

    methods

        function obj = DGDk3d(degree,mesh)
            %DGDk3d - Create a DGDk3d finite element object
            if degree < 0
                error('fem4m:DGDk3d:invalidDegree','DGDk3d degree must be at least 0');
            end
            obj.degree = degree;
            obj = obj.dofCounts();
            obj.nrBasisFunc = 4*obj.nrBasis4n + 6*obj.nrBasis4ed + 4*obj.nrBasis4f + obj.nrBasis4e;
            obj = obj.dofManager(mesh);
            obj.transformation = mesh.transformation;
            obj.mesh = mesh;
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,X,PARTS) evaluates the DGDk3d
            %   basis functions at the reference point X (1-by-3)
            %   and replicates the result for all elements, so VAL
            %   is nrElems-by-nrBasisFunc. PTS and PARTS are
            %   ignored.

            val = ones(size(obj.transformation,1),1).*permute(obj.getElemBasis(obj.degree,x,obj.nrBasis4e),[3 1 2]);

        end

        function val = getGradBasis(obj,~,x,~)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,X,PARTS) evaluates the
            %   gradients of the DGDk3d basis at the reference
            %   point X (1-by-3), mapped to the physical elements
            %   by the saved transformation. VAL is
            %   nrElems-by-nrBasisFunc-by-3. PTS and PARTS are
            %   ignored.
            grad = obj.getElemBasisGrad(obj.degree,x,obj.nrBasis4e);
            val = FE.mapDerivatives(obj.transformation,permute(grad,[3 1 2]),1);
        end

        function val = getD2Basis(obj,~,x,~)
            %getD2Basis - Evaluate the Hessians of the basis
            %   VAL = getD2Basis(OBJ,PTS,X,PARTS) evaluates the
            %   Hessians [dxx,dyy,dzz,dxy,dxz,dyz] of the DGDk3d
            %   basis at the reference point X (1-by-3), mapped to
            %   the physical elements. PTS and PARTS are ignored.
            d2 = obj.getElemBasisD2(obj.degree,x,obj.nrBasis4e);
            val = FE.mapDerivatives(obj.transformation,permute(d2,[3 1 2]),2);
        end

        function val = getD3Basis(obj,~,x,~)
            %getD3Basis - Evaluate third derivatives of the basis
            %   VAL = getD3Basis(OBJ,PTS,X,PARTS) evaluates the ten
            %   third derivatives [xxx,yyy,zzz,xxy,xxz,xyy,yyz,
            %   xzz,yzz,xyz] of the DGDk3d basis at the reference
            %   point X (1-by-3), mapped to the physical elements.
            %   PTS and PARTS are ignored.
            d3 = obj.getElemBasisD3(obj.degree,x,obj.nrBasis4e);
            val = FE.mapDerivatives(obj.transformation,permute(d3,[3 1 2]),3);
        end

        function val = getD4Basis(obj,~,x,~)
            %getD4Basis - Evaluate fourth derivatives of the basis
            %   VAL = getD4Basis(OBJ,PTS,X,PARTS) evaluates the 15
            %   fourth derivatives [xxxx,yyyy,zzzz,xxxy,xxxz,
            %   xyyy,yyyz,xzzz,yzzz,xxyy,xxzz,yyzz,xxyz,xyyz,
            %   xyzz] of the DGDk3d basis at the reference point
            %   X (1-by-3), mapped to the physical elements. PTS
            %   and PARTS are ignored.
            d4 = obj.getElemBasisD4(obj.degree,x,obj.nrBasis4e);
            val = FE.mapDerivatives(obj.transformation,permute(d4,[3 1 2]),4);
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
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 0;
            obj.nrBasis4f = 0;
            if obj.degree == 0
                obj.nrBasis4e = 1;
            else
                obj.nrBasis4e = 4 + 6*(obj.degree-1) + 4*(obj.degree-1)*(obj.degree-2)/2 + (obj.degree-1)*(obj.degree-2)*(obj.degree-3)/6;
            end
        end

        function obj = dofManager(obj,mesh)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ,MESH) enumerates the element-local
            %   degrees of freedom with FE.dofManager and keeps all of
            %   them free.
            obj.mesh = mesh;
            obj = dofManager@FE(obj);
            obj.fixedNodes = [];
            obj.freeNodes = 1:obj.nrDoF;
        end

    end

    methods(Static)

        % Dubiner basis functions

        function val = getElemBasis(p,x,nrBasisFunc4e)
            %getElemBasis - Evaluate the Dubiner basis on a tetrahedron
            %   VAL = getElemBasis(P,X,NRBASISFUNC4E) evaluates
            %   the Dubiner basis of total degree P at the
            %   reference points X (N-by-3) in local coordinates.
            %   NRBASISFUNC4E = (P+1)(P+2)(P+3)/6 is the number of
            %   basis functions, and VAL is NRBASISFUNC4E-by-1-by-N.
            val = DGDk3d.elementDerivatives(p,x,nrBasisFunc4e,[0 0 0]);
        end

        function val = getElemBasisGrad(p,x,nrBasisFunc4e)
            %getElemBasisGrad - Evaluate the gradient of the Dubiner basis
            %   VAL = getElemBasisGrad(P,X,NRBASISFUNC4E)
            %   evaluates the gradients [dx,dy,dz] of the
            %   Dubiner basis of total degree P at the reference
            %   points X (N-by-3) in local coordinates. VAL is
            %   NRBASISFUNC4E-by-3-by-N.
            val = DGDk3d.elementDerivatives(p,x,nrBasisFunc4e,[1 0 0; 0 1 0; 0 0 1]);
        end

        function val = getElemBasisD2(p,x,nrBasisFunc4e)
            %getElemBasisD2 - Evaluate the Hessians of the Dubiner basis
            %   VAL = getElemBasisD2(P,X,NRBASISFUNC4E) evaluates
            %   the Hessians [dxx,dyy,dzz,dxy,dxz,dyz] of the
            %   Dubiner basis of total degree P at the reference
            %   points X (N-by-3) in local coordinates. VAL is
            %   NRBASISFUNC4E-by-6-by-N.
            val = DGDk3d.elementDerivatives(p,x,nrBasisFunc4e, ...
                [2 0 0; 0 2 0; 0 0 2; 1 1 0; 1 0 1; 0 1 1]);
        end

        function val = getElemBasisD3(p,x,nrBasisFunc4e)
            %getElemBasisD3 - Evaluate third derivatives of the basis
            %   VAL = getElemBasisD3(P,X,NRBASISFUNC4E) evaluates
            %   the third derivatives [xxx,yyy,zzz,xxy,xxz,xyy,
            %   yyz,xzz,yzz,xyz] of the Dubiner basis of total
            %   degree P at the reference points X (N-by-3) in
            %   local coordinates. VAL is NRBASISFUNC4E-by-10-by-N.
            val = DGDk3d.elementDerivatives(p,x,nrBasisFunc4e, ...
                [3 0 0; 0 3 0; 0 0 3; 2 1 0; 2 0 1; 1 2 0; 0 2 1; 1 0 2; 0 1 2; 1 1 1]);
        end

        function val = getElemBasisD4(p,x,nrBasisFunc4e)
            %getElemBasisD4 - Evaluate fourth derivatives of the basis
            %   VAL = getElemBasisD4(P,X,NRBASISFUNC4E) evaluates
            %   the fourth derivatives [xxxx,yyyy,zzzz,xxxy,
            %   xxxz,xyyy,yyyz,xzzz,yzzz,xxyy,xxzz,yyzz,xxyz,
            %   xyyz,xyzz] of the Dubiner basis of total degree P
            %   at the reference points X (N-by-3) in local
            %   coordinates. VAL is NRBASISFUNC4E-by-15-by-N.
            val = DGDk3d.elementDerivatives(p,x,nrBasisFunc4e, ...
                [4 0 0; 0 4 0; 0 0 4; 3 1 0; 3 0 1; 1 3 0; 0 3 1; 1 0 3; 0 1 3; ...
                2 2 0; 2 0 2; 0 2 2; 2 1 1; 1 2 1; 1 1 2]);
        end

        %% P1 basis functions
        function val = basisN(x)
            %basisN - Evaluate the P1 basis functions
            %   VAL = basisN(X) evaluates the four linear basis
            %   functions 1-x-y-z, x, y and z at the reference
            %   point X (1-by-3) in local coordinates. VAL is
            %   4-by-1.
            val = zeros(4,1);
            val(1) = 1 - x(1) - x(2) - x(3);
            val(2) = x(1);
            val(3) = x(2);
            val(4) = x(3);
        end

    end

    methods (Static, Access = private)

        function val = elementDerivatives(p,x,nrBasisFunc4e,D)
            %elementDerivatives - Evaluate derivatives of the Dubiner basis
            %   VAL = elementDerivatives(P,X,NRBASISFUNC4E,D) evaluates
            %   the derivatives d^(a+b+c)/(dx^a dy^b dz^c), [a b c] =
            %   D(l,:), of the Dubiner basis phi_ijk = Q_i(x,1-y-z)
            %   Q_j^(2i+1,0)(y,1-z) P_k^(2i+2j+2,0)(2z-1) of total degree
            %   P at the reference points X, where Q_n^(alpha,beta)(x,s)
            %   = s^n P_n^(alpha,beta)(2x/s-1) are the scaled Jacobi
            %   polynomials. VAL is NRBASISFUNC4E-by-size(D,1)-by-N for
            %   the N-by-3 points X.
            m = max(sum(D,2));
            nrPoints = size(x,1);
            z = x(:,3);
            y = x(:,2);
            x = x(:,1);
            % derivative pairs (a,c) with a+c <= m, numbered by pair4ac
            pair4ac = zeros(m+1);
            pairs = zeros(0,2);
            for a = 0 : m
                for c = 0 : m-a
                    pairs(end+1,:) = [a c]; %#ok<AGROW>
                    pair4ac(a+1,c+1) = size(pairs,1);
                end
            end
            % Leibniz rule, where 1-y-z and 1-z have the derivative -1:
            % each term of the derivative D(l,:) is coeff times the
            % derivative pair iA of Q_i in (x,s), the pair iB of Q_j in
            % (y,s), and d^c3/dz^c3 of P_k
            fact = factorial(0:m);
            coeff = []; iA = []; iB = []; c3 = []; l = [];
            for d = 1 : size(D,1)
                a = D(d,1); b = D(d,2); c = D(d,3);
                for b1 = 0 : b
                    for c1 = 0 : c
                        for c2 = 0 : c-c1
                            coeff(end+1,1) = (-1)^(b1+c1+c2)*fact(b+1)/(fact(b1+1)*fact(b-b1+1)) ...
                                *fact(c+1)/(fact(c1+1)*fact(c2+1)*fact(c-c1-c2+1)); %#ok<AGROW>
                            iA(end+1,1) = pair4ac(a+1,b1+c1+1); %#ok<AGROW>
                            iB(end+1,1) = pair4ac(b-b1+1,c2+1); %#ok<AGROW>
                            c3(end+1,1) = c-c1-c2; %#ok<AGROW>
                            l(end+1,1) = d; %#ok<AGROW>
                        end
                    end
                end
            end
            nrTerms = numel(coeff);
            pairsA = unique(iA)';
            pairsB = unique(iB)';
            % sum of the terms with derivative c3 of P_k for each D(l,:)
            E3 = double(c3 == 0:m);
            El = double(l == 1:size(D,1));
            % QA(i+1,n,:) = d^a/dx^a d^c/ds^c Q_i(x,s), s = 1-y-z, (a,c) = pairs(n,:)
            QA = zeros(p+1,size(pairs,1),nrPoints);
            for n = pairsA
                QA(:,n,:) = reshape(JacobiPolynomials.scaledPoly(x,1-y-z,p,0,0,pairs(n,1),pairs(n,2)),p+1,1,nrPoints);
            end
            val = zeros(nrBasisFunc4e,size(D,1),nrPoints);
            indx = 1;
            for i = 0 : p
                % QB(j+1,n,:) = d^a/dy^a d^c/ds^c Q_j^(2i+1,0)(y,s), s = 1-z
                QB = zeros(p-i+1,size(pairs,1),nrPoints);
                for n = pairsB
                    QB(:,n,:) = reshape(JacobiPolynomials.scaledPoly(y,1-z,p-i,2*i+1,0,pairs(n,1),pairs(n,2)),p-i+1,1,nrPoints);
                end
                for j = 0 : p-i
                    % R(k+1,c+1,:) = d^c/dz^c P_k^(2i+2j+2,0)(2z-1)
                    R = zeros(p-i-j+1,m+1,nrPoints);
                    for c = 0 : m
                        R(:,c+1,:) = reshape(2^c*JacobiPolynomials.poly(2*z-1,p-i-j,2*i+2*j+2,0,c),p-i-j+1,1,nrPoints);
                    end
                    w = coeff.*reshape(QA(i+1,iA,:),nrTerms,nrPoints).*reshape(QB(j+1,iB,:),nrTerms,nrPoints);
                    M = pagemtimes(E3,'transpose',reshape(w,nrTerms,1,nrPoints).*El,'none');
                    val(indx+(0:p-i-j),:,:) = pagemtimes(R,M);
                    indx = indx + p-i-j+1;
                end
            end
        end
    end
end