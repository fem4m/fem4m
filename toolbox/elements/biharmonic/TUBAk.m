classdef TUBAk < FE
    %TUBAk - Argyris-type H2 element of degree k with additional dofs
    %   OBJ = TUBAk(DEGREE,MESH) creates an Argyris-type H2-conforming
    %   finite element of polynomial degree DEGREE on the triangular
    %   mesh MESH, using nodal values, gradients and Hessians together
    %   with edge and element-interior degrees of freedom.
    %   C1-conforming elements on triangles require DEGREE >= 5;
    %   TUBAk(5,MESH) spans the space of Argyris(5,MESH). The basis
    %   coefficients depend on the element geometry
    %   (isAffineEquivalent is false).
    %
    %   Each edge carries 2*DEGREE-9 degrees of freedom at the N-1
    %   points k/N, k = 1,...,N-1, of the edge, with N = DEGREE-3 for
    %   odd and N = DEGREE-2 for even DEGREE: for odd DEGREE the normal
    %   derivative along normal4ed at the edge midpoint, for even DEGREE
    %   the value at the edge midpoint and the normal derivatives at
    %   the two neighbouring points, and the x- and y-derivatives at
    %   all other points. Each element carries (DEGREE-4)*(DEGREE-5)/2
    %   interior degrees of freedom, the values at the interior points
    %   of the lattice with spacing 1/(DEGREE-3) on the reference
    %   triangle. The nrDoF = 6*nrNodes+(2*DEGREE-9)*nrEdges+
    %   (DEGREE-4)*(DEGREE-5)/2*nrElems global degrees of freedom list
    %   the nodal ones as in Argyris, then the edge and the interior
    %   ones.
    %
    %   fixedNodes holds the degrees of freedom of clamped conditions on
    %   the Dirichlet boundary as in Argyris, together with all edge
    %   degrees of freedom of the Dirichlet edges. This assumes
    %   Dirichlet edges parallel to the axes. The Neumann boundary fixes
    %   no degrees of freedom.
    %
    %   TUBAk functions:
    %       getBasis     - Evaluate the basis functions
    %       getGradBasis - Evaluate the gradients of the basis
    %       getD2Basis   - Evaluate the Hessians of the basis
    %       getD3Basis   - Evaluate third derivatives of the basis
    %       getD4Basis   - Evaluate fourth derivatives of the basis
    %
    %   TUBAk properties:
    %       mesh              - Mesh on which the element is defined
    %       nrDoF             - Total number of degrees of freedom
    %       dof4e             - Local-to-global dofs per element
    %       pk                - DGDk element used for the basis
    %       nrBasisFunc       - Number of basis functions per element
    %       degree            - Polynomial degree of the element
    %       fixedNodes        - Indices of the fixed degrees of freedom
    %       freeNodes         - Indices of the free degrees of freedom
    %       isRefElemBased    - Basis lives on the reference element
    %       basisCoefficients - Expansion coefficients in the DGDk basis
    %
    %   See also Argyris, DGDk, HHJk, FE

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = public)
        mesh              % mesh on which the element is defined
        nrDoF             % total number of degrees of freedom
        dof4e             % local-to-global dofs per element
        pk                % DGDk element used for the basis
        nrBasisFunc       % number of basis functions per element
        degree            % polynomial degree of the element
        fixedNodes        % indices of the fixed degrees of freedom
        freeNodes         % indices of the free degrees of freedom
        isRefElemBased = true % basis lives on the reference element
    end

    properties(Access = public)
        basisCoefficients % expansion coefficients in the DGDk basis
    end

    methods
        function obj = TUBAk(degree,mesh)
            %TUBAk - Create an Argyris-type finite element object
            if degree < 5
                error('fem4m:TUBAk:invalidDegree', ...
                    'TUBA degree must be >= 5, C1-conforming elements on triangles need degree 5');
            end
            obj.isAffineEquivalent = false;
            obj.degree = degree;
            obj.pk = DGDk(degree,mesh);
            obj.mesh = mesh;
            obj.nrBasisFunc = obj.pk.nrBasisFunc;
            obj = obj.getBasisCoefficients;
            obj = obj.dofManager;
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the TUBA
            %   basis functions at the reference point PTS_REF (1-by-2) in
            %   all elements by combining the DGDk basis with the stored
            %   expansion coefficients. VAL is nrElems-by-nrBasisFunc. PTS
            %   and PARTS are ignored.
            val = matMul(obj.basisCoefficients,obj.pk.getBasis([],x,[]));
        end

        function val = getGradBasis(obj,~,x,~)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   gradients of the TUBA basis functions at the reference
            %   point PTS_REF (1-by-2) in all elements. VAL is
            %   nrElems-by-nrBasisFunc-by-2. PTS and PARTS are ignored.
            val = matMul(obj.basisCoefficients,obj.pk.getGradBasis([],x,[]));
        end

        function val = getD2Basis(obj,~,x,~)
            %getD2Basis - Evaluate the Hessians of the basis
            %   VAL = getD2Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   second derivatives [xx yy xy] of the TUBA basis functions
            %   at the reference point PTS_REF (1-by-2) in all elements.
            %   VAL is nrElems-by-nrBasisFunc-by-3. PTS and PARTS are
            %   ignored.
            val = matMul(obj.basisCoefficients,obj.pk.getD2Basis([],x,[]));
        end

        function val = getD3Basis(obj,~,x,~)
            %getD3Basis - Evaluate third derivatives of the basis
            %   VAL = getD3Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   third derivatives [xxx yyy xxy xyy] of the TUBA basis
            %   functions at the reference point PTS_REF (1-by-2) in all
            %   elements. VAL is nrElems-by-nrBasisFunc-by-4. PTS and
            %   PARTS are ignored.
            val = matMul(obj.basisCoefficients,obj.pk.getD3Basis([],x,[]));
        end

        function val = getD4Basis(obj,~,x,~)
            %getD4Basis - Evaluate fourth derivatives of the basis
            %   VAL = getD4Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   fourth derivatives [xxxx yyyy xxyy xxxy xyyy] of the TUBA
            %   basis functions at the reference point PTS_REF (1-by-2) in
            %   all elements. VAL is nrElems-by-nrBasisFunc-by-5. PTS and
            %   PARTS are ignored.
            val = matMul(obj.basisCoefficients,obj.pk.getD4Basis([],x,[]));
        end

    end

    methods (Access = private)

        function obj = getBasisCoefficients(obj)
            %getBasisCoefficients - Compute the basis expansion coefficients

            nrElems = obj.mesh.nrElems;
            area4e = obj.mesh.area4e;
            % Initialisation for base functions
            obj.basisCoefficients = zeros(nrElems,obj.nrBasisFunc,obj.nrBasisFunc);
            % basis in nodes
            valueN1 = obj.pk.getBasis([],[0,0],[]);
            valueN2 = obj.pk.getBasis([],[1,0],[]);
            valueN3 = obj.pk.getBasis([],[0,1],[]);
            % gradBasis in nodes
            gradN1 = obj.pk.getGradBasis([],[0,0],[]);
            gradN2 = obj.pk.getGradBasis([],[1,0],[]);
            gradN3 = obj.pk.getGradBasis([],[0,1],[]);
            % D2Basis in nodes
            d2N1 = obj.pk.getD2Basis([],[0,0],[]);
            d2N2 = obj.pk.getD2Basis([],[1,0],[]);
            d2N3 = obj.pk.getD2Basis([],[0,1],[]);
            % normal directions
            normals4ed = obj.mesh.normal4ed;
            normals4e(:,1,:) = normals4ed(obj.mesh.ed4e(:,1),:)';
            normals4e(:,2,:) = normals4ed(obj.mesh.ed4e(:,2),:)';
            normals4e(:,3,:) = normals4ed(obj.mesh.ed4e(:,3),:)';
            valueE1 = obj.pk.getBasis([],[1/2,0],[])';
            valueE2 = obj.pk.getBasis([],[1/2,1/2],[])';
            valueE3 = obj.pk.getBasis([],[0,1/2],[])';
            if mod(obj.degree,2)==1
                N = obj.degree - 3;
            else
                N = obj.degree - 2;
            end
            gradE1 = zeros(obj.nrBasisFunc,2,obj.mesh.nrElems,N-1);
            gradE2 = zeros(obj.nrBasisFunc,2,obj.mesh.nrElems,N-1);
            gradE3 = zeros(obj.nrBasisFunc,2,obj.mesh.nrElems,N-1);
            for k = 1 : N-1
                gradE1(:,:,:,k) = permute(obj.pk.getGradBasis([],[k/N,0],[]),[2 3 1]);
                gradE2(:,:,:,k) = permute(obj.pk.getGradBasis([],[1-k/N,k/N],[]),[2 3 1]);
                gradE3(:,:,:,k) = permute(obj.pk.getGradBasis([],[0,1-k/N],[]),[2 3 1]);
            end
            indx = 1;
            valueT = zeros((obj.degree-4)*(obj.degree-5)/2,obj.nrBasisFunc,obj.mesh.nrElems);
            M = obj.degree - 3;
            for k = 1 : obj.degree - 4
                for j = 1 : obj.degree - 4 - k
                    valueT(indx,:,:) =  obj.pk.getBasis([],[k/M,j/M],[])';
                    indx = indx + 1;
                end
            end
            % Element specific transformation
            for curElem = 1 : nrElems
                C = zeros(obj.nrBasisFunc,obj.nrBasisFunc);
                H = sqrt(area4e(curElem));
                % Nodal degrees of freedom phi(n_j)
                C(1:3,:) = [valueN1(curElem,:);valueN2(curElem,:);valueN3(curElem,:)];
                % grad DoF
                C(4:6,:) = H.*[gradN1(curElem,:,1);gradN2(curElem,:,1);gradN3(curElem,:,1)];
                C(7:9,:) = H.*[gradN1(curElem,:,2);gradN2(curElem,:,2);gradN3(curElem,:,2)];
                % d2 DoF
                C(10:12,:) = H.^2.*[d2N1(curElem,:,1);d2N2(curElem,:,1);d2N3(curElem,:,1)];
                C(13:15,:) = H.^2.*[d2N1(curElem,:,2);d2N2(curElem,:,2);d2N3(curElem,:,2)];
                C(16:18,:) = H.^2.*[d2N1(curElem,:,3);d2N2(curElem,:,3);d2N3(curElem,:,3)];
                % right hand side
                I = ones(obj.nrBasisFunc,1);
                I(4:9) = H;
                I(10:18) = H.^2;
                indx = 18;
                for k = 1 : N-1
                    if k == N/2
                        if mod(obj.degree,2)==1
                            % Normal directions in mid(E)
                            C(indx+1,:) = H.*gradE1(:,:,curElem,N/2)*normals4e(:,1,curElem);
                            C(indx+2,:) = H.*gradE2(:,:,curElem,N/2)*normals4e(:,2,curElem);
                            C(indx+3,:) = H.*gradE3(:,:,curElem,N/2)*normals4e(:,3,curElem);
                            I(indx+(1:3)) = H;
                        else
                            % value in mid(E)
                            C(indx+1,:) = valueE1(:,curElem);
                            C(indx+2,:) = valueE2(:,curElem);
                            C(indx+3,:) = valueE3(:,curElem);
                            I(indx+(1:3)) = 1;
                        end
                        indx = indx + 3;
                    elseif k == N/2-1 && mod(obj.degree,2)==0
                        % Normal directions in mid(E)-1/N
                        C(indx+1,:) = H.*gradE1(:,:,curElem,N/2-1)*normals4e(:,1,curElem);
                        C(indx+2,:) = H.*gradE2(:,:,curElem,N/2-1)*normals4e(:,2,curElem);
                        C(indx+3,:) = H.*gradE3(:,:,curElem,N/2-1)*normals4e(:,3,curElem);
                        I(indx+(1:3)) = H;
                        indx = indx + 3;
                    elseif k == N/2+1 && mod(obj.degree,2)==0
                        % Normal directions in mid(E)+1/N
                        C(indx+1,:) = H.*gradE1(:,:,curElem,N/2+1)*normals4e(:,1,curElem);
                        C(indx+2,:) = H.*gradE2(:,:,curElem,N/2+1)*normals4e(:,2,curElem);
                        C(indx+3,:) = H.*gradE3(:,:,curElem,N/2+1)*normals4e(:,3,curElem);
                        I(indx+(1:3)) = H;
                        indx = indx + 3;
                    else
                        C(indx+1,:) = H.*gradE1(:,1,curElem,k);
                        C(indx+2,:) = H.*gradE2(:,1,curElem,k);
                        C(indx+3,:) = H.*gradE3(:,1,curElem,k);
                        C(indx+4,:) = H.*gradE1(:,2,curElem,k);
                        C(indx+5,:) = H.*gradE2(:,2,curElem,k);
                        C(indx+6,:) = H.*gradE3(:,2,curElem,k);
                        I(indx+(1:6)) = H;
                        indx = indx + 6;
                    end
                end
                % elem
                C(indx+1:end,:) = valueT(:,:,curElem);
                %right-hand side
                I = diag(I);
                % Assembly
                x = C\I;
                obj.basisCoefficients(curElem,:,:) = x';
            end

        end

    end

    methods (Access = protected)

        function obj = dofManager(obj)
            %dofManager - Assemble the degrees of freedom of the element
            nrNodes = obj.mesh.nrNodes;
            nrEdges = obj.mesh.nrEdges;
            nrElems = obj.mesh.nrElems;
            n4e = obj.mesh.n4e;
            ed4e = obj.mesh.ed4e;
            % degrees of freedom
            obj.dof4e = zeros(nrElems,obj.nrBasisFunc);
            % 3*6 node based
            obj.dof4e(:,1:3) = n4e;
            obj.dof4e(:,4:9) = [nrNodes + n4e,2*nrNodes + n4e];
            obj.dof4e(:,10:18) = [3*nrNodes + n4e,4*nrNodes + n4e,5*nrNodes + n4e];
            % 3*(2*degree - 9) edge based
            for k = 1 : 2*obj.degree-9
                obj.dof4e(:,18+(k-1)*3+(1:3)) = 6*nrNodes + (k-1)*nrEdges + ed4e;
            end
            % reverse the order of the edge dofs of every local edge
            % whose orientation is opposite to the global one, keeping
            % the x/y order of the gradient pairs
            for curEd = 1 : 3
                indx = obj.mesh.n4ed(ed4e(:,curEd),1) ~= obj.mesh.n4e(:,curEd);
                p = 18+((0:2*obj.degree-9-1)*3);
                q = p(end:-1:1);
                if obj.degree >= 7
                    if mod(obj.degree,2)==1 % odd
                        indy = setdiff(1:numel(p),(numel(p)+1)/2);
                    else %even
                        indy = setdiff(1:numel(p),[(numel(p)+1)/2,(numel(p)+1)/2-1,(numel(p)+1)/2+1]);
                    end
                    for j = indy(1:2:end)
                        q([j j+1]) = q([j+1 j]);
                    end
                end
                obj.dof4e(indx,p+curEd) = obj.dof4e(indx,q+curEd);
            end
            % elem based
            for k = 1:(obj.degree-4)*(obj.degree-5)/2
                obj.dof4e(:,18+3*(2*obj.degree-9)+k) = 6*nrNodes + (2*obj.degree-9)*nrEdges + (k-1)*nrElems + (1:nrElems)';
            end
            % nrDoF
            obj.nrDoF = 6*nrNodes + (2*obj.degree-9)*nrEdges + (obj.degree-4)*(obj.degree-5)/2*nrElems;
            % fix & free
            obj.fixedNodes = unique([obj.mesh.n4Db(:);...
                nrNodes + obj.mesh.n4Db(:);2*nrNodes + obj.mesh.n4Db(:);...
                5*nrNodes + obj.mesh.n4Db(:)]);
            for k = 1 : 2*obj.degree-9
                obj.fixedNodes = [obj.fixedNodes; 6*nrNodes + (k-1)*nrEdges + obj.mesh.ed4Db(:)];
            end
            % xx and yy derivatives only for Dirichlet edges parallel to
            % the axes (!!!)
            normals = obj.mesh.normal4ed;
            normals = normals(obj.mesh.ed4Db,:);
            xEdges = obj.mesh.n4Db(normals(:,1)==0,:);
            yEdges = obj.mesh.n4Db(normals(:,2)==0,:);
            obj.fixedNodes = unique([obj.fixedNodes(:);3*nrNodes + xEdges(:)] );
            obj.fixedNodes = unique([obj.fixedNodes(:);4*nrNodes + yEdges(:)] );
            obj.freeNodes = setdiff(1:obj.nrDoF, obj.fixedNodes);
        end

    end

end