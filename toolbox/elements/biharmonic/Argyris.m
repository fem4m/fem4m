classdef Argyris < FE
    %ARGYRIS - Argyris quintic finite element for biharmonic problems
    %   OBJ = ARGYRIS(DEGREE,MESH) creates the Argyris P5 H2-conforming
    %   finite element on the triangular mesh MESH with 21 degrees of
    %   freedom per element: values, gradients and Hessians at the nodes
    %   and normal derivatives at the edge midpoints. DEGREE must be 5;
    %   TUBAk provides all degrees DEGREE >= 5. The basis coefficients
    %   depend on the element geometry (isAffineEquivalent is false).
    %
    %   The nrDoF = 6*nrNodes+nrEdges global degrees of freedom are the
    %   values, the x- and y-derivatives and the second derivatives
    %   [xx yy xy] at all nodes, one block of nrNodes each, followed by
    %   the normal derivatives along normal4ed at the edge midpoints.
    %
    %   fixedNodes holds the degrees of freedom of clamped conditions
    %   u = grad(u)*n = 0 on the Dirichlet boundary: the values, the
    %   gradients and the xy derivatives at the Dirichlet nodes, the
    %   normal derivatives on the Dirichlet edges, and the xx (yy)
    %   derivatives at the nodes of Dirichlet edges whose normal has a
    %   zero x (y) component. This assumes Dirichlet edges parallel to
    %   the axes. The Neumann boundary fixes no degrees of freedom.
    %
    %   ARGYRIS functions:
    %       getBasis     - Evaluate the basis functions
    %       getGradBasis - Evaluate the gradients of the basis
    %       getD2Basis   - Evaluate the Hessians of the basis
    %       enrich       - Average a function on the same mesh into the space
    %
    %   ARGYRIS properties:
    %       mesh              - Mesh on which the element is defined
    %       nrDoF             - Total number of degrees of freedom
    %       dof4e             - Local-to-global dofs per element
    %       p5                - Quintic Pk element used for the basis
    %       nrBasisFunc       - Number of basis functions per element
    %       degree            - Polynomial degree of the element
    %       fixedNodes        - Indices of the fixed degrees of freedom
    %       freeNodes         - Indices of the free degrees of freedom
    %       isRefElemBased    - Basis lives on the reference element
    %       basisCoefficients - Expansion coefficients in the Pk basis
    %
    %   See also Pk, HHJk, TUBAk, FE

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = public)
        mesh              % mesh on which the element is defined
        nrDoF             % total number of degrees of freedom
        dof4e             % local-to-global dofs per element
        p5                % quintic Pk element used for the basis
        nrBasisFunc       % number of basis functions per element
        degree            % polynomial degree of the element
        fixedNodes        % indices of the fixed degrees of freedom
        freeNodes         % indices of the free degrees of freedom
        isRefElemBased = true % basis lives on the reference element
    end

    properties(Access = public)
        basisCoefficients % expansion coefficients in the Pk basis
    end

    methods
        function obj = Argyris(degree,mesh)
            %ARGYRIS - Create an Argyris finite element object
            if nargin < 2 || ~isequal(degree,5)
                error('fem4m:Argyris:invalidDegree', ...
                    'Argyris degree must be 5, TUBAk has the other degrees');
            end
            obj.isAffineEquivalent = false;
            obj.degree = degree;
            obj.p5 = Pk(degree,mesh);
            obj.mesh = mesh;
            obj.nrBasisFunc = obj.p5.nrBasisFunc;
            obj = getBasisCoefficients(obj);
            obj = obj.dofManager;
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the Argyris
            %   basis functions at the reference point PTS_REF (1-by-2) in
            %   all elements by combining the quintic Pk basis with the
            %   stored expansion coefficients. VAL is nrElems-by-21. PTS
            %   and PARTS are ignored.
            val = matMul(obj.basisCoefficients,obj.p5.getBasis([],x,[]));
        end

        function val = getGradBasis(obj,~,x,~)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   gradients of the Argyris basis functions at the reference
            %   point PTS_REF (1-by-2) in all elements. VAL is
            %   nrElems-by-21-by-2. PTS and PARTS are ignored.
            val = matMul(obj.basisCoefficients,obj.p5.getGradBasis([],x,[]));
        end

        function val = getD2Basis(obj,~,x,~)
            %getD2Basis - Evaluate the Hessians of the basis
            %   VAL = getD2Basis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   second derivatives [xx yy xy] of the Argyris basis
            %   functions at the reference point PTS_REF (1-by-2) in all
            %   elements. VAL is nrElems-by-21-by-3. PTS and PARTS are
            %   ignored.
            val = matMul(obj.basisCoefficients,obj.p5.getD2Basis([],x,[]));
        end

        function val = enrich(obj,c,fe)
            %ENRICH - Average a function on the same mesh into the space
            %   VAL = ENRICH(OBJ,C,FE) returns the coefficients VAL
            %   (nrDoF-by-1) of the Argyris function whose degrees of
            %   freedom are the values, gradients and second derivatives
            %   at the nodes and the normal derivatives at the edge
            %   midpoints of the function with coefficients C in the
            %   finite element FE, evaluated in each element and averaged
            %   over the elements that share a degree of freedom (an
            %   enriching operator). FE must be defined on the mesh of
            %   OBJ.
            n4e = obj.mesh.n4e;
            c4n = obj.mesh.c4n;
            ed4e = obj.mesh.ed4e;
            mid = obj.mesh.point4ed(1/2);
            normals4ed = obj.mesh.normal4ed;
            normals4e(:,:,1) = normals4ed(ed4e(:,1),:);
            normals4e(:,:,2) = normals4ed(ed4e(:,2),:);
            normals4e(:,:,3) = normals4ed(ed4e(:,3),:);

            val4dof = zeros(size(obj.dof4e));
            refC = [0 0; 1 0; 0 1];
            % dofs in vertices
            for i = 1 : 3
                nodes = n4e(:,i);
                pts = c4n(nodes,:);
                val4dof(:,i) = fe.getValue(pts,refC(i,:),1:obj.mesh.nrElems,c);
                val4dof(:,[3+i,6+i]) = fe.getGradValue(pts,refC(i,:),1:obj.mesh.nrElems,c);
                val4dof(:,[9+i,12+i,15+i]) = fe.getD2Value(pts,refC(i,:),1:obj.mesh.nrElems,c);
            end
            % dofs on edges
            refC = [1/2 0; 1/2 1/2; 0 1/2];
            for i = 1 : 3
                edges = ed4e(:,i);
                pts = mid(edges,:);
                grad = fe.getGradValue(pts,refC(i,:),1:obj.mesh.nrElems,c);
                val4dof(:,18+i) = sum(grad.*normals4e(:,:,i),2);
            end

            N = accumarray(obj.dof4e(:),ones(numel(obj.dof4e),1),[obj.nrDoF 1]);
            val = accumarray(obj.dof4e(:),val4dof(:),[obj.nrDoF 1]);
            val = val./N;
        end


    end

    methods (Access = private)

        function obj = getBasisCoefficients(obj)
            %getBasisCoefficients - Compute the basis expansion coefficients

            nrElems = obj.mesh.nrElems;
            area4e = obj.mesh.area4e;
            % Initialisation for base functions
            obj.basisCoefficients = zeros(nrElems,21,21);
            % gradBasis in nodes
            gradN1 = obj.p5.getGradBasis([],[0,0],[]);
            gradN2 = obj.p5.getGradBasis([],[1,0],[]);
            gradN3 = obj.p5.getGradBasis([],[0,1],[]);
            % D2Basis in nodes
            d2N1 = obj.p5.getD2Basis([],[0,0],[]);
            d2N2 = obj.p5.getD2Basis([],[1,0],[]);
            d2N3 = obj.p5.getD2Basis([],[0,1],[]);
            % normal directions
            normals4ed = obj.mesh.normal4ed;
            normals4e(:,1,:) = normals4ed(obj.mesh.ed4e(:,1),:)';
            normals4e(:,2,:) = normals4ed(obj.mesh.ed4e(:,2),:)';
            normals4e(:,3,:) = normals4ed(obj.mesh.ed4e(:,3),:)';
            gradE1 = permute(obj.p5.getGradBasis([],[1/2,0],[]),[2 3 1]);
            gradE2 = permute(obj.p5.getGradBasis([],[1/2,1/2],[]),[2 3 1]);
            gradE3 = permute(obj.p5.getGradBasis([],[0,1/2],[]),[2 3 1]);
            % Element specific transformation
            for curElem = 1 : nrElems
                C = zeros(21,21);
                H = sqrt(area4e(curElem));
                % Nodal degrees of freedom phi(n_j)
                C(1:3,1:3) = eye(3);
                % grad DoF
                C(4:6,:) = H.*[gradN1(curElem,:,1);gradN2(curElem,:,1);gradN3(curElem,:,1)];
                C(7:9,:) = H.*[gradN1(curElem,:,2);gradN2(curElem,:,2);gradN3(curElem,:,2)];
                % d2 DoF
                C(10:12,:) = H.^2.*[d2N1(curElem,:,1);d2N2(curElem,:,1);d2N3(curElem,:,1)];
                C(13:15,:) = H.^2.*[d2N1(curElem,:,2);d2N2(curElem,:,2);d2N3(curElem,:,2)];
                C(16:18,:) = H.^2.*[d2N1(curElem,:,3);d2N2(curElem,:,3);d2N3(curElem,:,3)];
                % Normal directions
                C(19,:) = H.*gradE1(:,:,curElem)*normals4e(:,1,curElem);
                C(20,:) = H.*gradE2(:,:,curElem)*normals4e(:,2,curElem);
                C(21,:) = H.*gradE3(:,:,curElem)*normals4e(:,3,curElem);
                %right-hand side
                I = ones(21,1);
                I(4:9) = H;
                I(10:18) = H.^2;
                I(19:21) = H;
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
            obj.dof4e = zeros(nrElems,21);
            obj.dof4e(:,1:3) = n4e;
            obj.dof4e(:,4:9) = [nrNodes + n4e,2*nrNodes + n4e];
            obj.dof4e(:,10:18) = [3*nrNodes + n4e,4*nrNodes + n4e,5*nrNodes + n4e];
            obj.dof4e(:,19:21) = 6*nrNodes + ed4e;
            % nrDoF
            obj.nrDoF = 6*nrNodes + nrEdges;
            % fix & free
            obj.fixedNodes = unique([obj.mesh.n4Db(:);...
                nrNodes + obj.mesh.n4Db(:);2*nrNodes + obj.mesh.n4Db(:);...
                5*nrNodes + obj.mesh.n4Db(:); 6*nrNodes+obj.mesh.ed4Db(:)]);
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