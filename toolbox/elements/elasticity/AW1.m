classdef AW1 < FE
    %AW1 - Arnold-Winther element of polynomial degree three
    %   OBJ = AW1(DEGREE,MESH) creates the cubic Arnold-Winther
    %   symmetric H(div)-conforming stress finite element on the
    %   triangular mesh MESH with 24 degrees of freedom per element.
    %   DEGREE must be 3; the 1 in the name is the degree of the
    %   displacements. The space consists of the symmetric cubic matrix
    %   fields with linear divergence; its degrees of freedom are the
    %   values of sigma at the vertices (3 per vertex), the moments of
    %   degree 0 and 1 of sigma*n on every edge (4 per edge), and the
    %   element means of the three components of sigma (3 per element).
    %   The basis coefficients depend on the element geometry
    %   (isAffineEquivalent is false). AWH hybridizes the vertex
    %   constraints at non-parallel Neumann corners.
    %
    %   AW1 functions:
    %       prepareCoefficients - Compute the basis expansion coefficients
    %       getBasis            - Evaluate the Arnold-Winther basis functions
    %       getDivBasis         - Evaluate the divergence of the basis
    %       basisP3             - Evaluate the cubic reference basis
    %       gradP3              - Evaluate the gradients of the cubic basis
    %       basisP1             - Evaluate the linear reference basis
    %       gradP1              - Evaluate the gradients of the linear basis
    %       boundaryConstraints - Build the vertex boundary constraints
    %       boundaryMoments     - Assemble the Neumann boundary moments
    %       interpolateBoundary - Interpolate Neumann data into dofs
    %       intMoment0          - Evaluate the zeroth order moment integrand
    %       intMoment1          - Evaluate the first order moment integrand
    %
    %   AW1 properties:
    %       mesh              - Mesh on which the element is defined
    %       degree            - Polynomial degree of the element
    %       nrBasisFunc       - Number of basis functions per element
    %       dof4e             - Local-to-global dofs per element
    %       nrDoF             - Total number of degrees of freedom
    %       basisCoefficients - Expansion coefficients in the P3 basis
    %       hybrid            - True for hybridized vertex constraints (AWH)
    %       fixedNodes        - Indices of the fixed degrees of freedom
    %       freeNodes         - Indices of the free degrees of freedom
    %       Constraints       - Constraint matrix of the boundary vertices
    %       scaling           - Scaling factors of the degrees of freedom
    %       dof4ed            - Local-to-global dofs per edge
    %       nrBasis4n         - Number of nodal basis functions
    %       nrBasis4ed        - Number of edge basis functions
    %       nrBasis4e         - Number of interior basis functions
    %       hybridNodes       - Number of hybridized corners, 3 dofs each
    %       isRefElemBased    - Basis lives on the reference element
    %       momentDegree      - Quadrature degree of the Neumann moments
    %
    %   See also AWk, AWH, AW0, HZ1, BDMk, FEMatrix

    % Copyright (c) 2026 by the fem4m authors

    properties
        mesh              % mesh on which the element is defined
        degree            % polynomial degree of the element
        nrBasisFunc       % number of basis functions per element
        dof4e             % local-to-global dofs per element
        nrDoF             % total number of degrees of freedom
        basisCoefficients % expansion coefficients in the P3 basis
        fixedNodes        % indices of the fixed degrees of freedom
        freeNodes         % indices of the free degrees of freedom
        Constraints       % constraint matrix of the boundary vertices
        scaling           % scaling factors of the degrees of freedom
        dof4ed            % local-to-global dofs per edge
        nrBasis4n         % number of nodal basis functions
        nrBasis4ed        % number of edge basis functions
        nrBasis4e         % number of interior basis functions
        hybridNodes       % number of hybridized corners, 3 dofs each
        isRefElemBased = true % basis lives on the reference element
    end

    properties (SetAccess = protected)
        hybrid = false    % true for hybridized vertex constraints (AWH)
    end

    properties (Constant)
        momentDegree = 12 % quadrature degree of the Neumann moments
    end

    methods

        function obj = AW1(degree,mesh)
            %AW1 - Create an Arnold-Winther finite element object
            % also for AWk, whose constructor calls AW1() without arguments
            obj.isAffineEquivalent = false;
            obj.degree = 3;
            if nargin > 0
                if nargin < 2 || ~isequal(degree,3)
                    error(['fem4m:',class(obj),':invalidDegree'], ...
                        '%s degree must be 3, the stresses are cubic',class(obj));
                end
                obj.mesh = mesh;
                obj = obj.prepareCoefficients();
                obj = obj.dofCounts();
                obj = obj.dofManager();
                obj = obj.boundaryConstraints();
            end
        end

        %% Calculation of the Arnold Winther transformation
        function obj = prepareCoefficients(obj)
            %prepareCoefficients - Compute the basis expansion coefficients
            %   OBJ = prepareCoefficients(OBJ) solves on every element the
            %   local 30x30 Vandermonde system given by the nodal values,
            %   the zeroth and first order edge moments of sigma*n, the
            %   element means of sigma, and the six constraints
            %   div(sigma) in P1 to obtain the expansion of the
            %   Arnold-Winther basis in the P3 basis (basisCoefficients).
            normals4e = obj.mesh.normal4e;
            nrElems = obj.mesh.nrElems;
            gradP14e = obj.gradP1;
            gradP14e = permute(matMul(obj.mesh.transformation,permute(gradP14e,[3 2 1])),[1 3 2]);
            ed4e = obj.mesh.ed4e;
            normal4ed = obj.mesh.normal4ed;
            length4ed = obj.mesh.length4ed;
            % Initialisation
            obj.basisCoefficients = zeros(nrElems,24,30);
            C = zeros(30,30);
            % Nodal degrees of freedom
            C(1:9,1:9) = 60*eye(9,9);
            % Element degrees of freedom
            I = [eye(3,3),eye(3,3),eye(3,3)];
            C(22:24,:) = [20*I,5*I,0*I,eye(3,3)];
            % Element specific transformation
            for curElem = 1 : nrElems
                % Edge degrees of freedom
                curNormal = normal4ed(ed4e(curElem,:),:);
                N1 = [diag(curNormal(1,:)),curNormal(1,[2 1])'];
                N2 = [diag(curNormal(2,:)),curNormal(2,[2 1])'];
                N3 = [diag(curNormal(3,:)),curNormal(3,[2 1])'];
                curNormal = permute(normals4e(curElem,:,:),[2 3 1]);
                M1 = [diag(curNormal(1,:)),curNormal(1,[2 1])'];
                M2 = [diag(curNormal(2,:)),curNormal(2,[2 1])'];
                M3 = [diag(curNormal(3,:)),curNormal(3,[2 1])'];
                O = zeros(2,3);
                R = [30*N1,30*N1,O;-5*M1,5*M1,O;
                    O,30*N2,30*N2;O,-5*M2,5*M2;
                    30*N3,O,30*N3;5*M3,O,-5*M3];
                S1 = [10*N1,O,O;O,O,O;
                    O,10*N2,O;O,O,O;
                    O,O,10*N3;O,O,O];
                S2 = [O,O,O;-M1,O,O;
                    O,O,O;O,-M2,O;
                    O,O,O;O,O,-M3];
                C(10:21,1:27) = [R,S1,S2];
                % Constraint div(sigma) in P1
                gradP1X = gradP14e(curElem,:,1);
                gradP1Y = gradP14e(curElem,:,2);
                dxW = [3*(gradP1X(1)-gradP1X(2)), -gradP1X(3), gradP1X(3), gradP1X(3);
                    gradP1X(1), 3*(gradP1X(2)-gradP1X(3)), -gradP1X(1), gradP1X(1);
                    -gradP1X(2), gradP1X(2), 3*(gradP1X(3)-gradP1X(1)), gradP1X(2)];
                dyW = [3*(gradP1Y(1)-gradP1Y(2)), -gradP1Y(3), gradP1Y(3), gradP1Y(3);
                    gradP1Y(1), 3*(gradP1Y(2)-gradP1Y(3)), -gradP1Y(1), gradP1Y(1);
                    -gradP1Y(2), gradP1Y(2), 3*(gradP1Y(3)-gradP1Y(1)), gradP1Y(2)];
                O = zeros(3,1);
                Q = [dxW(:,1),O,dyW(:,1),dxW(:,2),O,dyW(:,2),dxW(:,3),O,dyW(:,3),dxW(:,4),O,dyW(:,4);
                    O,dyW(:,1),dxW(:,1),O,dyW(:,2),dxW(:,2),O,dyW(:,3),dxW(:,3),O,dyW(:,4),dxW(:,4)];
                C(25:30,19:30) = 60*max(length4ed(ed4e(curElem,:)))*Q;
                % Assembly
                x = (C/60)\[eye(24);zeros(6,24)];
                x(abs(x)<1e-12) = 0;
                obj.basisCoefficients(curElem,:,:) = x';
            end
        end

        %% Stress basis functions
        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the Arnold-Winther basis functions
            %   VAL = getBasis(OBJ,PTS,X,PARTS) evaluates the symmetric
            %   cubic Arnold-Winther stress basis functions at the point
            %   X (1-by-2) of the reference triangle in all elements. VAL
            %   is nrElems-by-nrBasisFunc-by-3 with the components
            %   [s11 s22 s12] (Voigt notation), each basis function
            %   divided by the scaling of its degree of freedom. PTS and
            %   PARTS are ignored.
            basisPk = obj.basisP3(x(1),x(2))';
            a = obj.basisCoefficients(:,:,1:3:end);
            b = obj.basisCoefficients(:,:,2:3:end);
            c = obj.basisCoefficients(:,:,3:3:end);
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,3);
            val(:,:,1) = matMul(a,basisPk);
            val(:,:,2) = matMul(b,basisPk);
            val(:,:,3) = matMul(c,basisPk);
            val = val./obj.scaling(obj.dof4e);
        end

        function val = getDivBasis(obj,~,x,~)
            %getDivBasis - Evaluate the divergence of the basis
            %   VAL = getDivBasis(OBJ,PTS,X,PARTS) evaluates the row-wise
            %   divergences of the Arnold-Winther stress basis functions
            %   at the point X (1-by-2) of the reference triangle in all
            %   elements. VAL is nrElems-by-nrBasisFunc-by-2, scaled as
            %   in getBasis. PTS and PARTS are ignored.
            gradBasisPk = obj.gradP3(x(1),x(2));
            gradBasisPk = permute(matMul(obj.mesh.transformation,permute(gradBasisPk,[3 2 1])),[1 3 2]);
            a = obj.basisCoefficients(:,:,1:3:end);
            b = obj.basisCoefficients(:,:,2:3:end);
            c = obj.basisCoefficients(:,:,3:3:end);
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,1);
            val(:,:,1) = matMul(a,gradBasisPk(:,:,1)) + matMul(c,gradBasisPk(:,:,2));
            val(:,:,2) = matMul(c,gradBasisPk(:,:,1)) + matMul(b,gradBasisPk(:,:,2));
            val = val./obj.scaling(obj.dof4e);
        end

        %% P3 basis functions
        function val = basisP3(obj,x,y)
            %basisP3 - Evaluate the cubic reference basis
            %   VAL = basisP3(OBJ,X,Y) evaluates the ten hierarchical
            %   cubic basis functions lambda_i, lambda_i*lambda_j and
            %   lambda_i*lambda_j*(lambda_i-lambda_j) for (i,j) = (1,2),
            %   (2,3), (3,1), and lambda_1*lambda_2*lambda_3 of the
            %   reference triangle (lambda of basisP1) at the point
            %   (X,Y). VAL is 10-by-1.
            basisP1 = obj.basisP1(x,y);
            val = zeros(10,1);
            val( 1) = basisP1(1);
            val( 2) = basisP1(2);
            val( 3) = basisP1(3);
            val( 4) = basisP1(1).*basisP1(2);
            val( 5) = basisP1(2).*basisP1(3);
            val( 6) = basisP1(3).*basisP1(1);
            val( 7) = basisP1(1).*basisP1(2).*(basisP1(1) - basisP1(2));
            val( 8) = basisP1(2).*basisP1(3).*(basisP1(2) - basisP1(3));
            val( 9) = basisP1(3).*basisP1(1).*(basisP1(3) - basisP1(1));
            val(10) = basisP1(1).*basisP1(2).*basisP1(3);
        end

        function val = gradP3(obj,x,y)
            %gradP3 - Evaluate the gradients of the cubic basis
            %   VAL = gradP3(OBJ,X,Y) evaluates the gradients of the ten
            %   cubic basis functions of basisP3 with respect to the
            %   reference coordinates at the point (X,Y). VAL is 10-by-2,
            %   one row per basis function.
            gradP1 = obj.gradP1;
            basisP1 = obj.basisP1(x,y);

            val = zeros(10,2);
            val(1,:) = gradP1(1,:);
            val(2,:) = gradP1(2,:);
            val(3,:) = gradP1(3,:);
            val(4,:) = gradP1(1,:)*basisP1(2) + gradP1(2,:)*basisP1(1);
            val(5,:) = gradP1(2,:)*basisP1(3) + gradP1(3,:)*basisP1(2);
            val(6,:) = gradP1(3,:)*basisP1(1) + gradP1(1,:)*basisP1(3);
            val(7,:) = gradP1(1,:)*basisP1(2)*(basisP1(1) - basisP1(2)) + ...
                basisP1(1)*gradP1(2,:)*(basisP1(1) - basisP1(2)) + ...
                basisP1(1)*basisP1(2)*(gradP1(1,:) - gradP1(2,:));
            val(8,:) = gradP1(2,:)*basisP1(3)*(basisP1(2) - basisP1(3)) + ...
                basisP1(2)*gradP1(3,:)*(basisP1(2) - basisP1(3)) + ...
                basisP1(2)*basisP1(3)*(gradP1(2,:) - gradP1(3,:));
            val(9,:) = gradP1(3,:)*basisP1(1)*(basisP1(3) - basisP1(1)) + ...
                basisP1(3)*gradP1(1,:)*(basisP1(3) - basisP1(1)) + ...
                basisP1(3)*basisP1(1)*(gradP1(3,:) - gradP1(1,:));
            val(10,:) = gradP1(1,:)*basisP1(2)*basisP1(3) + ...
                basisP1(1)*gradP1(2,:)*basisP1(3) + ...
                basisP1(1)*basisP1(2)*gradP1(3,:);
        end

        %% P1 basis functions
        function val = basisP1(~,x,y)
            %basisP1 - Evaluate the linear reference basis
            %   VAL = basisP1(OBJ,X,Y) evaluates the three linear nodal
            %   basis functions 1-X-Y, X, and Y of the reference triangle
            %   at the point (X,Y). VAL is 3-by-1.
            val = zeros(3,1);
            val(1) = 1 - x - y;
            val(2) = x;
            val(3) = y;
        end

        function val = gradP1(~)
            %gradP1 - Evaluate the gradients of the linear basis
            %   VAL = gradP1(OBJ) returns the constant gradients of the
            %   three linear nodal basis functions of the reference
            %   triangle. VAL is 3-by-2, one row per basis function.
            val = zeros(3,2);
            val(1,:) = [-1,-1];
            val(2,:) = [1,0];
            val(3,:) = [0,1];
        end

        function obj = boundaryConstraints(obj)
            %boundaryConstraints - Build the vertex boundary constraints
            %   OBJ = boundaryConstraints(OBJ) assembles the sparse
            %   constraint matrix enforcing the Neumann and symmetry
            %   conditions at boundary vertices, fixing nodal dofs or
            %   adding hybrid dofs at non-parallel Neumann corners. It
            %   sets Constraints, fixedNodes, freeNodes, and hybridNodes;
            %   with hybrid dofs (AWH) it also extends nrDoF, scaling, and
            %   dof4e.
            % bd conditions at the vertices
            % vertices are constrained via multiplier
            C = [];

            nrNodes = obj.mesh.nrNodes;
            n4Nb = obj.mesh.n4Nb;
            n4Sb = obj.mesh.n4Sb;
            obj.hybridNodes = 0;

            if ~isempty(n4Nb)
                n4Nb = obj.mesh.n4Nb;

                NbEd4n = obj.mesh.NbEd4n;
                normal4ed = obj.mesh.normal4ed;
                e4n = obj.mesh.e4n;
                e4ed = obj.mesh.e4ed;

                NbNodes = unique(n4Nb(:));
                NbNodesDof = zeros(nrNodes,2);
                NbNodesDof(NbNodes,:) = reshape(1:2*length(NbNodes),[],2);
                fixNbNodes = [];
                % nodes
                C = sparse(2*length(NbNodes),obj.nrDoF);
                for j = 1 : length(NbNodes)
                    curNbEdges = NbEd4n(NbNodes(j),:);
                    % one Nb one Db
                    if nnz(curNbEdges) < 2
                        if curNbEdges(2) == 0
                            curNormal = normal4ed(curNbEdges(1),:);
                        else
                            curNormal = normal4ed(curNbEdges(2),:);
                        end
                        C(NbNodesDof(NbNodes(j),1),NbNodes(j)) = curNormal(1); %#ok<SPRIX>
                        C(NbNodesDof(NbNodes(j),2),nrNodes+NbNodes(j)) = curNormal(2); %#ok<SPRIX>
                        C(NbNodesDof(NbNodes(j),:),2*nrNodes+NbNodes(j)) = curNormal([2 1]); %#ok<SPRIX>
                    else % two adjacent Nb Edges
                        curNormals = normal4ed(curNbEdges,:);
                        % Parallel
                        if abs(abs(sum(curNormals(1,:).*curNormals(2,:)))-1) <= eps
                            curNormal = curNormals(1,:);
                            C(NbNodesDof(NbNodes(j),1),NbNodes(j)) = curNormal(1); %#ok<SPRIX>
                            C(NbNodesDof(NbNodes(j),2),nrNodes+NbNodes(j)) = curNormal(2); %#ok<SPRIX>
                            C(NbNodesDof(NbNodes(j),:),2*nrNodes+NbNodes(j)) = curNormal([2 1]); %#ok<SPRIX>
                        else % non parallel
                            if ~obj.hybrid
                                fixNbNodes = [fixNbNodes,NbNodesDof(NbNodes(j),:)]; %#ok<AGROW>
                                dof = [NbNodes(j),nrNodes+NbNodes(j),2*nrNodes+NbNodes(j)];
                                obj.fixedNodes = [obj.fixedNodes;dof(:)];
                            else
                                fixNbNodes = [fixNbNodes,NbNodesDof(NbNodes(j),:)]; %#ok<AGROW>
                                dof = [NbNodes(j),nrNodes+NbNodes(j),2*nrNodes+NbNodes(j)];
                                newDof = obj.nrDoF+(1:3);
                                obj.hybridNodes = obj.hybridNodes + 1;
                                obj.scaling(newDof) = obj.scaling(dof);
                                obj.nrDoF = obj.nrDoF + 3;
                                C(:,end:end+3) = 0; %#ok<SPRIX>
                                if numel(find(e4n(NbNodes(j),:))) ~= 2; error('fem4m:AW1:invalidCorner','corner without two triangles');end
                                T = e4ed(curNbEdges,1);
                                indx = find(obj.dof4e(T(2),:)==dof(1)|obj.dof4e(T(2),:)==dof(2)...
                                    |obj.dof4e(T(2),:)==dof(3));
                                obj.dof4e(T(2),indx) = newDof;
                                obj.fixedNodes = [obj.fixedNodes;dof(:);newDof(:)];
                            end
                        end
                    end
                end
                freeNbNodes = setdiff(1:2*length(NbNodes),fixNbNodes);
                C = C(freeNbNodes,:);
            end

            % Symmetric bd conditions - only for separate axis-parallel lines
            % normal is assumed to be constant along each line of boundary part
            % symmetry condition is constant along each separate boundary part
            if ~isempty(n4Sb)
                % node values
                SbNodes = n4Sb(:,1:2);
                [SbNodes,indx] = unique(SbNodes(:));
                SbType4Nodes = [n4Sb(:,3);n4Sb(:,3)];
                SbType4Nodes = SbType4Nodes(indx);
                BSb = sparse(length(SbNodes),obj.nrDoF);
                SbEd4n = obj.mesh.SbEd4n;
                normal4ed = obj.mesh.normal4ed;
                fixNodesSb = [];
                for j = 1 : length(SbNodes)
                    curSbEdges = SbEd4n(SbNodes(j),:);
                    curNbEdges = obj.mesh.NbEd4n(SbNodes(j),:);
                    if nnz(curNbEdges)>0
                        %case treated by Neumann bdy cond.
                        fixNodesSb = [fixNodesSb,j]; %#ok<AGROW>
                    else
                        if curSbEdges(2) == 0
                            curNormal = normal4ed(curSbEdges(1),:);
                        else
                            curNormal = normal4ed(curSbEdges(2),:);
                        end
                        curType = SbType4Nodes(j);
                        if curType == 2
                            BSb(j,SbNodes(j)) = curNormal(1); %#ok<SPRIX> %alpha
                            BSb(j,2*nrNodes+SbNodes(j)) = curNormal(2); %#ok<SPRIX> %gamma
                        else
                            BSb(j,2*nrNodes+SbNodes(j)) = curNormal(1); %#ok<SPRIX> % gamma
                            BSb(j,nrNodes+SbNodes(j)) = curNormal(2); %#ok<SPRIX> % beta
                        end
                    end
                end
                BSb(fixNodesSb,:)=[];
                C = [C;BSb];
            end

            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
            obj.Constraints = C;
        end

        function x = boundaryMoments(obj,g)
            %boundaryMoments - Assemble the Neumann boundary moments
            %   X = boundaryMoments(OBJ,G) computes the zeroth and first
            %   order edge moments of the traction function handle
            %   G(PTS,NORMALS) on all Neumann boundary edges and returns
            %   them at the Neumann edge dofs of the vector X
            %   (nrDoF-by-1, zero elsewhere).
            x = zeros(obj.nrDoF,1);
            if ~isempty(obj.mesh.n4Nb)
                % moments
                ed4Nb = obj.mesh.ed4Nb;
                length4ed = obj.mesh.length4ed;
                % order 0
                moments0 = (1./length4ed(obj.mesh.ed4Nb)*ones(1,2)).*integrate(...
                    @(points,refp,ind)obj.intMoment0(points,refp,ind,g,obj.mesh),obj.mesh.n4Nb,obj.mesh.ed4Nb,obj.mesh, ...
                    obj.momentDegree);
                fixNodes0 = 3*obj.mesh.nrNodes+[ed4Nb,obj.mesh.nrEdges+ed4Nb];
                x(fixNodes0(:)) = moments0(:);
                % order 1
                moments1 = (1./length4ed(obj.mesh.ed4Nb)).^2*ones(1,2).*integrate(...
                    @(points,refp,ind)obj.intMoment1(points,refp,ind,g,obj.mesh),obj.mesh.n4Nb,obj.mesh.ed4Nb,obj.mesh, ...
                    obj.momentDegree);
                fixNodes1 = 3*obj.mesh.nrNodes+[2*obj.mesh.nrEdges+ed4Nb,3*obj.mesh.nrEdges+ed4Nb];
                x(fixNodes1(:)) = moments1(:);
            end
        end

        function [x,bN] = interpolateBoundary(obj,g)
            %interpolateBoundary - Interpolate Neumann data into dofs
            %   X = interpolateBoundary(OBJ,G) interpolates the traction
            %   function handle G(PTS,NORMALS) into the boundary degrees
            %   of freedom (the Neumann edge moments of boundaryMoments
            %   and the nodal values at non-parallel Neumann corners),
            %   solving the local corner systems. X (nrDoF-by-1) is
            %   multiplied by scaling, so it holds the coefficients of
            %   the scaled basis of getBasis.
            %
            %   [X,bN] = interpolateBoundary(OBJ,G) also returns the
            %   right-hand side bN of the vertex constraint system
            %   Constraints*X = bN.
            x = obj.boundaryMoments(g);
            bN = [];
            % node constraint right-hand side
            nrNodes = obj.mesh.nrNodes;
            n4Nb = obj.mesh.n4Nb;
            n4Sb = obj.mesh.n4Sb;
            % count all original dofs and then add hybrid ones, one by one
            curNrDoF = obj.nrDoF-3*obj.hybridNodes;

            if ~isempty(n4Nb)
                n4Nb = obj.mesh.n4Nb;

                NbEd4n = obj.mesh.NbEd4n;
                normal4ed = obj.mesh.normal4ed;
                e4n = obj.mesh.e4n;
                e4ed = obj.mesh.e4ed;
                ed4e = obj.mesh.ed4e;
                c4n = obj.mesh.c4n;

                NbNodes = unique(n4Nb(:));
                NbNodesDof = zeros(nrNodes,2);
                NbNodesDof(NbNodes,:) = reshape(1:2*length(NbNodes),[],2);
                fixNbNodes = [];
                % nodes
                bN = zeros(2*length(NbNodes),1);
                for j = 1 : length(NbNodes)
                    curNbEdges = NbEd4n(NbNodes(j),:);
                    % one Nb one Db
                    if nnz(curNbEdges) < 2
                        if curNbEdges(2) == 0
                            curNormal = normal4ed(curNbEdges(1),:);
                        else
                            curNormal = normal4ed(curNbEdges(2),:);
                        end
                        bN(NbNodesDof(NbNodes(j),:)) = g(c4n(NbNodes(j),:),curNormal).*obj.scaling(NbNodes(j));
                    else % two adjacent Nb Edges
                        curNormals = normal4ed(curNbEdges,:);
                        % Parallel
                        if abs(abs(sum(curNormals(1,:).*curNormals(2,:)))-1) <= eps
                            curNormal = curNormals(1,:);
                            bN(NbNodesDof(NbNodes(j),:)) = g(c4n(NbNodes(j),:),curNormal).*obj.scaling(NbNodes(j));
                        else % non parallel
                            if ~obj.hybrid
                                curG = g([c4n(NbNodes(j),:);c4n(NbNodes(j),:)],curNormals)';
                                curA = [curNormals(1,1),0,curNormals(1,2);
                                    0, curNormals(1,2), curNormals(1,1);
                                    curNormals(2,1),0,curNormals(2,2);
                                    0, curNormals(2,2), curNormals(2,1)];
                                z =  pinv(curA)*curG(:);
                                fixNbNodes = [fixNbNodes,NbNodesDof(NbNodes(j),:)]; %#ok<AGROW>
                                dof = [NbNodes(j),nrNodes+NbNodes(j),2*nrNodes+NbNodes(j)];
                                x(dof) = z;
                            else
                                fixNbNodes = [fixNbNodes,NbNodesDof(NbNodes(j),:)]; %#ok<AGROW>
                                dof = [NbNodes(j),nrNodes+NbNodes(j),2*nrNodes+NbNodes(j)];
                                if numel(find(e4n(NbNodes(j),:))) ~= 2; error('fem4m:AW1:invalidCorner','corner without two triangles');end
                                T = e4ed(curNbEdges,1);
                                % solve local problem
                                curG = g([c4n(NbNodes(j),:);c4n(NbNodes(j),:)],curNormals)';
                                curG = [curG(1:2),0,0,curG(3:4)]';
                                % third normal
                                E3 = intersect(ed4e(T(1),:),ed4e(T(2),:));
                                curNormals(3,:) = normal4ed(E3,:);
                                curA = [curNormals(1,1),0,curNormals(1,2),0,0,0;
                                    0, curNormals(1,2), curNormals(1,1),0,0,0;
                                    curNormals(3,1),0,curNormals(3,2),-curNormals(3,1),0,-curNormals(3,2);
                                    0, curNormals(3,2), curNormals(3,1),0, -curNormals(3,2), -curNormals(3,1);
                                    0,0,0,curNormals(2,1),0,curNormals(2,2);
                                    0,0,0,0, curNormals(2,2), curNormals(2,1)];
                                z =  curA\curG;
                                % count hybrid dofs
                                newDof = curNrDoF+(1:3);
                                curNrDoF = curNrDoF + 3;
                                % assign values
                                x([dof,newDof]) = z;
                            end
                        end
                    end
                end
                freeNbNodes = setdiff(1:length(bN),fixNbNodes);
                bN = bN(freeNbNodes,:);
            end

            if ~isempty(n4Sb)
                bN = [bN; zeros(size(obj.Constraints,1)-numel(bN),1)];
            end

            x = x.*obj.scaling;
        end

        function val = intMoment0(~,points,~,~,g,mesh)
            %intMoment0 - Evaluate the zeroth order moment integrand
            %   VAL = intMoment0(OBJ,POINTS,REFP,IND,G,MESH) evaluates
            %   the traction G(POINTS,NORMALS) at the points POINTS (one
            %   row per Neumann edge) with the normals of the Neumann
            %   edges of MESH, the integrand of the zeroth order edge
            %   moments. REFP and IND are ignored.
            normals = mesh.normal4ed;
            normals = normals(mesh.ed4Nb,:);
            val = g(points,normals);
        end

        function val = intMoment1(~,points,~,~,g,mesh)
            %intMoment1 - Evaluate the first order moment integrand
            %   VAL = intMoment1(OBJ,POINTS,REFP,IND,G,MESH) evaluates
            %   the traction G(POINTS,NORMALS) on the Neumann edges of
            %   MESH weighted by the signed tangential distance
            %   (POINTS-MID)*T of POINTS from the edge midpoints MID, the
            %   integrand of the first order edge moments. REFP and IND
            %   are ignored.
            tangents = mesh.tangent4ed;
            tangents = tangents(mesh.ed4Nb,:);
            normals = mesh.normal4ed;
            normals = normals(mesh.ed4Nb,:);
            mid = mesh.point4ed(1/2);
            mid = mid(mesh.ed4Nb,:);
            curG = g(points,normals);
            val = (sum((points-mid).*tangents,2)*ones(1,2)).*curG;
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of the 3 nodal stress values, 4
            %   edge moments, and 3 element moments.
            obj.nrBasis4n = 3;
            obj.nrBasis4ed = 4;
            obj.nrBasis4e = 3;
        end

        function byEntity = dofOrder(~)
            %dofOrder - Local order of the dofs per entity type
            %   BYENTITY = dofOrder(OBJ) is [true true false]: the nodal
            %   and edge dofs of each element are listed entity by entity.
            byEntity = [true true false];
        end

        %% Organisation of the degrees of freedom
        function obj = dofManager(obj)
            %dofManager - Assemble the stress degrees of freedom
            %   OBJ = dofManager(OBJ) enumerates the nodal, edge, and element
            %   stress degrees of freedom of dofCounts with FE.dofManager
            %   (the local dofs listed entity by entity as given by
            %   dofOrder), fixes all edge moments on the Neumann boundary
            %   and those of symmetryDoFs on the symmetry boundary, and
            %   computes the scaling of each degree of freedom.
            obj.nrBasisFunc = 3*obj.nrBasis4n + 3*obj.nrBasis4ed + obj.nrBasis4e;
            obj = dofManager@FE(obj,false,false,obj.dofOrder());
            % edge moments of each edge
            edgeDoF = obj.dof4ed(:,2*obj.nrBasis4n+1:end);
            obj.fixedNodes = reshape(edgeDoF(obj.mesh.ed4Nb,:),[],1);
            if ~isempty(obj.mesh.n4Sb)
                obj.fixedNodes = [obj.fixedNodes;obj.symmetryDoFs(edgeDoF)];
            end
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
            % scale dofs properly
            obj.scaling = [repmat(sqrt(obj.mesh.area4n),obj.nrBasis4n,1);...
                repmat(sqrt(obj.mesh.area4ed),obj.nrBasis4ed,1);...
                zeros(obj.nrBasis4e*obj.mesh.nrElems,1)];
            obj.scaling(obj.dof4e(:,end-obj.nrBasis4e+1:end)) = sqrt(obj.mesh.area4e)*ones(1,obj.nrBasis4e);
        end

        function dofs = symmetryDoFs(obj,edgeDoF)
            %symmetryDoFs - Edge moments fixed on the symmetry boundary
            %   DOFS = symmetryDoFs(OBJ,EDGEDOF) returns the edge moments
            %   EDGEDOF(ED,J) that vanish on the symmetry boundary edges
            %   ED: the odd columns J (the x components of the moments of
            %   every degree) on the edges with fixed displacement
            %   component 2 (n4Sb(:,3) = 2), where (sigma*n)_x = 0, and
            %   the even columns J (the y components) on the edges with
            %   fixed component 1, where (sigma*n)_y = 0.
            ed4Sb = obj.mesh.ed4Sb;
            flag = obj.mesh.n4Sb(:,3);
            dofs = [];
            for j = 1 : size(edgeDoF,2)
                if mod(j,2) == 1
                    dofs = [dofs;edgeDoF(ed4Sb(flag==2),j)]; %#ok<AGROW>
                else
                    dofs = [dofs;edgeDoF(ed4Sb(flag==1),j)]; %#ok<AGROW>
                end
            end
        end

    end
end