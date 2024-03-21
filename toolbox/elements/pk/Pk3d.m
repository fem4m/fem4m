classdef Pk3d < FE
    %Pk3d - Continuous Lagrange Pk element of degree k on tetrahedra
    %   OBJ = Pk3d(DEGREE,MESH) creates the H1-conforming scalar
    %   Lagrange finite element of polynomial degree DEGREE >= 1 on
    %   the tetrahedral mesh MESH. The degrees of freedom are the
    %   values at the Lagrange points of the tetrahedra (barycentric
    %   coordinates in multiples of 1/DEGREE): one per vertex,
    %   DEGREE-1 per edge, (DEGREE-1)*(DEGREE-2)/2 per face, and
    %   (DEGREE-1)*(DEGREE-2)*(DEGREE-3)/6 in the interior of each
    %   tetrahedron. The basis is defined on the reference
    %   tetrahedron. The degrees of freedom on the Dirichlet
    %   boundary are the fixedNodes, all others the freeNodes.
    %
    %   Pk3d functions:
    %       getBasis              - Evaluate the basis functions
    %       getGradBasis          - Evaluate the gradients of the basis
    %       getD2Basis            - Evaluate the Hessians of the basis
    %       getD3Basis            - Evaluate third derivatives of the basis
    %       getD4Basis            - Evaluate fourth derivatives of the basis
    %       local2global          - Map reference dofs to global points
    %
    %   Pk3d properties:
    %       nrDoF                 - Total number of degrees of freedom
    %       dof4e                 - Local-to-global dofs per element
    %       dof4f                 - Local-to-global dofs per face
    %       dof4ed                - Local-to-global dofs per edge
    %       c4dof                 - Coordinates of the dof points
    %       freeNodes             - Indices of the free degrees of freedom
    %       fixedNodes            - Indices of the fixed degrees of freedom
    %       dof4Db                - Dofs on the Dirichlet boundary
    %       degree                - Polynomial degree of the element
    %       mesh                  - Mesh on which the element is defined
    %       nrBasisFunc           - Number of basis functions per element
    %       nrBasis4n             - Number of basis functions per vertex
    %       nrBasis4ed            - Number of basis functions per edge
    %       nrBasis4f             - Number of basis functions per face
    %       nrBasis4e             - Number of interior bubble functions
    %       isRefElemBased        - Basis lives on the reference element
    %       monomBasis            - Monomial basis functions
    %       refBasisCoeff         - Coefficients of the basis functions
    %       gradBasisCoefficients - Gradient coefficients of the basis
    %       d2BasisCoefficients   - Hessian coefficients of the basis
    %       d3BasisCoefficients   - Third derivative coefficients
    %       d4BasisCoefficients   - Fourth derivative coefficients
    %       c4refDof              - Reference coordinates of the dof points
    %       refPoints             - Not set (always empty)
    %       nrNodePoints          - Number of interpolation points per vertex
    %       nrEdgePoints          - Number of interpolation points per edge
    %       nrFacePoints          - Number of interpolation points per face
    %       nrInnerPoints         - Number of interior interpolation points
    %       fperm                 - Face permutations for dof transformations
    %       findices              - Face index tables for dof transformations
    %       transformation        - Inverse Jacobians of the element maps
    %
    %   See also Pk, DGPk3d, HPk3d, PkVector, FE, Monomials

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = public)
        nrDoF          % total number of degrees of freedom
        dof4e          % local-to-global dofs per element
        dof4f          % local-to-global dofs per face
        dof4ed         % local-to-global dofs per edge
        c4dof          % coordinates of the dof points
        freeNodes      % indices of the free degrees of freedom
        fixedNodes     % indices of the fixed degrees of freedom
        dof4Db         % dofs on the Dirichlet boundary
        degree         % polynomial degree of the element
        mesh           % mesh on which the element is defined
        nrBasisFunc    % number of basis functions per element
        nrBasis4n      % number of basis functions per vertex
        nrBasis4ed     % number of basis functions per edge
        nrBasis4f      % number of basis functions per face
        nrBasis4e      % number of interior bubble functions
        isRefElemBased = true % basis lives on the reference element
    end

    properties(Access = public)
        c4refDof       % reference coordinates of the dof points
        refPoints      % not set (always empty)
        nrNodePoints   % number of interpolation points per vertex
        nrEdgePoints   % number of interpolation points per edge
        nrFacePoints   % number of interpolation points per face
        nrInnerPoints  % number of interior interpolation points
        fperm          % face permutations for dof transformations
        findices       % face index tables for dof transformations
        transformation % inverse Jacobians of the element maps
    end

    properties (SetAccess = private)
        monomBasis            % monomial basis functions
        refBasisCoeff         % coefficients of the basis functions
        gradBasisCoefficients % coefficients of the basis gradients
        d2BasisCoefficients   % coefficients of the basis Hessians
        d3BasisCoefficients   % coefficients of third derivatives
        d4BasisCoefficients   % coefficients of fourth derivatives
    end

    methods
        function obj = Pk3d(degree,mesh)
            %Pk3d - Create a Pk3d finite element object
            if degree < 1
                error('fem4m:Pk3d:invalidDegree','Pk3d degree must be at least 1');
            end
            obj.degree = degree;
            obj = obj.dofCounts();
            obj.nrBasisFunc = 4*obj.nrBasis4n + 6*obj.nrBasis4ed + 4*obj.nrBasis4f + obj.nrBasis4e;
            obj = getMonomBasis(obj);
            obj = getBasisCoefficients(obj);
            obj = obj.dofManager(mesh);
            obj.transformation = mesh.transformation;
            obj.mesh = mesh;
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,X,PARTS) evaluates the Pk3d
            %   basis functions at the reference point X (1-by-3) and
            %   replicates the result for all elements
            %   (nrElems-by-nrBasisFunc). PTS and PARTS are ignored.
            val = obj.refBasisCoeff*obj.monomBasis(x(1),x(2),x(3))';
            val = reshape((val(:)*ones(1,size(obj.transformation,1)))',size(obj.transformation,1),[]);
        end

        function val = getGradBasis(obj,~,x,~)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,X,PARTS) evaluates the
            %   gradients of the Pk3d basis at the reference point X
            %   (1-by-3), mapped to the physical elements by the saved
            %   transformation (nrElems-by-nrBasisFunc-by-3). PTS and
            %   PARTS are ignored.
            m = obj.monomBasis(x(1),x(2),x(3)).';
            grad = permute(reshape(obj.gradBasisCoefficients*m,obj.nrBasisFunc,3),[3 1 2]);
            val = FE.mapDerivatives(obj.transformation,grad,1);
        end

        function val = getD2Basis(obj,~,x,~)
            %getD2Basis - Evaluate the Hessians of the basis
            %   VAL = getD2Basis(OBJ,PTS,X,PARTS) evaluates the
            %   Hessians [u_xx,u_yy,u_zz,u_xy,u_xz,u_yz] of the Pk3d
            %   basis at the reference point X (1-by-3), mapped to the
            %   physical elements. PTS and PARTS are ignored.
            m = obj.monomBasis(x(1),x(2),x(3)).';
            d2 = permute(reshape(obj.d2BasisCoefficients*m,obj.nrBasisFunc,6),[3 1 2]);
            val = FE.mapDerivatives(obj.transformation,d2,2);
        end

        function val = getD3Basis(obj,~,x,~)
            %getD3Basis - Evaluate third derivatives of the basis
            %   VAL = getD3Basis(OBJ,PTS,X,PARTS) evaluates the third
            %   derivatives [xxx,yyy,zzz,xxy,xxz,xyy,yyz,xzz,yzz,xyz]
            %   of the Pk3d basis at the reference point X (1-by-3),
            %   mapped to the physical elements. PTS and PARTS are
            %   ignored.
            m = obj.monomBasis(x(1),x(2),x(3)).';
            d3 = permute(reshape(obj.d3BasisCoefficients*m,obj.nrBasisFunc,10),[3 1 2]);
            val = FE.mapDerivatives(obj.transformation,d3,3);
        end

        function val = getD4Basis(obj,~,x,~)
            %getD4Basis - Evaluate fourth derivatives of the basis
            %   VAL = getD4Basis(OBJ,PTS,X,PARTS) evaluates the 15
            %   fourth derivatives [xxxx,yyyy,zzzz,xxxy,xxxz,xyyy,yyyz,
            %   xzzz,yzzz,xxyy,xxzz,yyzz,xxyz,xyyz,xyzz] of the Pk3d
            %   basis at the reference point X (1-by-3), mapped to the
            %   physical elements. PTS and PARTS are ignored.
            m = obj.monomBasis(x(1),x(2),x(3)).';
            d4 = permute(reshape(obj.d4BasisCoefficients*m,obj.nrBasisFunc,15),[3 1 2]);
            val = FE.mapDerivatives(obj.transformation,d4,4);
        end

    end

    methods (Access = private)

        function obj = getBasisCoefficients(obj)
            %getBasisCoefficients - Build the dual basis coefficients

            % Lagrange points 1d
            points = zeros(1,obj.degree-1);
            for j = 1 : (obj.degree-1)
                points(j) = j/obj.degree;
            end

            x = zeros(obj.nrBasisFunc,1);
            y = zeros(obj.nrBasisFunc,1);
            z = zeros(obj.nrBasisFunc,1);

            % nodal points
            x(1:4) = [ 0 1 0 0];
            y(1:4) = [ 0 0 1 0];
            z(1:4) = [ 0 0 0 1];
            obj.nrNodePoints = 4;

            % edge points
            for i = 1 : (obj.degree-1)
                curPoint = points(i);
                x(5+6*(i-1))   = curPoint;
                y(5+6*(i-1))   = 0;
                z(5+6*(i-1))   = 0;
                x(5+6*(i-1)+1)   = 1-curPoint;
                y(5+6*(i-1)+1)   = curPoint;
                z(5+6*(i-1)+1)   = 0;
                x(5+6*(i-1)+2)   = 0;
                y(5+6*(i-1)+2)   = 1-curPoint;
                z(5+6*(i-1)+2)   = 0;
                x(5+6*(i-1)+3)   = 0;
                y(5+6*(i-1)+3)   = 0;
                z(5+6*(i-1)+3)   = curPoint;
                x(5+6*(i-1)+4)   = 1-curPoint;
                y(5+6*(i-1)+4)   = 0;
                z(5+6*(i-1)+4)   = curPoint;
                x(5+6*(i-1)+5)   = 0;
                y(5+6*(i-1)+5)   = 1-curPoint;
                z(5+6*(i-1)+5)   = curPoint;
            end
            obj.nrEdgePoints = 6*(obj.degree-1);

            % Face Points
            curPoint = obj.nrNodePoints+obj.nrEdgePoints+1;
            obj.nrFacePoints = 0;
            for j = 1 : (obj.degree-1)
                for k = 1 : (obj.degree-1)
                    if points(k)+points(j) < 1
                        x(curPoint) = points(k);
                        y(curPoint) = points(j);
                        z(curPoint) = 0;
                        curPoint = curPoint+1;
                        x(curPoint) = points(k);
                        y(curPoint) = 0;
                        z(curPoint) = points(j);
                        curPoint = curPoint+1;
                        x(curPoint) = points(k);
                        y(curPoint) = points(j);
                        z(curPoint) = 1-points(k)-points(j);
                        curPoint = curPoint+1;
                        x(curPoint) = 0;
                        y(curPoint) = points(k);
                        z(curPoint) = points(j);
                        curPoint = curPoint+1;
                        obj.nrFacePoints = obj.nrFacePoints + 4;
                    end
                end
            end

            %inner Points
            obj.nrInnerPoints = 0;
            for j = 1 : (obj.degree-1)
                for k = 1 : (obj.degree-1)
                    for l = 1 : (obj.degree-1)
                        if points(k)+points(j)+points(l) < 1-eps
                            x(curPoint) = points(k);
                            y(curPoint) = points(j);
                            z(curPoint) = points(l);
                            curPoint = curPoint+1;
                            obj.nrInnerPoints = obj.nrInnerPoints + 1;
                        end
                    end
                end
            end

            % get coefficients
            A = obj.monomBasis(x,y,z).';
            obj.refBasisCoeff = A\eye(obj.nrBasisFunc);
            obj.c4refDof = [x,y,z];

            % coefficients of the derivatives in the monomials
            C = obj.refBasisCoeff;
            E = Monomials.exponents(3,obj.degree);
            obj.gradBasisCoefficients = Monomials.basisCoefficients(C,E,1);
            obj.d2BasisCoefficients = Monomials.basisCoefficients(C,E,2);
            obj.d3BasisCoefficients = Monomials.basisCoefficients(C,E,3);
            obj.d4BasisCoefficients = Monomials.basisCoefficients(C,E,4);

            obj.nrNodePoints = obj.nrNodePoints/4;
            obj.nrEdgePoints = obj.nrEdgePoints/6;
            obj.nrFacePoints = obj.nrFacePoints/4;
            obj.nrInnerPoints = obj.nrInnerPoints;

            % transformation of dof on faces
            dof = zeros(obj.degree-2,obj.degree-2);
            j=1;
            for k = obj.degree-2:-1:1
                for l = 1 : k
                    dof(k,l) = j;
                    j = j+1;
                end
            end
            perm = zeros(6,3);
            indices = zeros(6,nnz(dof));
            perm(1,:) = [1 2 3];
            indices(1,:) = 1:nnz(dof);

            perm(2,:)=[3 1 2];
            j = 1;
            for l = obj.degree-2:-1:1
                for k = obj.degree-2:-1:1
                    if l+k-obj.degree+2 > 0
                        indices(2,j) = dof(k,l+k-obj.degree+2);
                        j = j+1;
                    end
                end
            end

            perm(3,:)=[2 3 1];
            j = 1;
            for l = 1:obj.degree-2
                for k = l:obj.degree-2
                    indices(3,j) = dof(k,l);
                    j = j+1;
                end
            end

            perm(4,:)=[2 1 3];
            j = 1;
            for k = obj.degree-2:-1:1
                for l = k:-1:1
                    indices(4,j) = dof(k,l);
                    j = j+1;
                end
            end

            perm(5,:)=[1 3 2];
            j = 1;
            for l = 1:obj.degree-2
                for k = obj.degree-2:-1:l
                    indices(5,j) = dof(k,l);
                    j = j+1;
                end
            end

            perm(6,:)=[3 2 1];
            j = 1;
            for l = 1:obj.degree-2
                for k = l:obj.degree-2
                    indices(6,j) = dof(k,k-l+1);
                    j = j+1;
                end
            end

            obj.fperm = perm;
            obj.findices = indices;

        end

        function obj = getMonomBasis(obj)
            %getMonomBasis - Build the monomial basis functions
            %   OBJ = getMonomBasis(OBJ) sets the handle monomBasis of
            %   the monomials x^i*y^j*z^k of total degree up to the
            %   degree of OBJ in the lexicographic order of
            %   Monomials.exponents, one row per point.
            obj.monomBasis = Monomials.basis(Monomials.exponents(3,obj.degree));
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, nrBasis4f, and nrBasis4e of nodal, edge,
            %   face, and interior degrees of freedom. The constructor
            %   calls it first and computes nrBasisFunc from these
            %   counts, so that subclasses may override it.
            obj.nrBasis4n = 1;
            obj.nrBasis4ed = obj.degree-1;
            obj.nrBasis4f = (obj.degree-1)*max(0,obj.degree-2)/2;
            obj.nrBasis4e = (obj.degree-1)*max(0,obj.degree-2)*max(0,obj.degree-3)/6;
        end

        function obj = dofManager(obj,mesh)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ,MESH) enumerates the nodal, edge,
            %   face, and interior degrees of freedom with FE.dofManager
            %   (edge dofs oriented along n4ed, face dofs by fperm and
            %   findices), fixes them on the Dirichlet boundary, and
            %   computes the coordinates of the dof points and the face
            %   connectivity dof4f.
            obj.mesh = mesh;
            obj = dofManager@FE(obj,true,true);
            obj.fixedNodes = obj.dof4Db;
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);

            n4e = mesh.n4e;
            c4n = mesh.c4n;
            nrNodes = mesh.nrNodes;
            nrEdges = mesh.nrEdges;
            nrFaces = mesh.nrFaces;
            ed4f = mesh.ed4f;
            n4f = mesh.n4f;
            n4ed = mesh.n4ed;
            e4f = mesh.e4f;

            %% c4dof
            obj.c4dof = zeros(obj.nrDoF,3);

            obj.c4dof(1:nrNodes,:) = c4n;

            % edge points
            for k = 1 : obj.nrEdgePoints
                obj.c4dof(nrNodes+(k-1)*nrEdges+1:nrNodes+k*nrEdges,:) = ...
                    obj.local2global(n4ed,obj.c4refDof(4*obj.nrNodePoints+(k-1)*6+1,:),mesh);
            end

            % face points
            for k = 1 : obj.nrFacePoints
                obj.c4dof(nrNodes+obj.nrEdgePoints*nrEdges+(k-1)*nrFaces+1:nrNodes+obj.nrEdgePoints*nrEdges+k*nrFaces,:) = ...
                    obj.local2global(n4f,obj.c4refDof(4*obj.nrNodePoints+6*obj.nrEdgePoints+(k-1)*4+1,:),mesh);
            end

            % inner Points
            for k = 1 : obj.nrInnerPoints
                obj.c4dof(obj.dof4e(:,4*obj.nrNodePoints+6*obj.nrEdgePoints+4*obj.nrFacePoints+k),:) = ...
                    obj.local2global(n4e,obj.c4refDof(4*obj.nrNodePoints+6*obj.nrEdgePoints+4*obj.nrFacePoints+k,:),mesh);
            end

            %% dof4f
            dof4n = n4f;
            d4ed = zeros(nrFaces,3*obj.nrEdgePoints);
            if ~isempty(d4ed) && obj.degree == 2
                d4ed = nrNodes+mesh.ed4f;
            elseif ~isempty(d4ed)
                indices = (nrNodes+reshape(1:obj.nrEdgePoints*nrEdges,nrEdges,[]));
                edges = n4f(:,[1 2 2 3 3 1]);
                for j = 1 : 3
                    indx = edges(:,(j-1)*2+1) < edges(:,j*2);
                    d4ed(indx,j:3:3*obj.nrEdgePoints) = indices(ed4f(indx,j),:);
                    d4ed(~indx,j:3:3*obj.nrEdgePoints) = indices(ed4f(~indx,j),end:-1:1);
                end
            end

            d4f = zeros(nrFaces,obj.nrFacePoints);
            if ~isempty(d4f) && obj.degree == 3
                d4f = nrNodes+obj.nrEdgePoints*nrEdges+(1:nrFaces)';
            elseif ~isempty(d4f)
                indices = nrNodes+obj.nrEdgePoints*nrEdges+reshape(1:obj.nrFacePoints*nrFaces,nrFaces,[]);
                [~,sortIdx] = sort(n4f,2);
                [~,perm] = ismember(sortIdx,obj.fperm,'rows');
                for j = 1 : 4
                    indx = mesh.fNr4f(:,1) == j;
                    d4f(indx,:) = ...
                        rowaddr(indices,find(indx)*ones(1,size(obj.findices,2)),obj.findices(perm(e4f(indx,1)),:));
                end
            end
            obj.dof4f = [dof4n,d4ed,d4f];
        end

    end

    methods(Static)

        function val = local2global(parts,curPoint,mesh)
            %local2global - Map reference dofs to global points
            %   VAL = local2global(PARTS,CURPOINT,MESH) maps the
            %   reference coordinates CURPOINT of one point to the
            %   physical points of the edges, triangles, or tetrahedra
            %   given by the node numbers PARTS of the mesh MESH, one
            %   row of VAL per row of PARTS: VAL = c1 +
            %   CURPOINT(1)*(c2-c1) + ... + CURPOINT(m-1)*(cm-c1) for
            %   the coordinates c1,...,cm of the m = 2, 3, or 4 nodes
            %   in a row of PARTS; further entries of CURPOINT are
            %   ignored.
            %1d or 2d
            if size(parts,2) == 2
                %1d ref = conv(0,1)
                c1 = mesh.c4n(parts(:,1),:);
                c2 = mesh.c4n(parts(:,2),:);
                val = c1 + curPoint(1).*(c2-c1);

            elseif size(parts,2) == 3
                % 2d ref = conv{ (0,0),(1,0),(0,1)}
                c1 = mesh.c4n(parts(:,1),:);
                c2 = mesh.c4n(parts(:,2),:);
                c3 = mesh.c4n(parts(:,3),:);

                val = c1 + curPoint(1)*(c2-c1) + ...
                    curPoint(2)*(c3-c1);
            elseif size(parts,2) == 4
                % 3d ref = conv{(0,0,0),(1,0,0),(0,1,0),(0,0,1)}
                c1 = mesh.c4n(parts(:,1),:);
                c2 = mesh.c4n(parts(:,2),:);
                c3 = mesh.c4n(parts(:,3),:);
                c4 = mesh.c4n(parts(:,4),:);

                val = c1 + curPoint(1)*(c2-c1) + ...
                    curPoint(2)*(c3-c1) + ...
                    curPoint(3)*(c4-c1);
            end
        end
    end

end