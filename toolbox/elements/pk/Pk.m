classdef Pk < FE
    %PK - Continuous Lagrange Pk element of degree k
    %   OBJ = PK(DEGREE,MESH) creates the H1-conforming scalar
    %   Lagrange finite element of polynomial degree DEGREE >= 1 on
    %   the triangular mesh MESH (a Mesh or a MeshSurface). The
    %   degrees of freedom are the values at the Lagrange points of
    %   the triangles (barycentric coordinates in multiples of
    %   1/DEGREE): one per vertex, DEGREE-1 per edge, and
    %   (DEGREE-1)*(DEGREE-2)/2 in the interior of each triangle.
    %   The basis is defined on the reference triangle. The degrees
    %   of freedom on the Dirichlet boundary are the fixedNodes, all
    %   others the freeNodes.
    %
    %   PK functions:
    %       getBasis              - Evaluate the basis functions
    %       getGradBasis          - Evaluate the gradients of the basis
    %       getD2Basis            - Evaluate the Hessians of the basis
    %       getD3Basis            - Evaluate third derivatives of the basis
    %       getD4Basis            - Evaluate fourth derivatives of the basis
    %       local2global          - Map reference dofs to global points
    %       getMonomBasis         - Build the monomial basis functions
    %       getBasisCoefficients  - Build the dual basis coefficients
    %
    %   PK properties:
    %       mesh                  - Mesh on which the element is defined
    %       nrDoF                 - Total number of degrees of freedom
    %       dof4e                 - Local-to-global dofs per element
    %       dof4ed                - Local-to-global dofs per edge
    %       c4dof                 - Coordinates of the dof points
    %       freeNodes             - Indices of the free degrees of freedom
    %       fixedNodes            - Indices of the fixed degrees of freedom
    %       dof4Db                - Dofs on the Dirichlet boundary
    %       degree                - Polynomial degree of the element
    %       isRefElemBased        - Basis lives on the reference element
    %       monomBasis            - Monomial basis functions
    %       refBasisCoeff         - Coefficients of the basis functions
    %       gradBasisCoefficients - Gradient coefficients of the basis
    %       d2BasisCoefficients   - Hessian coefficients of the basis
    %       d3BasisCoefficients   - Third derivative coefficients
    %       d4BasisCoefficients   - Fourth derivative coefficients
    %       c4refDof              - Reference coordinates of the dof points
    %       refPoints             - Lagrange points on the reference edges
    %       nrBasisFunc           - Number of basis functions per element
    %       nrBasis4n             - Number of basis functions per vertex
    %       nrBasis4ed            - Number of basis functions per edge
    %       nrBasis4e             - Number of interior bubble functions
    %
    %   See also Pk3d, DGPk, PkVector, HPk, CRk, FE, Monomials

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = public)
        mesh                   % mesh on which the element is defined
        nrDoF                  % total number of degrees of freedom
        dof4e                  % local-to-global dofs per element
        dof4ed                 % local-to-global dofs per edge
        c4dof                  % coordinates of the dof points
        freeNodes              % indices of the free degrees of freedom
        fixedNodes             % indices of the fixed degrees of freedom
        dof4Db                 % dofs on the Dirichlet boundary
        degree                 % polynomial degree of the element
        isRefElemBased = true  % basis lives on the reference element
    end

    properties(Access = public)
        c4refDof               % reference coordinates of the dof points
        refPoints              % Lagrange points on the reference edges
        nrBasisFunc            % number of basis functions per element
        nrBasis4n              % number of basis functions per vertex
        nrBasis4ed             % number of basis functions per edge
        nrBasis4e              % number of interior bubble functions
    end

    properties (SetAccess = private)
        monomBasis             % monomial basis functions
        refBasisCoeff          % coefficients of the basis functions
        gradBasisCoefficients  % coefficients of the basis gradients
        d2BasisCoefficients    % coefficients of the basis Hessians
        d3BasisCoefficients    % coefficients of third derivatives
        d4BasisCoefficients    % coefficients of fourth derivatives
    end

    methods
        function obj = Pk(degree,mesh)
            %PK - Create a Pk finite element object
            if degree < 1
                error('fem4m:Pk:invalidDegree','Pk degree must be at least 1');
            end
            obj.degree = degree;
            obj.mesh = mesh;
            obj = getMonomBasis(obj);
            obj = getBasisCoefficients(obj);
            obj = obj.dofCounts();
            obj = obj.dofManager;
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,X,PARTS) evaluates the Pk basis
            %   functions at the reference point X (1-by-2) and
            %   replicates the result for all elements
            %   (nrElems-by-nrBasisFunc). PTS and PARTS are ignored.
            val = obj.refBasisCoeff*obj.monomBasis(x(:,1),x(:,2))';
            val = reshape((val(:)*ones(1,obj.mesh.nrElems))',obj.mesh.nrElems,[]);
        end

        function val = getGradBasis(obj,~,x,~)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,X,PARTS) evaluates the
            %   gradients of the Pk basis functions at the reference
            %   point X (1-by-2), mapped to the physical elements by
            %   the mesh transformation (nrElems-by-nrBasisFunc-by-d,
            %   d = 3 on a MeshSurface). PTS and PARTS are ignored.
            y = x(:,2);
            x = x(:,1);
            m = obj.monomBasis(x,y).';
            grad = permute(reshape(obj.gradBasisCoefficients*m,obj.nrBasisFunc,2),[3 1 2]);
            c = obj.mesh.transformation([x,y]);
            val = matMul(grad,permute(c,[1 3 2]));
        end

        function val = getD2Basis(obj,~,x,~)
            %getD2Basis - Evaluate the Hessians of the basis
            %   VAL = getD2Basis(OBJ,PTS,X,PARTS) evaluates the
            %   Hessians [u_xx,u_yy,u_xy] of the Pk basis functions at
            %   the reference point X (1-by-2), mapped to the physical
            %   elements; on a MeshSurface the components are
            %   [u_xx,u_yy,u_zz,u_xy,u_xz,u_yz]. On a curved surface
            %   mesh (MeshSurface of degree > 1) FE.curvedSurfaceD2
            %   adds the derivative of the transformation applied to
            %   the reference gradients and projects the result onto
            %   the tangent plane. PTS and PARTS are ignored.
            y = x(:,2);
            x = x(:,1);
            m = obj.monomBasis(x,y).';
            d2 = permute(reshape(obj.d2BasisCoefficients*m,obj.nrBasisFunc,3),[3 1 2]);
            c = obj.mesh.transformation([x,y]);
            val = FE.mapDerivatives(c,d2,2);
            if isa(obj.mesh,'MeshSurface') && obj.mesh.degree > 1
                grad = permute(reshape(obj.gradBasisCoefficients*m,obj.nrBasisFunc,2),[3 1 2]);
                val = FE.curvedSurfaceD2(obj.mesh,[x,y],c,val,grad);
            end
        end

        function val = getD3Basis(obj,~,x,~)
            %getD3Basis - Evaluate third derivatives of the basis
            %   VAL = getD3Basis(OBJ,PTS,X,PARTS) evaluates the third
            %   derivatives [u_xxx,u_yyy,u_xxy,u_xyy] of the Pk basis
            %   functions at the reference point X (1-by-2), mapped to
            %   the physical elements. Surface meshes raise the error
            %   fem4m:Pk:surfaceDerivative. PTS and PARTS are ignored.
            if isa(obj.mesh,'MeshSurface')
                error('fem4m:Pk:surfaceDerivative', ...
                    'Third and fourth derivatives of Pk are not supported on surface meshes');
            end
            y = x(:,2);
            x = x(:,1);
            m = obj.monomBasis(x,y).';
            d3 = permute(reshape(obj.d3BasisCoefficients*m,obj.nrBasisFunc,4),[3 1 2]);
            val = FE.mapDerivatives(obj.mesh.transformation([x,y]),d3,3);
        end

        function val = getD4Basis(obj,~,x,~)
            %getD4Basis - Evaluate fourth derivatives of the basis
            %   VAL = getD4Basis(OBJ,PTS,X,PARTS) evaluates the fourth
            %   derivatives [u_xxxx,u_yyyy,u_xxyy,u_xxxy,u_xyyy] of the
            %   Pk basis functions at the reference point X (1-by-2),
            %   mapped to the physical elements. Surface meshes raise
            %   the error fem4m:Pk:surfaceDerivative. PTS and PARTS are
            %   ignored.
            if isa(obj.mesh,'MeshSurface')
                error('fem4m:Pk:surfaceDerivative', ...
                    'Third and fourth derivatives of Pk are not supported on surface meshes');
            end
            y = x(:,2);
            x = x(:,1);
            m = obj.monomBasis(x,y).';
            d4 = permute(reshape(obj.d4BasisCoefficients*m,obj.nrBasisFunc,5),[3 1 2]);
            val = FE.mapDerivatives(obj.mesh.transformation([x,y]),d4,4);
        end


    end

    methods (Access = public)

        function val = local2global(obj,parts,curPoint)
            %local2global - Map reference dofs to global points
            %   VAL = local2global(OBJ,PARTS,CURPOINT) maps the
            %   reference coordinates CURPOINT of one point to the
            %   physical points of the edges or triangles given by the
            %   node numbers PARTS, one row of VAL per row of PARTS.
            %   For edges (two columns) CURPOINT is the scalar t and
            %   VAL = c1 + t*(c2-c1); for triangles (three columns)
            %   CURPOINT = [x y] and VAL = c1 + x*(c2-c1) + y*(c3-c1),
            %   with the coordinates c1, c2, c3 of the nodes in PARTS.
            %1d or 2d
            if size(parts,2) == 2
                %1d ref = conv(0,1)
                c1 = obj.mesh.c4n(parts(:,1),:);
                c2 = obj.mesh.c4n(parts(:,2),:);
                val = c1 + curPoint.*(c2-c1);

            elseif size(parts,2) == 3
                % 2d ref = conv{ (0,0),(1,0),(0,1)}
                c1 = obj.mesh.c4n(parts(:,1),:);
                c2 = obj.mesh.c4n(parts(:,2),:);
                c3 = obj.mesh.c4n(parts(:,3),:);

                val = c1 + curPoint(1)*(c2-c1) + ...
                    curPoint(2)*(c3-c1);
            end
        end

        function obj = getMonomBasis(obj)
            %getMonomBasis - Build the monomial basis functions
            %   OBJ = getMonomBasis(OBJ) sets the number nrBasisFunc
            %   of the monomials x^i*y^j of total degree up to the
            %   degree of OBJ and the handle monomBasis that evaluates them
            %   in the lexicographic order of Monomials.exponents, one
            %   row per point.
            E = Monomials.exponents(2,obj.degree);
            obj.nrBasisFunc = size(E,1);
            obj.monomBasis = Monomials.basis(E);
        end

        function obj = getBasisCoefficients(obj)
            %getBasisCoefficients - Build the dual basis coefficients
            %   OBJ = getBasisCoefficients(OBJ) places the Lagrange
            %   interpolation points on the reference triangle
            %   (c4refDof; refPoints holds the edge parameters j/DEGREE,
            %   j = 1,...,DEGREE-1), computes the coefficients
            %   refBasisCoeff turning monomials into dual basis
            %   functions, and the stacked coefficients
            %   gradBasisCoefficients, d2BasisCoefficients,
            %   d3BasisCoefficients, and d4BasisCoefficients of their
            %   first to fourth derivatives with one block of rows per
            %   derivative component.
            x = zeros(obj.nrBasisFunc,1);
            y = zeros(obj.nrBasisFunc,1);

            % nodal points
            x(1:3) = [ 0 1 0];
            y(1:3) = [0 0 1];

            % Lagrange points
            points = zeros(1,obj.degree-1);
            for j = 1 : (obj.degree-1)
                points(j) = j/obj.degree;
            end
            obj.refPoints = points;

            % edge points
            for i = 1 : (obj.degree-1)
                curPoint = points(i);
                x(4+3*(i-1))   = curPoint;
                y(4+3*(i-1))   = 0;
                x(4+3*(i-1)+1)   = 1-curPoint;
                y(4+3*(i-1)+1)   = curPoint;
                x(4+3*(i-1)+2)   = 0;
                y(4+3*(i-1)+2)   = 1-curPoint;
            end

            % inner Points
            curPoint = 3+3*(obj.degree-1)+1;
            for j = 1 : (obj.degree-1)
                for k = 1 : (obj.degree-1)
                    if points(k)+points(j) < 1
                        x(curPoint) = points(k);
                        y(curPoint) = points(j);
                        curPoint = curPoint+1;
                    end
                end
            end

            % get coefficients
            A = obj.monomBasis(x,y).';
            C = A\eye(obj.nrBasisFunc);
            E = Monomials.exponents(2,obj.degree);
            obj.refBasisCoeff = C;
            obj.gradBasisCoefficients = Monomials.basisCoefficients(C,E,1);
            obj.d2BasisCoefficients = Monomials.basisCoefficients(C,E,2);
            obj.d3BasisCoefficients = Monomials.basisCoefficients(C,E,3);
            obj.d4BasisCoefficients = Monomials.basisCoefficients(C,E,4);
            obj.c4refDof = [x,y];
        end

    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of nodal, edge, and interior
            %   degrees of freedom. The constructor calls it right
            %   before dofManager, so that subclasses may override it.
            obj.nrBasis4n = 1;
            obj.nrBasis4ed = obj.degree-1;
            obj.nrBasis4e = max(0,(obj.degree-1)*(obj.degree-2)/2);
        end

        function obj = dofManager(obj)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ) enumerates the nodal, edge, and
            %   interior degrees of freedom with FE.dofManager for the
            %   counts nrBasis4n, nrBasis4ed, and nrBasis4e set by
            %   dofCounts (edge dofs oriented along n4ed), fixes them
            %   on the Dirichlet boundary, and computes the coordinates
            %   of the dof points.
            obj = dofManager@FE(obj,true);
            obj.fixedNodes = obj.dof4Db;
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);

            nrNodes = obj.mesh.nrNodes;
            nrEdges = obj.mesh.nrEdges;
            n4e = obj.mesh.n4e;
            c4n = obj.mesh.c4n;
            points = obj.refPoints;
            %% c4dof
            obj.c4dof = zeros(obj.nrDoF,size(c4n,2));
            obj.c4dof(1:nrNodes,:) = c4n;
            % edge points
            for j = 1 : (obj.degree-1)
                obj.c4dof(nrNodes+(j-1)*nrEdges+1:nrNodes+j*nrEdges,:) = obj.local2global(obj.mesh.n4ed,points(j));
            end
            % inner Points
            l = 1;
            for j = 1 : (obj.degree-1)
                for k = 1 : (obj.degree-1)
                    if points(k)+points(j) < 1
                        obj.c4dof(obj.dof4e(:,3+3*(obj.degree-1)+l),:) = ...
                            obj.local2global(n4e,[points(k),points(j)]);
                        l = l+1;
                    end
                end
            end
        end

    end

end