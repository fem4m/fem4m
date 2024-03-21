classdef CRk < FE
    %CRk - Crouzeix-Raviart element of degree 1 or 3
    %   OBJ = CRk(DEGREE,MESH) creates the nonconforming
    %   Crouzeix-Raviart-type finite element of degree DEGREE on the
    %   triangular mesh MESH. For DEGREE 1 the degrees of freedom are
    %   the mean values on the edges; for DEGREE 3 they are the mean
    %   value and the first and second moments on every edge and the
    %   mean value on every element (three per edge, one per element).
    %   The basis functions are the dual basis of these degrees of
    %   freedom, expressed in the barycentric P3 functions lambda_i,
    %   lambda_i*lambda_j, lambda_i*lambda_j*(lambda_i-lambda_j) and
    %   lambda_1*lambda_2*lambda_3 (basisCoefficients). DEGREE 2 is
    %   rejected because the edge moments are not unisolvent for
    %   quadratic polynomials, and degrees above 3 are not implemented.
    %
    %   CRk functions:
    %       getBasis            - Evaluate the basis functions
    %       getGradBasis        - Evaluate the gradients of the basis
    %       getD2Basis          - Evaluate the Hessians of the basis
    %       interpolateBoundary - Interpolate data on the boundary
    %
    %   CRk properties:
    %       mesh              - Mesh on which the element is defined
    %       nrDoF             - Total number of degrees of freedom
    %       dof4e             - Local-to-global dofs per element
    %       dof4ed            - Local-to-global dofs per edge
    %       degree            - Polynomial degree of the element
    %       nrBasisFunc       - Number of basis functions per element
    %       fixedNodes        - Indices of the fixed degrees of freedom
    %       freeNodes         - Indices of the free degrees of freedom
    %       dof4Db            - Dofs on the Dirichlet boundary
    %       nrBasis4n         - Number of nodal basis functions
    %       nrBasis4ed        - Number of edge basis functions
    %       nrBasis4e         - Number of interior basis functions
    %       isRefElemBased    - Basis lives on the reference element
    %       basisCoefficients - Basis coefficients in barycentric P3 functions
    %
    %   See also Pk, Bk, HPk, FE

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = public)
        mesh              % mesh on which the element is defined
        nrDoF             % total number of degrees of freedom
        dof4e             % local-to-global dofs per element
        dof4ed            % local-to-global dofs per edge
        degree            % polynomial degree of the element
        nrBasisFunc       % number of basis functions per element
        fixedNodes        % indices of the fixed degrees of freedom
        freeNodes         % indices of the free degrees of freedom
        dof4Db            % dofs on the Dirichlet boundary
        nrBasis4n         % number of nodal basis functions
        nrBasis4ed        % number of edge basis functions
        nrBasis4e         % number of interior basis functions
        isRefElemBased = true % basis lives on the reference element
    end

    properties(Access = public)
        basisCoefficients % basis coefficients in barycentric P3 functions
    end

    methods
        function obj = CRk(degree,mesh)
            %CRk - Create a Crouzeix-Raviart finite element object
            obj.degree = degree;
            if degree < 1
                error('fem4m:CRk:invalidDegree','CR degree must be at least 1');
            elseif degree == 2
                error('fem4m:CRk:invalidDegree','CR of degree 2 is not unisolvent');
            elseif degree > 3
                error('fem4m:CRk:degreeNotImplemented','CR of degree > 3 is not implemented');
            end
            obj.mesh = mesh;
            obj.nrBasisFunc = (3+(obj.degree-1)*3 + (obj.degree==3));
            obj = getBasisCoefficients(obj);
            obj = obj.dofCounts();
            obj = obj.dofManager;
        end

        function val = getBasis(obj,~,pts,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the basis
            %   functions in all elements at the reference point PTS_REF
            %   (1-by-2), one row per element; PTS and PARTS are
            %   ignored. A PTS_REF with exactly three rows is taken as
            %   one point per local edge: row j is evaluated on the edges
            %   with the local number j in their first element
            %   e4ed(:,1), and VAL has one row per edge.
            x = pts(:,1);
            y = pts(:,2);
            if length(x) == 3
                e4ed = obj.mesh.e4ed;
                edNr4ed = obj.mesh.edNr4ed;
                ind1 = find(edNr4ed(:,1) == 1);
                ind2 = find(edNr4ed(:,1) == 2);
                ind3 = find(edNr4ed(:,1) == 3);
                val = zeros(size(edNr4ed,1),obj.nrBasisFunc);
                C = obj.basisCoefficients;
                basisP3 = obj.basisP3(x(1),y(1));
                basisPk = basisP3(1:obj.nrBasisFunc);
                basisPk = (basisPk(:)*ones(1,obj.mesh.nrElems))';
                val1 = matMul(C,basisPk);
                val(ind1,:) = val1(e4ed(ind1,1),:);
                basisP3 = obj.basisP3(x(2),y(2));
                basisPk = basisP3(1:obj.nrBasisFunc);
                basisPk = (basisPk(:)*ones(1,obj.mesh.nrElems))';
                val2 = matMul(C,basisPk);
                val(ind2,:) = val2(e4ed(ind2,1),:);
                basisP3 = obj.basisP3(x(3),y(3));
                basisPk = basisP3(1:obj.nrBasisFunc);
                basisPk = (basisPk(:)*ones(1,obj.mesh.nrElems))';
                val3 = matMul(C,basisPk);
                val(ind3,:) = val3(e4ed(ind3,1),:);
            else
                basisP3 = obj.basisP3(x,y);
                basisPk = basisP3(1:obj.nrBasisFunc);
                basisPk = (basisPk(:)*ones(1,obj.mesh.nrElems))';
                C = obj.basisCoefficients;
                val = matMul(C,basisPk);
            end
        end

        function val = getGradBasis(obj,~,pts_ref,~)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   gradients of the basis functions in all elements at the
            %   reference point PTS_REF (1-by-2), mapped to the physical
            %   elements, one row per element and the components [x y]
            %   along the third dimension; PTS and PARTS are ignored. A
            %   PTS_REF with exactly three rows is taken as one point per
            %   local edge as in getBasis, and VAL has one row per edge.

            x_ref = pts_ref(:,1);
            y_ref = pts_ref(:,2);
            area4e = obj.mesh.area4e;
            c4n = obj.mesh.c4n;
            n4e = obj.mesh.n4e;
            edNr4ed = obj.mesh.edNr4ed;

            if length(x_ref) == 3
                parts = obj.mesh.e4ed(:,1);
                x = x_ref(edNr4ed(:,1));
                y = y_ref(edNr4ed(:,1));
            else
                parts = 1:obj.mesh.nrElems;
                x = x_ref;
                y = y_ref;
            end

            gradP3 = obj.gradP3(x,y);
            gradPk = gradP3(1:obj.nrBasisFunc,:);
            gradPk = reshape((gradPk(:)*ones(1,obj.mesh.nrElems))',obj.mesh.nrElems,[],2);
            C = obj.basisCoefficients;
            grad = matMul(C,gradPk);

            P1 = c4n(n4e(parts,1),:);
            P2 = c4n(n4e(parts,2),:);
            P3 = c4n(n4e(parts,3),:);

            ax = (P3(:,2)-P1(:,2))./(2*area4e);
            bx = (-P2(:,2)+P1(:,2))./(2*area4e);
            ay = (-P3(:,1)+P1(:,1))./(2*area4e);
            by = (P2(:,1)-P1(:,1))./(2*area4e);

            val = zeros(length(ax),size(grad,2),size(grad,3));
            for i = 1 : size(grad,2)
                val(:,i,1) = grad(:,i,1).*ax+grad(:,i,2).*bx;
                val(:,i,2) = grad(:,i,1).*ay+grad(:,i,2).*by;
            end
        end

        function val = getD2Basis(obj,~,~,~)
            %getD2Basis - Evaluate the Hessians of the basis
            %   VAL = getD2Basis(OBJ,PTS,PTS_REF,PARTS) returns the zero
            %   Hessians [xx yy xy] of the basis functions for DEGREE 1,
            %   an nrElems-by-3-by-3 array, and raises the error
            %   fem4m:CRk:d2BasisNotImplemented for DEGREE 3. PTS,
            %   PTS_REF and PARTS are ignored.
            if obj.degree == 1
                val = zeros(obj.mesh.nrElems,3,3);
            else
                error('fem4m:CRk:d2BasisNotImplemented','Not implemented');
            end
        end

        function val = interpolateBoundary(obj,u4Db)
            %interpolateBoundary - Interpolate data on the boundary
            %   VAL = interpolateBoundary(OBJ,U4DB) returns the vector of
            %   length nrDoF whose entries at the degrees of freedom of
            %   the Dirichlet edges are those of the boundary data U4DB,
            %   a function handle @(x) of the points x (N-by-2): the
            %   mean values of U4DB on the edges and, for DEGREE 3, the
            %   mean values of s*U4DB and s^2*U4DB with s = (x-m)*t/|E|
            %   (m the midpoint, t the unit tangent and |E| the length
            %   of the edge). All other entries are zero.
            val = zeros(obj.nrDoF,1);
            length4ed = obj.mesh.length4ed;
            ed4Db = obj.mesh.ed4Db;
            n4Db = obj.mesh.n4Db;

            moments = 1./length4ed(ed4Db).*integrate(...
                @(points,refp,ind)obj.intMoment0(points,refp,ind,u4Db,obj.mesh),n4Db,ed4Db,obj.mesh,4);
            if obj.degree > 1
                moments1 = 1./length4ed(ed4Db).^2.*integrate(...
                    @(points,refp,ind)obj.intMoment1(points,refp,ind,u4Db,obj.mesh),n4Db,ed4Db,obj.mesh,5);
                moments = [moments,moments1];
            end
            if obj.degree > 2
                moments2 = 1./length4ed(ed4Db).^3.*integrate(...
                    @(points,refp,ind)obj.intMoment2(points,refp,ind,u4Db,obj.mesh),n4Db,ed4Db,obj.mesh,6);
                moments = [moments,moments2];
            end
            val(obj.dof4ed(ed4Db,1:size(moments,2))) = moments;
        end

    end

    methods (Access = private)

        function obj = getBasisCoefficients(obj)
            %getBasisCoefficients - Build the dual basis coefficients

            nrElems = obj.mesh.nrElems;
            ed4e = obj.mesh.ed4e;
            n4e = obj.mesh.n4e;
            n4ed = obj.mesh.n4ed;
            % Initialisation of the basis coefficients
            obj.basisCoefficients = zeros(nrElems,obj.nrBasisFunc,obj.nrBasisFunc);
            C = zeros(10,10);
            % Element degree of freedom (mean value)
            C(10,:) = [40,40,40,10,10,10,0,0,0,2];
            % Edge degrees of freedom (mean value, second moment)
            C(1:3,:) = [60, 60, 0, 20, 0, 0, 0, 0, 0, 0;
                0, 60, 60, 0, 20, 0, 0, 0, 0, 0;
                60, 0, 60, 0, 0, 20, 0, 0, 0, 0];
            C(7:9,:) = [5, 5, 0, 1, 0, 0, 0, 0, 0, 0;
                0, 5, 5, 0, 1, 0, 0, 0, 0, 0;
                5, 0, 5, 0, 0, 1, 0, 0, 0, 0];
            % First edge moments, signed by the edge orientation
            for curElem = 1 : nrElems
                sig1 = -1; sig2 = -1; sig3 = -1;
                if(n4e(curElem,1)==n4ed(ed4e(curElem,1),1));sig1 = 1;end
                if(n4e(curElem,2)==n4ed(ed4e(curElem,2),1));sig2 = 1;end
                if(n4e(curElem,3)==n4ed(ed4e(curElem,3),1));sig3 = 1;end
                C(4,:) = [-10, 10, 0, 0, 0, 0, -2, 0, 0, 0]*sig1;
                C(5,:) = [0, -10, 10, 0, 0, 0, 0, -2, 0, 0]*sig2;
                C(6,:) = [10, 0, -10, 0, 0, 0, 0, 0, -2, 0]*sig3;
                % Dual basis: invert the matrix of the dofs
                x = (C(1:obj.nrBasisFunc,1:obj.nrBasisFunc)/120)\eye(obj.nrBasisFunc);
                obj.basisCoefficients(curElem,:,:) = x';
            end

        end

        %% P3 basis functions
        function val = basisP3(obj,x,y)
            %basisP3 - Evaluate the barycentric P3 functions
            basisP1 = obj.basisP1(x,y);
            val = zeros(size(x,1),10);
            val(:, 1) = basisP1(:,1);
            val(:, 2) = basisP1(:,2);
            val(:, 3) = basisP1(:,3);
            val(:, 4) = basisP1(:,1).*basisP1(:,2);
            val(:, 5) = basisP1(:,2).*basisP1(:,3);
            val(:, 6) = basisP1(:,3).*basisP1(:,1);
            val(:, 7) = basisP1(:,1).*basisP1(:,2).*(basisP1(:,1) - basisP1(:,2));
            val(:, 8) = basisP1(:,2).*basisP1(:,3).*(basisP1(:,2) - basisP1(:,3));
            val(:, 9) = basisP1(:,3).*basisP1(:,1).*(basisP1(:,3) - basisP1(:,1));
            val(:,10) = basisP1(:,1).*basisP1(:,2).*basisP1(:,3);
        end

        function val = gradP3(obj,x,y)
            %gradP3 - Evaluate gradients of the barycentric P3 functions
            gradP1 = obj.gradP1(x,y);
            basisP1 = obj.basisP1(x,y);

            val = zeros(10,2);
            val(1,:) = gradP1(:,1,:);
            val(2,:) = gradP1(:,2,:);
            val(3,:) = gradP1(:,3,:);
            val(4,:) = gradP1(:,1,:)*basisP1(:,2) + gradP1(:,2,:)*basisP1(:,1);
            val(5,:) = gradP1(:,2,:)*basisP1(:,3) + gradP1(:,3,:)*basisP1(:,2);
            val(6,:) = gradP1(:,3,:)*basisP1(:,1) + gradP1(:,1,:)*basisP1(:,3);
            val(7,:) = gradP1(:,1,:)*basisP1(:,2)*(basisP1(:,1) - basisP1(:,2)) + ...
                basisP1(:,1)*gradP1(:,2,:)*(basisP1(:,1) - basisP1(:,2)) + ...
                basisP1(:,1)*basisP1(:,2)*(gradP1(:,1,:) - gradP1(:,2,:));
            val(8,:) = gradP1(:,2,:)*basisP1(:,3)*(basisP1(:,2) - basisP1(:,3)) + ...
                basisP1(:,2)*gradP1(:,3,:)*(basisP1(:,2) - basisP1(:,3)) + ...
                basisP1(:,2)*basisP1(:,3)*(gradP1(:,2,:) - gradP1(:,3,:));
            val(9,:) = gradP1(:,3,:)*basisP1(:,1)*(basisP1(:,3) - basisP1(:,1)) + ...
                basisP1(:,3)*gradP1(:,1,:)*(basisP1(:,3) - basisP1(:,1)) + ...
                basisP1(:,3)*basisP1(:,1)*(gradP1(:,3,:) - gradP1(:,1,:));
            val(10,:) = gradP1(:,1,:)*basisP1(:,2)*basisP1(:,3) + ...
                basisP1(:,1)*gradP1(:,2,:)*basisP1(:,3) + ...
                basisP1(:,1)*basisP1(:,2)*gradP1(:,3,:);
        end

        %% P1 basis functions
        function val = basisP1(~,x,y)
            %basisP1 - Evaluate the P1 basis functions
            lambda1 = 1 - x - y;
            lambda2 = x;
            lambda3 = y;
            val = [lambda1,lambda2,lambda3];
        end

        function val = gradP1(~,x,y)
            %gradP1 - Evaluate gradients of the P1 basis
            val = zeros(size(x,1),3,2);
            val(:,1,:) = [-ones(size(x,1)),-ones(size(y,1))];
            val(:,2,:) = [ones(size(x,1)),zeros(size(y,1))];
            val(:,3,:) = [zeros(size(x,1)),ones(size(y,1))];
        end

        function val = intMoment0(~,points,refPoints,indices,u4Db,mesh)
            %intMoment0 - Integrand of the zeroth edge moment
            if isempty(points)
                n4edges = mesh.n4ed(indices,:);
                c1 = mesh.c4n(n4edges(:,1),:);
                c2 = mesh.c4n(n4edges(:,2),:);
                points = c1 + refPoints(1,1)*(c2-c1);
            end
            val = u4Db(points);
        end

        function val = intMoment1(~,points,refPoints,indices,u4Db,mesh)
            %intMoment1 - Integrand of the first edge moment
            if isempty(points)
                n4edges = mesh.n4ed(indices,:);
                c1 = mesh.c4n(n4edges(:,1),:);
                c2 = mesh.c4n(n4edges(:,2),:);
                points = c1 + refPoints(1,1)*(c2-c1);
            end
            tangents = mesh.tangent4ed;
            tangents = tangents(mesh.ed4Db,:);
            mid = mesh.point4ed(1/2);
            mid = mid(mesh.ed4Db,:);
            val = u4Db(points);
            val = (sum((points-mid).*tangents,2)).*val;
        end

        function val = intMoment2(~,points,refPoints,indices,u4Db,mesh)
            %intMoment2 - Integrand of the second edge moment
            if isempty(points)
                n4edges = mesh.n4ed(indices,:);
                c1 = mesh.c4n(n4edges(:,1),:);
                c2 = mesh.c4n(n4edges(:,2),:);
                points = c1 + refPoints(1,1)*(c2-c1);
            end
            tangents = mesh.tangent4ed;
            tangents = tangents(mesh.ed4Db,:);
            mid = mesh.point4ed(1/2);
            mid = mid(mesh.ed4Db,:);
            val = u4Db(points);
            val = (sum((points-mid).*tangents,2)).^2.*val;
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
            obj.nrBasis4e = double(obj.degree == 3);
        end

        function obj = dofManager(obj)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ) enumerates the edge and interior
            %   degrees of freedom with FE.dofManager and fixes them on
            %   the Dirichlet boundary.
            obj = dofManager@FE(obj);
            obj.fixedNodes = obj.dof4Db;
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
        end

    end

end