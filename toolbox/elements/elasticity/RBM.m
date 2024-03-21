classdef RBM < FE
    %RBM - Rigid body motion element for elasticity
    %   OBJ = RBM(DEGREE,MESH) creates the piecewise rigid body motion
    %   finite element on the triangular mesh MESH with three degrees of
    %   freedom per element: the two translations (1,0) and (0,1) and
    %   the infinitesimal rotation (-(y-y_T),x-x_T)/sqrt(|T|) about the
    %   centroid (x_T,y_T) of the element T. DEGREE must be 1. The basis
    %   depends on the element geometry (isAffineEquivalent is false).
    %
    %   RBM functions:
    %       getBasis            - Evaluate the rigid body basis functions
    %       getGradBasis        - Evaluate the gradients of the basis functions
    %       getSymGradBasis     - Evaluate the symmetric gradients
    %       getD2Basis          - Evaluate second derivatives of the basis
    %       plotDisplacement    - Plot the displaced mesh
    %
    %   RBM properties:
    %       mesh           - Mesh on which the element is defined
    %       freeNodes      - Indices of the free degrees of freedom
    %       fixedNodes     - Indices of the fixed degrees of freedom
    %       degree         - Polynomial degree of the element
    %       nrBasisFunc    - Number of basis functions per element
    %       dof4e          - Local-to-global dofs per element
    %       nrDoF          - Total number of degrees of freedom
    %       nrBasis4n      - Number of nodal basis functions
    %       nrBasis4ed     - Number of edge basis functions
    %       nrBasis4e      - Number of interior basis functions
    %       isRefElemBased - Basis lives on the reference element
    %
    %   See also FE, Pk, AW0, AWNC0

    % Copyright (c) 2026 by the fem4m authors

    properties
        mesh             % mesh on which the element is defined
        freeNodes        % indices of the free degrees of freedom
        fixedNodes       % indices of the fixed degrees of freedom
        degree           % polynomial degree of the element
        nrBasisFunc      % number of basis functions per element
        dof4e            % local-to-global dofs per element
        nrDoF            % total number of degrees of freedom
        nrBasis4n        % number of nodal basis functions
        nrBasis4ed       % number of edge basis functions
        nrBasis4e        % number of interior basis functions
        isRefElemBased = true % basis lives on the reference element
    end

    methods

        function obj = RBM(degree,mesh)
            %RBM - Create a rigid body motion finite element object
            if nargin < 2 || ~isequal(degree,1)
                error('fem4m:RBM:invalidDegree', ...
                    'RBM degree must be 1, rigid body motions are affine');
            end
            obj.isAffineEquivalent = false;
            obj.degree = degree;
            obj.nrBasisFunc = 3;
            obj.mesh = mesh;
            obj = obj.dofCounts();
            obj = obj.dofManager;
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the rigid body basis functions
            %   VAL = getBasis(OBJ,PTS,X,PARTS) evaluates the piecewise
            %   rigid body basis functions of two translations and one
            %   scaled rotation about the element centroids in all
            %   elements at the reference point X (1-by-2), an
            %   nrElems-by-3-by-2 array. PTS and PARTS are ignored.
            d = (obj.mesh.point4e(x) - obj.mesh.point4e([1 1]/3))./sqrt(obj.mesh.area4e);
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,2);
            val(:,1,1) = 1;
            val(:,2,2) = 1;
            val(:,3,1) = -d(:,2);
            val(:,3,2) = d(:,1);
        end

        function val = getGradBasis(obj,~,~,~)
            %getGradBasis - Evaluate the gradients of the basis functions
            %   VAL = getGradBasis(OBJ,PTS,X,PARTS) returns the constant
            %   gradients of the rigid body basis functions as an
            %   nrElems-by-3-by-4 array, the entries J(i,j) = d_j u_i
            %   column by column [d_1u_1 d_1u_2 d_2u_1 d_2u_2]. PTS, X
            %   and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,2,2);
            val(:,3,1,2) = -1./sqrt(obj.mesh.area4e);
            val(:,3,2,1) = 1./sqrt(obj.mesh.area4e);
            % flatten the Jacobian (component,derivative) column by column
            val = val(:,:,:);
        end

        function val = getSymGradBasis(obj,~,~,~)
            %getSymGradBasis - Evaluate the symmetric gradients
            %   VAL = getSymGradBasis(OBJ,PTS,X,PARTS) returns zeros
            %   (nrElems-by-3-by-3) since the symmetric gradients
            %   (e_11,e_22,e_12) of rigid body motions vanish. PTS, X and
            %   PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,3);
        end

        function val = getD2Basis(obj,~,~,~)
            %getD2Basis - Evaluate second derivatives of the basis
            %   VAL = getD2Basis(OBJ,PTS,X,PARTS) returns zeros
            %   (nrElems-by-3-by-3) since the rigid body basis functions
            %   are affine. PTS, X and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,3);
        end

        function plotDisplacement(obj,c,scale)
            %plotDisplacement - Plot the displaced mesh
            %   plotDisplacement(OBJ,C,SCALE) plots, as white triangles in
            %   the current axes, the elements with their vertices
            %   displaced by the piecewise rigid body motion with the
            %   coefficients C (nrDoF-by-1), amplified by the factor
            %   SCALE.
            %
            %   plotDisplacement(OBJ,C) uses SCALE = 1.
            if nargin < 3 || isempty(scale)
                scale = 1;
            end
            c4n = obj.mesh.c4n;
            n4e = obj.mesh.n4e;
            u = scale*[obj.getValue([],[0 0],[],c),...
                obj.getValue([],[1 0],[],c),...
                obj.getValue([],[0 1],[],c)];
            coord = [c4n(n4e(:,1),:),c4n(n4e(:,2),:),c4n(n4e(:,3),:)];
            patch(coord(:,[1 3 5])'+u(:,[1 3 5])',coord(:,[2 4 6])'+u(:,[2 4 6])','w')
        end
    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) sets the numbers nrBasis4n,
            %   nrBasis4ed, and nrBasis4e of nodal, edge, and interior
            %   degrees of freedom. The constructor calls it right before
            %   dofManager, so that subclasses may override it.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 0;
            obj.nrBasis4e = obj.nrBasisFunc;
        end

        function obj = dofManager(obj)
            %dofManager - Assemble the element degrees of freedom
            %   OBJ = dofManager(OBJ) enumerates three element-local
            %   rigid body degrees of freedom per triangle with
            %   FE.dofManager and keeps all of them free.
            obj = dofManager@FE(obj);
            obj.fixedNodes = [];
            obj.freeNodes = 1:obj.nrDoF;
        end

    end
end