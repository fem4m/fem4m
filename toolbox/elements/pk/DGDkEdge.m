classdef DGDkEdge < FEEdge
    %DGDkEdge - Discontinuous polynomial edge space of degree k
    %   OBJ = DGDkEdge(DEGREE,MESH) creates the discontinuous
    %   polynomial edge space of degree DEGREE >= 0 on the triangular
    %   mesh MESH with the DEGREE+1 Legendre polynomials (Jacobi
    %   polynomials with alpha = beta = 0) of the edge parameter as
    %   basis functions per edge, replicated on both sides of each
    %   edge (zero on the second side of a boundary edge).
    %
    %   DGDkEdge functions:
    %       getEdgeBasis    - Evaluate the edge basis functions
    %       dofManager      - Assemble the degrees of freedom
    %       getEdgeRefBasis - Evaluate the edge reference basis
    %
    %   DGDkEdge properties:
    %       mesh         - Mesh on which the element is defined
    %       degree       - Polynomial degree of the element
    %       nrBasisFunc  - Number of basis functions per edge
    %       dof4ed       - Local-to-global dofs per edge
    %       nrDoF        - Total number of degrees of freedom
    %       freeNodes    - Indices of the free degrees of freedom
    %       fixedNodes   - Indices of the fixed degrees of freedom
    %
    %   See also FEEdge, DGDkFace, DGDk, DGDk3d

    % Copyright (c) 2026 by the fem4m authors

    properties
        mesh          % mesh on which the element is defined
        degree        % polynomial degree of the element
        nrBasisFunc   % number of basis functions per edge
        dof4ed        % local-to-global dofs per edge
        nrDoF         % total number of degrees of freedom
        freeNodes     % indices of the free degrees of freedom
        fixedNodes    % indices of the fixed degrees of freedom
    end

    methods

        function obj = DGDkEdge(degree,mesh)
            %DGDkEdge - Create a discontinuous polynomial edge space
            obj.degree = degree;
            obj.nrBasisFunc = obj.degree+1;
            obj.mesh = mesh;
            obj = obj.dofManager;
        end

        function val = getEdgeBasis(obj,~,x,parts)
            %getEdgeBasis - Evaluate the edge basis functions
            %   VAL = getEdgeBasis(OBJ,PTS,X,PARTS) evaluates the
            %   Legendre edge basis of degree OBJ.degree at the
            %   point X (a scalar in [0,1]) on the edges PARTS and
            %   replicates the result on both sides of each edge,
            %   the elements e4ed(PARTS,1) and e4ed(PARTS,2); the
            %   second side of a boundary edge is zero. VAL is
            %   numel(PARTS)-by-nrBasisFunc-by-2. PTS is ignored.
            val = obj.getEdgeRefBasis(obj.degree,x);
            val = ones(numel(parts),1,2).*val;
            val(obj.mesh.e4ed(parts,2)==0,:,2)=0;
        end

        function obj = dofManager(obj)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ) numbers the degrees of
            %   freedom index by index (the first basis function on
            %   all edges, then the second, ...), fixes the degrees
            %   of freedom on the Dirichlet edges ed4Db, and frees
            %   all others.
            % dofs
            obj.dof4ed = zeros(obj.mesh.nrEdges,obj.degree+1);
            for i = 1 : obj.degree+1
                obj.dof4ed(:,i) = (i-1)*obj.mesh.nrEdges+(1:obj.mesh.nrEdges);
            end
            % nrDoF
            obj.nrDoF = max(obj.dof4ed(:));
            obj.fixedNodes = zeros(numel(obj.mesh.ed4Db),obj.degree+1);
            for i = 1 : obj.degree+1
                obj.fixedNodes(:,i) =  (i-1)*obj.mesh.nrEdges+obj.mesh.ed4Db;
            end
            obj.fixedNodes  = obj.fixedNodes(:);
            obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
        end

    end

    methods (Static)
        %% Edge basis helper functions
        function val = getEdgeRefBasis(p,x)
            %getEdgeRefBasis - Evaluate the edge reference basis
            %   VAL = getEdgeRefBasis(P,X) evaluates the Legendre
            %   polynomials of degree 0 to P in 2*X-1 at the points
            %   X of the reference interval [0,1]. VAL is
            %   numel(X)-by-(P+1).
            val = JacobiPolynomials.poly(2*x-1,p,0,0,0)';
        end
    end

end