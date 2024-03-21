classdef FEEdgeVector < FEEdge
    %FEEdgeVector - Vector-valued edge element built from scalar elements
    %   OBJ = FEEdgeVector(FE1) creates a vector-valued edge element from
    %   the scalar edge element FE1.
    %
    %   OBJ = FEEdgeVector(FE1,FE2,...) stacks several scalar edge
    %   elements, e.g. DGDkEdge, into one vector-valued element on the
    %   mesh of FE1 and numbers their degrees of freedom component after
    %   component.
    %
    %   FEEdgeVector functions:
    %       dofManager   - Assemble the combined degrees of freedom
    %       getEdgeBasis - Evaluate the edge basis functions
    %
    %   FEEdgeVector properties:
    %       mesh        - Mesh on which the elements are defined
    %       fe          - Cell array of the scalar edge elements
    %       dof4ed      - Local-to-global degrees of freedom per edge
    %       nrDoF       - Total number of degrees of freedom
    %       freeNodes   - Indices of the free degrees of freedom
    %       fixedNodes  - Indices of the fixed degrees of freedom
    %       nrBasisFunc - Number of basis functions per edge
    %       degree      - Largest polynomial degree of the components
    %
    %   See also FEEdge, FEVector, DGDkEdge

    % Copyright (c) 2026 by the fem4m authors

    properties
        mesh          % mesh on which the elements are defined
        fe            % cell array of the scalar edge elements
        dof4ed        % local-to-global degrees of freedom per edge
        nrDoF         % total number of degrees of freedom
        freeNodes     % indices of the free degrees of freedom
        fixedNodes    % indices of the fixed degrees of freedom
        nrBasisFunc   % number of basis functions per edge
        degree        % largest polynomial degree of the components
    end

    methods
        function obj = FEEdgeVector(varargin)
            %FEEdgeVector - Create a vector-valued edge element object
            obj.fe = varargin(:);
            obj.mesh = obj.fe{1}.mesh;
            obj.degree = obj.fe{1}.degree;
            for curFE = 2 : numel(obj.fe)
                obj.degree = max(obj.degree,obj.fe{curFE}.degree);
            end
            obj = obj.dofManager;
        end

        function obj = dofManager(obj)
            %dofManager - Assemble the combined degrees of freedom
            %   OBJ = dofManager(OBJ) concatenates the local-to-global maps,
            %   the free and fixed degrees of freedom, and the basis counts
            %   of all scalar elements in OBJ.fe with shifted indices.
            obj.dof4ed = obj.fe{1}.dof4ed;
            obj.nrDoF = obj.fe{1}.nrDoF;
            obj.freeNodes = obj.fe{1}.freeNodes;
            obj.fixedNodes = obj.fe{1}.fixedNodes;
            obj.nrBasisFunc = obj.fe{1}.nrBasisFunc;
            for curFE = 2 : numel(obj.fe)
                obj.dof4ed = [obj.dof4ed,obj.nrDoF+obj.fe{curFE}.dof4ed];
                obj.freeNodes = [obj.freeNodes(:);obj.nrDoF+obj.fe{curFE}.freeNodes(:)];
                obj.fixedNodes = [obj.fixedNodes(:);obj.nrDoF+obj.fe{curFE}.fixedNodes(:)];
                obj.nrDoF = obj.nrDoF + obj.fe{curFE}.nrDoF;
                obj.nrBasisFunc = obj.nrBasisFunc+obj.fe{curFE}.nrBasisFunc;
            end
        end

        function val = getEdgeBasis(obj,pts,pts_ref,parts)
            %getEdgeBasis - Evaluate the edge basis functions
            %   VAL = getEdgeBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the edge
            %   basis functions of all scalar elements in OBJ.fe at the
            %   point PTS_REF (1-by-1) of the reference edge [0,1] on the
            %   edges PARTS. VAL(:,J,K,S) is component K of the J-th basis
            %   function on side S of the edges: the basis functions of
            %   OBJ.fe{K} form the K-th block of columns and vanish in the
            %   other components.
            val = permute(obj.fe{1}.getEdgeBasis(pts,pts_ref,parts),[1 2 4 3]);
            for curFE = 2 : size(obj.fe,1)
                val(:,end+(1:obj.fe{curFE}.nrBasisFunc),curFE,:) = obj.fe{curFE}.getEdgeBasis(pts,pts_ref,parts);
            end
        end

    end
end