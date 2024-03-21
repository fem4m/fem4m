classdef Nd1DGk3d < Nd1k3d
    %Nd1DGk3d - Discontinuous Nedelec first family element of degree k on tetrahedra
    %   OBJ = Nd1DGk3d(DEGREE,MESH) creates the discontinuous variant of
    %   the tetrahedral Nedelec first family element Nd1k3d of polynomial
    %   degree DEGREE >= 1 on the tetrahedral mesh MESH. It uses the
    %   basis functions of Nd1k3d, but all edge, face and cell degrees
    %   of freedom are element local, so no tangential continuity is
    %   enforced and the space is not H(curl)-conforming. Only dof4e is
    %   set (dof4ed, dof4f and dof4Db are empty), and all degrees of
    %   freedom are free; curlfreeNodes keeps the Nd1k3d numbering.
    %
    %   See also Nd1k3d, Nd1k, Nd2k3d

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = Nd1DGk3d(degree,mesh)
            %Nd1DGk3d - Create a discontinuous Nedelec element object
            obj = obj@Nd1k3d(degree,mesh);
            obj = obj.dofManager();
        end
    end

    methods (Access = protected)

        function obj = dofManager(obj,mesh)
            %dofManager - Manage the degrees of freedom
            %   OBJ = dofManager(OBJ) enumerates purely element-local
            %   degrees of freedom for the edge, face and cell basis
            %   functions without imposing any continuity between elements
            %   and keeps only dof4e: the tables dof4ed, dof4f and dof4Db of
            %   the Nd1k3d numbering are cleared. No degrees of freedom are
            %   fixed: fixedNodes is empty and freeNodes lists all degrees of
            %   freedom. curlfreeNodes is not updated and still holds the
            %   Nd1k3d numbers.
            %
            %   OBJ = dofManager(OBJ,MESH), as called by the Nd1k3d
            %   constructor, assembles the Nd1k3d degrees of freedom.
            if nargin > 1
                obj = dofManager@Nd1k3d(obj,mesh);
                return
            end
            nrElems = obj.mesh.nrElems;
            % nrDoF
            obj.nrDoF = 6*nrElems*obj.nrBasis4ed + 4*nrElems*obj.nrBasis4f + nrElems*obj.nrBasis4e;
            % elems
            obj.dof4e = reshape(1:obj.nrDoF,nrElems,obj.nrBasisFunc);
            obj.dof4ed = [];
            obj.dof4f = [];
            obj.dof4Db = [];
            % fixed/freeNodes
            obj.fixedNodes = [];
            obj.freeNodes = 1:obj.nrDoF;
        end

    end
end