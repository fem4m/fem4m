classdef DGHPk < HPk & FE
    %DGHPk - hp-version discontinuous Pk element with elementwise degrees
    %   OBJ = DGHPk(DEGREE,MESH) creates the discontinuous hp-version
    %   Pk element on the triangular mesh MESH. DEGREE >= 1 is a scalar
    %   or an nrElems-by-1 vector with one polynomial degree per
    %   element. Each element carries the basis of HPk of the maximum
    %   degree without orientation factors (basisCoefficients are all
    %   ones); all degrees of freedom are element-local, and those above
    %   the degree of the element are unused (nonDoFNodes).
    %
    %   DGHPk inherits the functions and properties of HPk, with these
    %   differences: there are no Dirichlet dofs (dof4Db and fixedNodes
    %   are empty, freeNodes are all used dofs), dof4ed is empty, and
    %   degree4ed is the maximum of the degrees of the two elements at
    %   each edge. The inherited projectDbBoundary does not apply.
    %
    %   See also HPk, DGDk, DGPk, DGDk3d

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = DGHPk(degree,mesh)
            %DGHPk - Create an hp-version discontinuous Pk element
            obj@HPk(degree,mesh);
            obj.basisCoefficients = ones(obj.mesh.nrElems,obj.nrBasisFunc);
        end
    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) treats all nrBasisFunc degrees of
            %   freedom as element-local ones.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 0;
            obj.nrBasis4e = obj.nrBasisFunc;
        end

        function obj = dofManager(obj)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ) numbers the element-local
            %   degrees of freedom of dofCounts with FE.dofManager
            %   (bypassing HPk.dofManager), keeps only dof4e (the side
            %   tables of HPk stay empty), records the polynomial
            %   degree of each local dof, sets degree4ed to the maximum
            %   degree of the two elements at each edge, and collects in
            %   nonDoFNodes the dofs above degree4e of each element.
            %   nrDoF counts the used dofs, and freeNodes are all of
            %   them.
            obj = dofManager@FE(obj);
            obj.dof4ed = [];
            obj.dof4Db = [];
            obj.degree4DoF(1) = 0;
            %first 3 nodes
            if obj.degree > 0
                obj.degree4DoF(1:3) = 1;
            end
            % 3*(p-1) for edges E1;E2;E3
            if obj.degree > 1
                obj.degree4DoF(4:3+3*(obj.degree-1)) = [2:obj.degree,2:obj.degree,2:obj.degree];
            end
            % (p-2)*(p-1)/2 interior, sorted by degree (HPk.elementDerivatives)
            if obj.degree > 2
                for p = 3 : obj.degree
                    obj.degree4DoF(:,3+3*(obj.degree-1)+(1+(p-3)*(p-2)/2:(p-2)*(p-1)/2)) = p;
                end
            end

            obj.allnrDoF = obj.nrDoF;
            obj.fixedNodes = [];
            e4ed = obj.mesh.e4ed;
            e4ed(e4ed(:,2)==0,2) = e4ed(e4ed(:,2)==0,1);
            obj.degree4ed = max(obj.degree4e(e4ed),[],2);
            nonDoFNodes = [];
            for curDegree = 1 : obj.degree
                %inner dofs
                nonDoFNodes = [nonDoFNodes;reshape(obj.dof4e(obj.degree4e==curDegree,...
                    3+3*(obj.degree-1)+(curDegree-2)*(curDegree-1)/2+1:end),[],1)]; %#ok<AGROW>
                % E1
                nonDoFNodes = [nonDoFNodes;reshape(obj.dof4e(obj.degree4e==curDegree,...
                    3+(curDegree-1)+1:3+(obj.degree-1)),[],1)]; %#ok<AGROW>
                % E2
                nonDoFNodes = [nonDoFNodes;reshape(obj.dof4e(obj.degree4e==curDegree,...
                    3+(obj.degree-1)+(curDegree-1)+1:3+2*(obj.degree-1)),[],1)]; %#ok<AGROW>
                % E3
                nonDoFNodes = [nonDoFNodes;reshape(obj.dof4e(obj.degree4e==curDegree,...
                    3+2*(obj.degree-1)+(curDegree-1)+1:3+3*(obj.degree-1)),[],1)]; %#ok<AGROW>
            end
            nonDoFNodes = unique(nonDoFNodes);
            obj.nonDoFNodes = nonDoFNodes;
            obj.nrDoF = obj.allnrDoF - numel(nonDoFNodes);
            obj.freeNodes = setdiff(1:obj.allnrDoF,nonDoFNodes);
        end

    end
end