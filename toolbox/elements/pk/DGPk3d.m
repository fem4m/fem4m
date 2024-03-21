classdef DGPk3d < Pk3d
    %DGPk3d - Discontinuous Lagrange Pk element of degree k on tetrahedra
    %   OBJ = DGPk3d(DEGREE,MESH) creates the discontinuous
    %   scalar Lagrange finite element of polynomial degree
    %   DEGREE >= 1 on the tetrahedral mesh MESH. It has the basis
    %   functions of Pk3d, but each of its Lagrange degrees of
    %   freedom, also those at the vertices, on the edges, and on
    %   the faces, belongs to one element only; the degrees of
    %   freedom are numbered element by element. No degree of
    %   freedom is fixed: freeNodes holds all of them.
    %
    %   See also Pk3d, DGPk, DGPkVector, DGDk3d

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = DGPk3d(degree,mesh)
            %DGPk3d - Create a discontinuous Pk3d finite element object
            obj = obj@Pk3d(degree,mesh);
        end
    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) treats all degrees of freedom as
            %   element-local ones.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 0;
            obj.nrBasis4f = 0;
            obj.nrBasis4e = 4 + 6*(obj.degree-1) + 4*(obj.degree-1)*(obj.degree-2)/2 + (obj.degree-1)*(obj.degree-2)*(obj.degree-3)/6;
        end

        function obj = dofManager(obj,mesh)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ,MESH) numbers the element-local
            %   degrees of freedom of dofCounts with Pk3d.dofManager,
            %   keeps only dof4e (it empties dof4ed, dof4f, and dof4Db)
            %   and computes the coordinates of the dof points.
            obj = dofManager@Pk3d(obj,mesh);
            obj.dof4ed = [];
            obj.dof4f = [];
            obj.dof4Db = [];

            %% c4dof
            obj.c4dof = zeros(obj.nrDoF,3);
            % nodal
            for k = 1 : 4
                obj.c4dof(obj.dof4e(:,k),:) = mesh.c4n(mesh.n4e(:,k),:);
            end
            % edge
            for k = 1 : obj.nrEdgePoints
                for j = 1 : 6
                    obj.c4dof(obj.dof4e(:,4+(k-1)*6+j),:) = ...
                        obj.local2global(obj.mesh.n4e,obj.c4refDof(4+(k-1)*6+j,:),mesh);
                end
            end
            % face points
            for k = 1 : obj.nrFacePoints
                for j = 1 : 4
                    obj.c4dof(obj.dof4e(:,4+6*obj.nrEdgePoints+(k-1)*4+j),:) = ...
                        obj.local2global(obj.mesh.n4e,obj.c4refDof(4+6*obj.nrEdgePoints+(k-1)*4+j,:),mesh);
                end
            end

            % inner Points
            for k = 1 : obj.nrInnerPoints
                obj.c4dof(obj.dof4e(:,4+6*obj.nrEdgePoints+4*obj.nrFacePoints+k),:) = ...
                    obj.local2global(obj.mesh.n4e,obj.c4refDof(4+6*obj.nrEdgePoints+4*obj.nrFacePoints+k,:),mesh);
            end
        end

    end
end