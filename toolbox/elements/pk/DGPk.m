classdef DGPk < Pk
    %DGPk - Discontinuous Lagrange Pk element of degree k
    %   OBJ = DGPk(DEGREE,MESH) creates the discontinuous scalar
    %   Lagrange finite element of polynomial degree DEGREE >= 1 on
    %   the triangular mesh MESH. It has the basis functions of Pk,
    %   but each of its Lagrange degrees of freedom, also those at
    %   the vertices and on the edges, belongs to one element only;
    %   the degrees of freedom are numbered element by element. No
    %   degree of freedom is fixed: freeNodes holds all of them.
    %
    %   DGPk functions:
    %       local2global - Map reference dofs to global points
    %
    %   See also Pk, DGPk3d, DGPkVector, DGDk

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = DGPk(degree,mesh)
            %DGPk - Create a discontinuous Pk finite element object
            obj@Pk(degree,mesh);
        end

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
    end

    methods (Access = protected)

        function obj = dofCounts(obj)
            %dofCounts - Set the numbers of dofs per entity
            %   OBJ = dofCounts(OBJ) treats all degrees of freedom as
            %   element-local ones.
            obj.nrBasis4n = 0;
            obj.nrBasis4ed = 0;
            obj.nrBasis4e = obj.nrBasisFunc;
        end

        function obj = dofManager(obj)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ) numbers the degrees of freedom
            %   with Pk.dofManager for the counts of dofCounts (all
            %   element-local), keeps only dof4e (it empties dof4ed and
            %   dof4Db) and sets the coordinates of the dof points.
            obj = dofManager@Pk(obj);
            obj.dof4ed = [];
            obj.dof4Db = [];

            n4e = obj.mesh.n4e;
            c4n = obj.mesh.c4n;
            points = obj.refPoints;

            %% c4dof
            obj.c4dof = zeros(obj.nrDoF,size(c4n,2));
            % node points
            for j = 1 : 3
                obj.c4dof(obj.dof4e(:,j),:) = c4n(n4e(:,j),:);
            end
            % edge points
            for j = 1 : (obj.degree-1)
                obj.c4dof(obj.dof4e(:,3+3*(j-1)+1),:) = obj.local2global(n4e,[points(j),0]);
                obj.c4dof(obj.dof4e(:,3+3*(j-1)+2),:) = obj.local2global(n4e,[1-points(j),points(j)]);
                obj.c4dof(obj.dof4e(:,3+3*(j-1)+3),:) = obj.local2global(n4e,[0,1-points(j)]);
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