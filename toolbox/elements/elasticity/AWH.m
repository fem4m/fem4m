classdef AWH < AWk
    %AWH - Hybridized Arnold-Winther element of degree k for elasticity
    %   OBJ = AWH(DEGREE,MESH) creates the Arnold-Winther element
    %   AWk(DEGREE,MESH) of polynomial degree DEGREE >= 3 with hybridized
    %   vertex constraints: at every non-parallel Neumann corner the
    %   second of the two triangles at the corner gets its own three
    %   nodal stress dofs, so the stress is double-valued there, and the
    %   six nodal values satisfy both traction conditions and the
    %   continuity of the normal stress across the interior edge. Every
    %   such corner needs exactly two triangles. AWH(3,MESH) is the
    %   hybridized space of AW1.
    %
    %   AWH functions:
    %       boundaryConstraints - Build the hybridized vertex constraints
    %
    %   See also AWk, AW1, FEMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = AWH(degree,mesh)
            %AWH - Create a hybridized Arnold-Winther element object
            if nargin < 2
                [degree,mesh] = deal([]);
            end
            obj = obj@AWk(degree,mesh);
        end

        function obj = boundaryConstraints(obj)
            %boundaryConstraints - Build the hybridized vertex constraints
            %   OBJ = boundaryConstraints(OBJ) switches the hybridization
            %   on and builds the vertex constraints of AW1, which add
            %   hybrid dofs at non-parallel Neumann corners.
            obj.hybrid = true;
            obj = boundaryConstraints@AWk(obj);
        end

    end
end