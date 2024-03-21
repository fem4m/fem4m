classdef PDGkSymMatrix < FESymMatrix
    %PDGkSymMatrix - Symmetric matrix-valued PDGk element
    %   OBJ = PDGkSymMatrix(DEGREE,MESH) stacks d*(d+1)/2 copies
    %   of the scalar element of polynomial degree DEGREE >= 0, PDGk on
    %   the polygonal mesh MESH (MeshPoly, d = 2) or PDGk3d on the
    %   tetrahedral mesh MESH (Mesh3d, d = 3), into one symmetric
    %   matrix-valued element with the components [11 22 12] or
    %   [11 22 33 12 13 23], and manages their combined degrees of
    %   freedom.
    %
    %   See also PDGk, PDGk3d, PDGkMatrix, PDGkVector, FESymMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = PDGkSymMatrix(degree,mesh)
            %PDGkSymMatrix - Create a symmetric matrix-valued element object
            if size(mesh.c4n,2) == 2
                pk = PDGk(degree,mesh);
            elseif size(mesh.c4n,2) == 3
                pk = PDGk3d(degree,mesh);
            end
            args = cell(size(mesh.c4n,2)*(size(mesh.c4n,2)+1)/2,1);
            for k = 1 : numel(args)
                args{k} = pk;
            end
            obj = obj@FESymMatrix(args{:});
        end

    end

end