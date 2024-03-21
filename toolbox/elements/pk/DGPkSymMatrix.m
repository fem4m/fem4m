classdef DGPkSymMatrix < FESymMatrix
    %DGPkSymMatrix - Symmetric matrix-valued DGPk element
    %   OBJ = DGPkSymMatrix(DEGREE,MESH) stacks d*(d+1)/2 copies
    %   of the scalar element of polynomial degree DEGREE >= 1, DGPk on
    %   the triangular mesh MESH (d = 2) or DGPk3d on the
    %   tetrahedral mesh MESH (d = 3), into one symmetric
    %   matrix-valued element with the components [11 22 12] or
    %   [11 22 33 12 13 23], and manages their combined degrees of
    %   freedom.
    %
    %   See also DGPk, DGPk3d, DGPkMatrix, FESymMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = DGPkSymMatrix(degree,mesh)
            %DGPkSymMatrix - Create a symmetric matrix-valued element object
            if size(mesh.c4n,2) == 2
                pk = DGPk(degree,mesh);
            elseif size(mesh.c4n,2) == 3
                pk = DGPk3d(degree,mesh);
            end
            args = cell(size(mesh.c4n,2)*(size(mesh.c4n,2)+1)/2,1);
            for k = 1 : numel(args)
                args{k} = pk;
            end
            obj = obj@FESymMatrix(args{:});
        end

    end

end