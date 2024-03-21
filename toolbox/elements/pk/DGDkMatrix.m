classdef DGDkMatrix < FEMatrix
    %DGDkMatrix - Matrix-valued element built from DGDk or DGDk3d elements
    %   OBJ = DGDkMatrix(DEGREE,MESH) stacks d^2 copies of the
    %   scalar element of total degree DEGREE, DGDk on the
    %   triangular mesh MESH (d = 2) or DGDk3d on the tetrahedral
    %   mesh MESH (d = 3), into one matrix-valued element with the
    %   components ordered column by column, and manages their
    %   combined degrees of freedom.
    %
    %   See also DGDk, DGDk3d, DGDkVector, DGPkMatrix, PDGkMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = DGDkMatrix(degree,mesh)
            %DGDkMatrix - Create a matrix-valued DGDk finite element object
            if size(mesh.c4n,2) == 2
                pk = DGDk(degree,mesh);
            elseif size(mesh.c4n,2) == 3
                pk = DGDk3d(degree,mesh);
            end
            args = cell(size(mesh.c4n,2)*size(mesh.c4n,2),1);
            for k = 1 : numel(args)
                args{k} = pk;
            end
            obj = obj@FEMatrix(args{:});
        end

    end

end