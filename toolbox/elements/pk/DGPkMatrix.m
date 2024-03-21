classdef DGPkMatrix < FEMatrix
    %DGPkMatrix - Matrix-valued DGPk element built from DGPk elements
    %   OBJ = DGPkMatrix(DEGREE,MESH) stacks d^2 copies of the
    %   scalar element of polynomial degree DEGREE >= 1, DGPk on the
    %   triangular mesh MESH (d = 2) or DGPk3d on the tetrahedral
    %   mesh MESH (d = 3), into one matrix-valued element with the
    %   components ordered column by column, and manages their
    %   combined degrees of freedom.
    %
    %   See also DGPk, DGPk3d, DGPkVector, DGPkSymMatrix, FEMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = DGPkMatrix(degree,mesh)
            %DGPkMatrix - Create a matrix-valued DGPk finite element object
            if size(mesh.c4n,2) == 2
                pk = DGPk(degree,mesh);
            elseif size(mesh.c4n,2) == 3
                pk = DGPk3d(degree,mesh);
            end
            args = cell(size(mesh.c4n,2)*size(mesh.c4n,2),1);
            for k = 1 : numel(args)
                args{k} = pk;
            end
            obj = obj@FEMatrix(args{:});
        end

    end

end