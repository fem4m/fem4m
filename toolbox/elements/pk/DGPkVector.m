classdef DGPkVector < FEVector
    %DGPkVector - Vector-valued DGPk element built from DGPk elements
    %   OBJ = DGPkVector(DEGREE,MESH) stacks d copies of the scalar
    %   element of polynomial degree DEGREE >= 1, DGPk on the
    %   triangular mesh MESH (d = 2) or DGPk3d on the tetrahedral
    %   mesh MESH (d = 3), into one vector-valued element.
    %
    %   See also DGPk, DGPk3d, PkVector, DGDkVector

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = DGPkVector(degree,mesh)
            %DGPkVector - Create a vector-valued DGPk finite element object
            if size(mesh.c4n,2) == 2
                pk = DGPk(degree,mesh);
            elseif size(mesh.c4n,2) == 3
                pk = DGPk3d(degree,mesh);
            end
            args = cell(size(mesh.c4n,2),1);
            for k = 1 : numel(args)
                args{k} = pk;
            end
            obj = obj@FEVector(args{:});
        end

    end

end