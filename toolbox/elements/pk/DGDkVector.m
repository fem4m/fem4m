classdef DGDkVector < FEVector
    %DGDkVector - Vector-valued element built from DGDk or DGDk3d elements
    %   OBJ = DGDkVector(DEGREE,MESH) stacks one scalar element of
    %   total degree DEGREE per spatial dimension, DGDk on the
    %   triangular mesh MESH (two components) or DGDk3d on the
    %   tetrahedral mesh MESH (three components), into one
    %   vector-valued element. The dimension is the number of
    %   columns of MESH.c4n.
    %
    %   See also DGDk, DGDk3d, FEVector, DGPkVector, DGDkMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = DGDkVector(degree,mesh)
            %DGDkVector - Create a vector-valued DGDk finite element object
            if size(mesh.c4n,2) == 2
                pk = DGDk(degree,mesh);
            elseif size(mesh.c4n,2) == 3
                pk = DGDk3d(degree,mesh);
            end
            args = cell(size(mesh.c4n,2),1);
            for k = 1 : size(mesh.c4n,2)
                args{k} = pk;
            end
            obj = obj@FEVector(args{:});
        end

    end

end