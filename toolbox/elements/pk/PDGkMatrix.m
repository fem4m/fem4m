classdef PDGkMatrix < FEMatrix
    %PDGkMatrix - Matrix-valued PDGk element built from PDGk elements
    %   OBJ = PDGkMatrix(DEGREE,MESH) stacks d*d copies of the
    %   scalar element of polynomial degree DEGREE >= 0, PDGk on the
    %   polygonal mesh MESH (MeshPoly, d = 2) or PDGk3d on the
    %   tetrahedral mesh MESH (Mesh3d, d = 3), into one matrix-valued
    %   element with the components ordered column by column, and
    %   manages their combined degrees of freedom.
    %
    %   See also PDGk, PDGk3d, PDGkVector, PDGkSymMatrix, FEMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = PDGkMatrix(degree,mesh)
            %PDGkMatrix - Create a matrix-valued PDGk finite element object
            if size(mesh.c4n,2) == 2
                pk = PDGk(degree,mesh);
            elseif size(mesh.c4n,2) == 3
                pk = PDGk3d(degree,mesh);
            end
            args = cell(size(mesh.c4n,2)*size(mesh.c4n,2),1);
            for k = 1 : numel(args)
                args{k} = pk;
            end
            obj = obj@FEMatrix(args{:});
        end

    end

end