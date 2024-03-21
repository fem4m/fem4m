classdef PDGkVector < FEVector
    %PDGkVector - Vector-valued PDGk element built from PDGk elements
    %   OBJ = PDGkVector(DEGREE,MESH) stacks one scalar element of
    %   polynomial degree DEGREE >= 0 per spatial dimension of the
    %   mesh MESH, PDGk on the polygonal mesh MESH (MeshPoly) or PDGk3d on
    %   the tetrahedral mesh MESH (Mesh3d), into one vector-valued
    %   element.
    %
    %   PDGkVector functions:
    %       plotDisplacement - Plot the displaced polygonal mesh
    %
    %   See also PDGk, PDGk3d, FEVector, PDGkMatrix, DGDkVector

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = PDGkVector(degree,mesh)
            %PDGkVector - Create a vector-valued PDGk finite element object
            if size(mesh.c4n,2) == 2
                pk = PDGk(degree,mesh);
            elseif size(mesh.c4n,2) == 3
                pk = PDGk3d(degree,mesh);
            end
            args = cell(size(mesh.c4n,2),1);
            for k = 1 : size(mesh.c4n,2)
                args{k} = pk;
            end
            obj = obj@FEVector(args{:});
        end

        function plotDisplacement(obj,c,scale,c2,fe2)
            %plotDisplacement - Plot the displaced polygonal mesh
            %   plotDisplacement(OBJ,C) plots the 2d polygonal mesh
            %   displaced by the vector-valued finite element
            %   function with coefficients C.
            %
            %   plotDisplacement(OBJ,C,SCALE) amplifies the
            %   displacement by the factor SCALE (default 1, also for
            %   an empty SCALE).
            %
            %   plotDisplacement(OBJ,C,SCALE,C2,FE2) also colors
            %   the displaced mesh by the finite element function
            %   with coefficients C2 on the element FE2 at the
            %   vertices: by its values for scalar data, by the von
            %   Mises stress for 2d matrix data (3 or 4 components),
            %   and by the norm otherwise.
            if nargin < 3 || isempty(scale)
                scale = 1;
            end
            c4n = obj.mesh.c4n;
            n4e = obj.mesh.n4e;
            for j = unique(obj.mesh.nrN4e)'
                elems = find(obj.mesh.nrN4e==j);
                curN4e = cell2mat(n4e(obj.mesh.nrN4e==j));
                u = [];
                coord = [];
                for k = 1 : j
                    u = [u,obj.getValue(c4n(curN4e(:,k),:),[],elems,c)]; %#ok<AGROW>
                    coord = [coord,c4n(curN4e(:,k),:)]; %#ok<AGROW>
                end
                u = scale*u;
                if nargin < 5
                    patch(coord(:,1:2:end)'+u(:,1:2:end)',coord(:,2:2:end)'+u(:,2:2:end)','w');
                else
                    color = zeros(numel(elems),j);
                    for k = 1 : j
                        value = fe2.getValue(c4n(curN4e(:,k),:),[],elems,c2);
                        if size(value,2) == 1
                            % value
                            color(:,k) = value;
                        elseif size(value,2) == 3
                            % von Mises stress
                            color(:,k) = sqrt(value(:,1).^2+value(:,2).^2-value(:,1).*value(:,2)+3*value(:,3).^2);
                        elseif size(value,2) == 4
                            % von Mises stress
                            color(:,k) = sqrt(value(:,1).^2+value(:,4).^2-value(:,1).*value(:,4)+3*(1/2*(value(:,2)+value(:,3))).^2);
                        else
                            % norm
                            color(:,k) = sqrt(sum(value.^2,[2 3]));
                        end
                    end
                    patch(coord(:,1:2:end)'+u(:,1:2:end)',coord(:,2:2:end)'+u(:,2:2:end)',color');
                end
                hold on;
            end
        end

    end

end