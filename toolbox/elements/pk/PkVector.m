classdef PkVector < FEVector
    %PkVector - Vector-valued Pk element built from scalar Pk elements
    %   OBJ = PkVector(DEGREE,MESH) stacks d copies of the scalar
    %   element of polynomial degree DEGREE >= 1, Pk on the
    %   triangular mesh MESH (d = 2) or Pk3d on the tetrahedral mesh
    %   MESH (d = 3), into one vector-valued element and manages
    %   their combined degrees of freedom. On a mesh with a symmetry
    %   boundary n4Sb, component k is also fixed (fixedNodes) on the
    %   symmetry sides whose last column of n4Sb is k; sides with the
    %   axis 0 (a general normal) stay free and, in 2d, are handled
    %   by curvedSymmetricBoundaryConstraints.
    %
    %   PkVector functions:
    %       curvedSymmetricBoundaryConstraints - Build boundary constraints
    %
    %   PkVector properties:
    %       dof4n       - Nodal dof numbers per component
    %       dof4ed      - Dof numbers per edge and component
    %       dof4f       - Dof numbers per face and component
    %       nrBasis4n   - Number of basis functions per vertex
    %       nrBasis4ed  - Number of basis functions per edge
    %       nrBasis4f   - Number of basis functions per face
    %       nrBasis4e   - Number of interior basis functions
    %       idof4ed     - Columns of dof4ed of the edge-interior dofs
    %       idof4f      - Columns of dof4f of the face-interior dofs
    %       idof4e      - Columns of dof4e of the interior dofs
    %
    %   See also FEVector, Pk, Pk3d

    % Copyright (c) 2026 by the fem4m authors

    properties
        dof4n        % nodal dof numbers per component
        dof4ed       % dof numbers per edge and component
        dof4f        % dof numbers per face and component
        nrBasis4n    % number of basis functions per vertex
        nrBasis4ed   % number of basis functions per edge
        nrBasis4f    % number of basis functions per face
        nrBasis4e    % number of interior basis functions
        idof4ed      % columns of dof4ed of the edge-interior dofs
        idof4f       % columns of dof4f of the face-interior dofs
        idof4e       % columns of dof4e of the interior dofs
    end

    methods
        function obj = PkVector(degree,mesh)
            %PkVector - Create a vector-valued Pk finite element object
            if size(mesh.c4n,2) == 2
                pk = Pk(degree,mesh);
            elseif size(mesh.c4n,2) == 3
                pk = Pk3d(degree,mesh);
            end
            args = cell(size(mesh.c4n,2),1);
            for k = 1 : size(mesh.c4n,2)
                args{k} = pk;
            end
            obj = obj@FEVector(args{:});
            % axis-parallel symmetry boundary: component k is fixed on
            % the sides with axis k (column dim+1 of n4Sb)
            if isprop(mesh,'n4Sb') && ~isempty(mesh.n4Sb)
                dim = size(mesh.c4n,2);
                axis4Sb = mesh.n4Sb(:,dim+1);
                fixedNodes4Sb = cell(numel(obj.fe),1);
                nrDoF = 0;
                for k = 1 : numel(obj.fe)
                    if size(mesh.c4n,2) == 2
                        dof4s = obj.fe{k}.dof4ed(mesh.ed4Sb(axis4Sb==k),:);
                    else
                        dof4s = obj.fe{k}.dof4f(mesh.f4Sb(axis4Sb==k),:);
                    end
                    fixedNodes4Sb{k} = nrDoF + unique(dof4s(:));
                    nrDoF = nrDoF + obj.fe{k}.nrDoF;
                end
                obj.fixedNodes = [obj.fixedNodes(:);vertcat(fixedNodes4Sb{:})];
                obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes);
            end
            % dof4n
            obj.nrBasis4n = numel(obj.fe);
            obj.dof4n = (1:obj.fe{1}.mesh.nrNodes)';
            nrDoF = obj.fe{1}.nrDoF;
            for curFE = 2 : numel(obj.fe)
                obj.dof4n = [obj.dof4n,nrDoF+(1:obj.fe{curFE}.mesh.nrNodes)'];
                nrDoF = nrDoF + obj.fe{curFE}.nrDoF;
            end
            % dof4ed
            obj.dof4ed = obj.fe{1}.dof4ed;
            obj.nrBasis4ed = obj.fe{1}.nrBasis4ed;
            obj.idof4ed = size(obj.fe{1}.dof4ed,2)-obj.fe{1}.nrBasis4ed+1:size(obj.fe{1}.dof4ed,2);
            nrDoF = obj.fe{1}.nrDoF;
            nrBasisFunc = size(obj.fe{1}.dof4ed,2);
            for curFE = 2 : numel(obj.fe)
                obj.dof4ed = [obj.dof4ed,nrDoF+obj.fe{curFE}.dof4ed];
                obj.nrBasis4ed = obj.nrBasis4ed + obj.fe{curFE}.nrBasis4ed;
                obj.idof4ed = [obj.idof4ed,nrBasisFunc+(size(obj.fe{curFE}.dof4ed,2)-obj.fe{curFE}.nrBasis4ed+1:size(obj.fe{1}.dof4ed,2))];
                nrDoF = nrDoF + obj.fe{curFE}.nrDoF;
                nrBasisFunc = nrBasisFunc + size(obj.fe{curFE}.dof4ed,2);
            end
            % dof4f
            if size(mesh.c4n,2) > 2
                obj.dof4f = obj.fe{1}.dof4f;
                obj.nrBasis4f = obj.fe{1}.nrBasis4f;
                obj.idof4f = size(obj.fe{1}.dof4f,2)-obj.fe{1}.nrBasis4f+1:size(obj.fe{1}.dof4f,2);
                nrDoF = obj.fe{1}.nrDoF;
                nrBasisFunc = size(obj.fe{1}.dof4f,2);
                for curFE = 2 : numel(obj.fe)
                    obj.dof4f = [obj.dof4f,nrDoF+obj.fe{curFE}.dof4f];
                    obj.nrBasis4f = obj.nrBasis4f + obj.fe{curFE}.nrBasis4f;
                    obj.idof4f = [obj.idof4f,nrBasisFunc+(size(obj.fe{curFE}.dof4f,2)-obj.fe{curFE}.nrBasis4f+1:size(obj.fe{1}.dof4f,2))];
                    nrDoF = nrDoF + obj.fe{curFE}.nrDoF;
                    nrBasisFunc = nrBasisFunc + size(obj.fe{curFE}.dof4f,2);
                end
            end
            obj.nrBasis4e = obj.fe{1}.nrBasis4e;
            obj.idof4e = size(obj.fe{1}.dof4e,2)-obj.fe{1}.nrBasis4e+1:size(obj.fe{1}.dof4e,2);
            nrBasisFunc = size(obj.fe{1}.dof4e,2);
            for curFE = 2 : numel(obj.fe)
                obj.nrBasis4e = obj.nrBasis4e + obj.fe{curFE}.nrBasis4e;
                obj.idof4e = [obj.idof4e,nrBasisFunc+(size(obj.fe{curFE}.dof4e,2)-obj.fe{curFE}.nrBasis4e+1:size(obj.fe{1}.dof4e,2))];
                nrBasisFunc = nrBasisFunc + size(obj.fe{curFE}.dof4e,2);
            end
        end

        function C = curvedSymmetricBoundaryConstraints(obj,geom)
            %curvedSymmetricBoundaryConstraints - Build boundary constraints
            %   C = curvedSymmetricBoundaryConstraints(OBJ) returns
            %   the sparse constraint matrix of the symmetry boundary
            %   edges with a general normal (axis 0 in the last column
            %   of n4Sb) of a 2d mesh. C has nrDoF rows and one column
            %   per first-component dof on these edges; the column of
            %   the dof i holds the edge normal normal4ed in the rows
            %   of the dof i of both components (at a vertex shared by
            %   two such edges, the sum of their normals). Hence
            %   C'*U = 0 makes the normal component U*n of the
            %   coefficient vector U vanish at these dofs, e.g. imposed
            %   with a Lagrange multiplier. Without such edges C is the
            %   empty 0-by-0 sparse matrix.
            %
            %   C = curvedSymmetricBoundaryConstraints(OBJ,GEOM) uses
            %   at the two vertices of each edge the boundary normals
            %   GEOM.boundaryNormals of the geometry GEOM instead of
            %   the mesh normal.
            %
            %   In 3d, C is empty and faces with a general normal (axis
            %   0 in n4Sb) raise the error
            %   fem4m:PkVector:generalSymmetry3d.
            C = sparse(0,0);
            dim = size(obj.mesh.c4n,2);
            if ~isempty(obj.mesh.n4Sb) && dim > 2
                if any(obj.mesh.n4Sb(:,dim+1) == 0)
                    error('fem4m:PkVector:generalSymmetry3d', ...
                        'Symmetry boundary faces with a general normal (axis 0) are not supported in 3d.');
                end
            elseif ~isempty(obj.mesh.n4Sb)
                ed4Sb0 = obj.mesh.ed4Sb(obj.mesh.n4Sb(:,dim+1)==0);
                n4Sb0 = obj.mesh.n4Sb(obj.mesh.n4Sb(:,dim+1)==0,[1 2]);
                if ~isempty(ed4Sb0)
                    % all normals of circle
                    %X = obj.fe{1}.c4dof(:,1); Y = obj.fe{1}.c4dof(:,2);
                    %normals4ed = [X(obj.fe{1}.dof4ed(ed4Sb0,:)); ...
                    %    Y(obj.fe{1}.dof4ed(ed4Sb0,:))];
                    normals = obj.mesh.normal4ed;
                    normals = normals(ed4Sb0,:);
                    normals = normals(:).*ones(1,size(obj.fe{1}.dof4ed,2));
                    if nargin > 1 && ~isempty(geom)
                        n1 = geom.boundaryNormals(obj.mesh.c4n(n4Sb0(:,1),:));
                        n2 = geom.boundaryNormals(obj.mesh.c4n(n4Sb0(:,2),:));
                        normals(:,[1,end]) = [n1(:),n2(:)];
                    end
                    nodes4ed = [obj.fe{1}.dof4ed(ed4Sb0,:);obj.fe{1}.dof4ed(ed4Sb0,:)];
                    nodes4ed2 = [obj.fe{1}.dof4ed(ed4Sb0,:);obj.fe{1}.nrDoF+obj.fe{1}.dof4ed(ed4Sb0,:)];
                    C = sparse(nodes4ed2(:),nodes4ed(:),normals(:),obj.nrDoF,obj.nrDoF);
                    C = C(:,unique(nodes4ed(:)));
                end
            end
        end

    end

end