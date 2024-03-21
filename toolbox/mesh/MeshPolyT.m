classdef MeshPolyT < MeshPoly
    %MeshPolyT - Triangular mesh stored as polygons with a triangulation
    %   OBJ = MeshPolyT(GEOM) creates the triangular mesh of the Geometry
    %   object or options struct GEOM as MeshPoly does; every element of
    %   the geometry must be a triangle (Geometry option type 'T', the
    %   default).
    %
    %   OBJ = MeshPolyT(NAME) reads the files NAME_*.dat as MeshPoly
    %   does, e.g. written by save, whose NAME_n4tri.dat and
    %   NAME_hnodes.dat restore meshes with hanging nodes. Without these
    %   files every element must be a triangle.
    %
    %   OBJ = MeshPolyT(c4n,n4e) creates the mesh with the node
    %   coordinates c4n and the node triples n4e of the elements; all
    %   boundary edges are Dirichlet edges.
    %
    %   OBJ = MeshPolyT(...,n4Db) also specifies the Dirichlet boundary.
    %
    %   OBJ = MeshPolyT(...,n4Nb) also specifies the Neumann boundary.
    %
    %   OBJ = MeshPolyT(...,n4Sb) also specifies the symmetry boundary.
    %
    %   OBJ = MeshPolyT(c4n,n4e,n4Db,n4Nb,n4Sb,n4tri,HNODES) also
    %   specifies the node triples n4tri of the underlying triangulation
    %   and the sparse matrix HNODES of the midpoint nodes, as refine
    %   passes them; n4tri and HNODES are given only together.
    %
    %   Every element is the triangle of the same row of n4tri, stored as
    %   the polygon of its vertices and of the hanging nodes on its
    %   edges. HNODES(i,j) = HNODES(j,i) is the id of the node at the
    %   midpoint of the triangle edge between the nodes i and j, for
    %   every edge split by refine.
    %
    %   MeshPolyT functions:
    %       mid4e  - Compute the centroids of all elements
    %       refine - Refine marked triangles of the polygon mesh
    %       save   - Write the mesh data to .dat files
    %
    %   MeshPolyT properties:
    %       n4tri  - Node triples of the underlying triangulation
    %       hnodes - Sparse matrix of midpoint node ids for node pairs
    %
    %   See also MeshPoly, MeshPolyQ, MeshPolyH

    % Copyright (c) 2026 by the fem4m authors

    properties
        n4tri  % node triples of the underlying triangulation
        hnodes % sparse matrix of midpoint node ids for node pairs
    end

    methods
        function obj = MeshPolyT(c4n,n4e,n4Db,n4Nb,n4Sb,n4tri,hnodes)
            %MeshPolyT - Create a triangular polygon mesh object
            % pass on only the given arguments (MeshPoly has defaults)
            given = max(1,min(nargin,5));
            if nargin < 2
                n4e = [];
            end
            if nargin < 3
                n4Db = [];
            end
            if nargin < 4
                n4Nb = [];
            end
            if nargin < 5
                n4Sb = [];
            end

            args = {c4n,n4e,n4Db,n4Nb,n4Sb};
            obj = obj@MeshPoly(args{1:given});

            if nargin < 6
                [obj.n4tri,obj.hnodes] = subMesh(obj,c4n,3,'n4tri');
            else
                obj.hnodes = hnodes;
                obj.n4tri = n4tri;
            end
        end

        function mid4e = mid4e(obj)
            %mid4e - Compute the centroids of all elements
            %   MID4E = mid4e(OBJ) computes the centroid of every triangle
            %   as the average of its three vertices and caches the result.
            mid4e = obj.mid4eSave;
            if isempty(mid4e)
                P1 = obj.c4n(obj.n4tri(:,1),:);
                P2 = obj.c4n(obj.n4tri(:,2),:);
                P3 = obj.c4n(obj.n4tri(:,3),:);
                mid4e = (P1+P2+P3)/3;
                obj.mid4eSave = mid4e;
            end
        end

        function obj = save(obj,filename)
            %SAVE - Write the mesh data to .dat files
            %   SAVE(OBJ) writes the files of MeshPoly.save and
            %   mesh_n4tri.dat and mesh_hnodes.dat with the underlying
            %   triangulation and the midpoint nodes (none: no file, and
            %   the file of an earlier save is deleted), which
            %   MeshPolyT('mesh') reads again.
            %
            %   SAVE(OBJ,FILENAME) also uses FILENAME as the prefix of the
            %   file names instead of 'mesh'.
            if nargin < 2 || isempty(filename)
                filename = 'mesh';
            end
            save@MeshPoly(obj,filename);
            writematrix(obj.n4tri,[filename,'_n4tri.dat'],'Delimiter',' ');
            [i,j,v] = find(obj.hnodes);
            file = [filename,'_hnodes.dat'];
            if ~isempty(v)
                writematrix([i,j,v],file,'Delimiter',' ');
            elseif isfile(file)
                % no older hanging nodes are read back
                delete(file);
            end
        end

        function obj = refine(obj,r4e)
            %REFINE - Refine marked triangles of the polygon mesh
            %   OBJ = REFINE(OBJ,R4E) inserts nodes at the midpoints of
            %   the edges of the triangles marked by R4E (nrElems-by-1,
            %   true or positive for the marked elements), unless an
            %   earlier refinement inserted them, splits these triangles
            %   into four (red refinement), adds the new midpoints as
            %   hanging nodes to the polygons of the neighboring elements,
            %   and returns the refined mesh object. The refined mesh lists
            %   the unmarked elements first (in their old order), followed
            %   by four blocks with the first to fourth children of the
            %   marked elements.
            r4e = r4e>0;
            % insert new nodes
            triEdges = [obj.n4tri(r4e,[1 2]);obj.n4tri(r4e,[2 3]);obj.n4tri(r4e,[3 1])];
            m4n = sparse(triEdges(:,1),triEdges(:,2),ones(3*nnz(r4e),1),obj.nrNodes,obj.nrNodes);
            m4n = triu(m4n + m4n' - obj.hnodes);
            [I,J] = find(m4n>0);
            m4n = sparse(I,J,obj.nrNodes + (1:numel(I)),obj.nrNodes,obj.nrNodes);
            m4n = m4n + m4n';
            c4n = [obj.c4n; 1/2*(obj.c4n(I,:)+obj.c4n(J,:))];
            % include new nodes into polygons
            for curElem = 1 : size(obj.n4e,1)
                elem = obj.n4e{curElem};
                edges = obj.ed4e{curElem};
                shift = 0;
                for curEd = 1 : numel(edges)
                    if full(m4n(elem(curEd+shift),elem(mod(curEd+shift,numel(elem))+1)))>0
                        elem = [elem(1:curEd+shift),...
                            full(m4n(elem(curEd+shift),elem(mod(curEd+shift,numel(elem))+1))),...
                            elem(curEd+shift+1:end)];
                        shift = shift + 1;
                    end
                end
                obj.n4e{curElem} = elem;
            end
            m4n = m4n + obj.hnodes;
            % refine elements
            n4triRef = obj.n4tri;
            m1 = full(m4n(uint64(n4triRef(r4e,1)) + size(m4n,1)*uint64(n4triRef(r4e,2)-1)));
            m2 = full(m4n(uint64(n4triRef(r4e,2)) + size(m4n,1)*uint64(n4triRef(r4e,3)-1)));
            m3 = full(m4n(uint64(n4triRef(r4e,3)) + size(m4n,1)*uint64(n4triRef(r4e,1)-1)));
            n4triRef = [n4triRef(~r4e,:);
                n4triRef(r4e,1),m1,m3;
                m1,n4triRef(r4e,2),m2;
                m2,m3,m1;
                m3,m2,n4triRef(r4e,3)];
            n4e = cell(size(n4triRef,1),1);
            nrOld = nnz(r4e==0);
            nrNew = nnz(r4e);
            n4e(1:nrOld) = obj.n4e(~r4e);
            r4e = find(r4e);
            for curElem = 1:nnz(r4e)
                elem = obj.n4e{r4e(curElem)};
                indm1 = find(elem==m1(curElem));
                indm2 = find(elem==m2(curElem));
                indm3 = find(elem==m3(curElem));
                n4e{nrOld+curElem} = [elem(1:indm1),elem(indm3:end)];
                n4e{nrOld+nrNew+curElem} = elem(indm1:indm2);
                n4e{nrOld+2*nrNew+curElem} = [elem(indm2),elem(indm3),elem(indm1)];
                n4e{nrOld+3*nrNew+curElem} = [elem(indm3),elem(indm2:indm3-1)];
            end
            % refine boundary
            n4Db = obj.n4Db;
            if ~isempty(n4Db)
                r4Db = full(m4n(uint64(n4Db(:,1)) + size(m4n,1)*uint64(n4Db(:,2)-1))) > 0;
                m = full(m4n(uint64(n4Db(r4Db,1)) + size(m4n,1)*uint64(n4Db(r4Db,2)-1)));
                n4Db = [n4Db(~r4Db,:);n4Db(r4Db,1),m; m,n4Db(r4Db,2)];
            end
            n4Nb = obj.n4Nb;
            if ~isempty(n4Nb)
                r4Nb = full(m4n(uint64(n4Nb(:,1)) + size(m4n,1)*uint64(n4Nb(:,2)-1))) > 0;
                m = full(m4n(uint64(n4Nb(r4Nb,1)) + size(m4n,1)*uint64(n4Nb(r4Nb,2)-1)));
                n4Nb = [n4Nb(~r4Nb,:);n4Nb(r4Nb,1),m; m,n4Nb(r4Nb,2)];
            end
            n4Sb = obj.n4Sb;
            if ~isempty(n4Sb)
                r4Sb = full(m4n(uint64(n4Sb(:,1)) + size(m4n,1)*uint64(n4Sb(:,2)-1))) > 0;
                m = full(m4n(uint64(n4Sb(r4Sb,1)) + size(m4n,1)*uint64(n4Sb(r4Sb,2)-1)));
                n4Sb = [n4Sb(~r4Sb,:);n4Sb(r4Sb,1),m,n4Sb(r4Sb,3); m,n4Sb(r4Sb,2),n4Sb(r4Sb,3)];
            end
            % expand hnodes
            m4n = [m4n,sparse(size(m4n,1),size(c4n,1)-size(m4n,1));
                sparse(size(c4n,1)-size(m4n,1),size(m4n,1)),sparse(size(c4n,1)-size(m4n,1),size(c4n,1)-size(m4n,1))];
            % refined mesh
            obj = MeshPolyT(c4n,n4e,n4Db,n4Nb,n4Sb,n4triRef,m4n);
        end

    end
end