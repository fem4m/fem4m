classdef MeshPolyQ < MeshPoly
    %MeshPolyQ - Quadrilateral mesh stored as polygons with a quad mesh
    %   OBJ = MeshPolyQ(GEOM) creates the quadrilateral mesh of the
    %   Geometry object or options struct GEOM as MeshPoly does, with the
    %   Geometry option type set to 'Q' (a Geometry object of another
    %   type is rebuilt).
    %
    %   OBJ = MeshPolyQ(NAME) reads the files NAMEQ_*.dat if
    %   NAMEQ_c4n.dat exists and otherwise NAME_*.dat, as MeshPoly does,
    %   e.g. written by save, whose NAME_n4quad.dat and NAME_hnodes.dat
    %   restore meshes with hanging nodes. Without these files every
    %   element must be a quadrilateral.
    %
    %   OBJ = MeshPolyQ(c4n,n4e) creates the mesh with the node
    %   coordinates c4n and the node quadruples n4e of the elements; all
    %   boundary edges are Dirichlet edges.
    %
    %   OBJ = MeshPolyQ(...,n4Db) also specifies the Dirichlet boundary.
    %
    %   OBJ = MeshPolyQ(...,n4Nb) also specifies the Neumann boundary.
    %
    %   OBJ = MeshPolyQ(...,n4Sb) also specifies the symmetry boundary.
    %
    %   OBJ = MeshPolyQ(c4n,n4e,n4Db,n4Nb,n4Sb,n4quad,HNODES) also
    %   specifies the node quadruples n4quad of the underlying
    %   quadrilateral mesh and the sparse matrix HNODES of the midpoint
    %   nodes, as refine passes them; n4quad and HNODES are given only
    %   together.
    %
    %   Every element is the quadrilateral of the same row of n4quad,
    %   stored as the polygon of its vertices and of the hanging nodes
    %   on its edges. HNODES(i,j) = HNODES(j,i) is the id of the node at
    %   the midpoint of the quadrilateral edge between the nodes i and j,
    %   for every edge split by refine.
    %
    %   MeshPolyQ functions:
    %       mid4e  - Compute the vertex averages of all elements
    %       refine - Refine marked quads of the polygon mesh
    %       save   - Write the mesh data to .dat files
    %
    %   MeshPolyQ properties:
    %       n4quad - Node quadruples of the underlying quad mesh
    %       hnodes - Sparse matrix of midpoint node ids for node pairs
    %
    %   See also MeshPoly, MeshPolyT, MeshPolyH

    % Copyright (c) 2026 by the fem4m authors

    properties
        n4quad % node quadruples of the underlying quad mesh
        hnodes % sparse matrix of midpoint node ids for node pairs
    end

    methods
        function obj = MeshPolyQ(c4n,n4e,n4Db,n4Nb,n4Sb,n4quad,hnodes)
            %MeshPolyQ - Create a quadrilateral polygon mesh object
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
            if ischar(c4n) && exist([c4n,'Q_c4n.dat'],'file')
                % files of the Q type, NAMEQ_c4n.dat
                c4n = [c4n,'Q'];
            elseif isstruct(c4n)
                c4n.type = 'Q';
            elseif strcmpi(class(c4n),'Geometry')
                if c4n.opts.type ~= 'Q'
                    opts = c4n.opts;
                    opts.type = 'Q';
                    c4n = Geometry(opts);
                end
            end

            args = {c4n,n4e,n4Db,n4Nb,n4Sb};
            obj = obj@MeshPoly(args{1:given});

            if nargin < 6
                [obj.n4quad,obj.hnodes] = subMesh(obj,c4n,4,'n4quad');
            else
                obj.hnodes = hnodes;
                obj.n4quad = n4quad;
            end
        end

        function mid4e = mid4e(obj)
            %mid4e - Compute the vertex averages of all elements
            %   MID4E = mid4e(OBJ) computes the average of the four
            %   vertices of every quad (its centroid for parallelograms)
            %   and caches the result.
            mid4e = obj.mid4eSave;
            if isempty(mid4e)
                P1 = obj.c4n(obj.n4quad(:,1),:);
                P2 = obj.c4n(obj.n4quad(:,2),:);
                P3 = obj.c4n(obj.n4quad(:,3),:);
                P4 = obj.c4n(obj.n4quad(:,4),:);
                mid4e = (P1+P2+P3+P4)/4;
                obj.mid4eSave = mid4e;
            end
        end

        function obj = save(obj,filename)
            %SAVE - Write the mesh data to .dat files
            %   SAVE(OBJ) writes the files of MeshPoly.save and
            %   mesh_n4quad.dat and mesh_hnodes.dat with the underlying
            %   quadrilateral mesh and the midpoint nodes (none: no file,
            %   and the file of an earlier save is deleted), which
            %   MeshPolyQ('mesh') reads again.
            %
            %   SAVE(OBJ,FILENAME) also uses FILENAME as the prefix of the
            %   file names instead of 'mesh'.
            if nargin < 2 || isempty(filename)
                filename = 'mesh';
            end
            save@MeshPoly(obj,filename);
            writematrix(obj.n4quad,[filename,'_n4quad.dat'],'Delimiter',' ');
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
            %REFINE - Refine marked quads of the polygon mesh
            %   OBJ = REFINE(OBJ,R4E) inserts nodes at the midpoints of
            %   the edges of the quads marked by R4E (nrElems-by-1, true
            %   or positive for the marked elements), unless an earlier
            %   refinement inserted them, and at the vertex average mid4e
            %   of each marked quad, splits these quads into four quads,
            %   adds the new midpoints as hanging nodes to the polygons of
            %   the neighboring elements, and returns the refined mesh
            %   object. The refined mesh lists the unmarked elements first
            %   (in their old order), followed by four blocks with the
            %   first to fourth children of the marked elements.
            r4e = r4e>0;
            % insert new nodes
            quadEdges = [obj.n4quad(r4e,[1 2]);obj.n4quad(r4e,[2 3]);obj.n4quad(r4e,[3 4]);obj.n4quad(r4e,[4 1])];
            m4n = sparse(quadEdges(:,1),quadEdges(:,2),ones(4*nnz(r4e),1),obj.nrNodes,obj.nrNodes);
            m4n = triu(m4n + m4n' - obj.hnodes);
            [I,J] = find(m4n>0);
            m4n = sparse(I,J,obj.nrNodes + (1:numel(I)),obj.nrNodes,obj.nrNodes);
            m4n = m4n + m4n';
            mid4e = obj.mid4e;
            c4n = [obj.c4n; 1/2*(obj.c4n(I,:)+obj.c4n(J,:));mid4e(r4e,:)];
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
            n4quadRef = obj.n4quad;
            m1 = full(m4n(uint64(n4quadRef(r4e,1)) + size(m4n,1)*uint64(n4quadRef(r4e,2)-1)));
            m2 = full(m4n(uint64(n4quadRef(r4e,2)) + size(m4n,1)*uint64(n4quadRef(r4e,3)-1)));
            m3 = full(m4n(uint64(n4quadRef(r4e,3)) + size(m4n,1)*uint64(n4quadRef(r4e,4)-1)));
            m4 = full(m4n(uint64(n4quadRef(r4e,4)) + size(m4n,1)*uint64(n4quadRef(r4e,1)-1)));
            m5 = (size(c4n,1)-nnz(r4e)+1:size(c4n,1))';
            n4quadRef = [n4quadRef(~r4e,:);
                n4quadRef(r4e,1),m1,m5,m4;
                m1,n4quadRef(r4e,2),m2,m5;
                m2,n4quadRef(r4e,3),m3,m5;
                m3,n4quadRef(r4e,4),m4,m5];
            n4e = cell(size(n4quadRef,1),1);
            nrOld = nnz(r4e==0);
            nrNew = nnz(r4e);
            n4e(1:nrOld) = obj.n4e(~r4e);
            r4e = find(r4e);
            for curElem = 1:nnz(r4e)
                elem = obj.n4e{r4e(curElem)};
                indm1 = find(elem==m1(curElem));
                indm2 = find(elem==m2(curElem));
                indm3 = find(elem==m3(curElem));
                indm4 = find(elem==m4(curElem));
                n4e{nrOld+curElem} = [elem(1:indm1),m5(curElem),elem(indm4:end)];
                n4e{nrOld+nrNew+curElem} = [elem(indm1:indm2),m5(curElem)];
                n4e{nrOld+2*nrNew+curElem} = [elem(indm2:indm3),m5(curElem)];
                n4e{nrOld+3*nrNew+curElem} = [elem(indm3:indm4),m5(curElem)];
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
                sparse(size(c4n,1)-size(m4n,1),size(m4n,1)),...
                sparse(size(c4n,1)-size(m4n,1),size(c4n,1)-size(m4n,1))];
            % refined mesh
            obj = MeshPolyQ(c4n,n4e,n4Db,n4Nb,n4Sb,n4quadRef,m4n);
        end

    end
end