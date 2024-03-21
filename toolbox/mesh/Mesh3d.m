classdef Mesh3d
    %Mesh3d - Tetrahedral volume mesh with recursive newest vertex bisection
    %   OBJ = Mesh3d(GEOM) creates a tetrahedral mesh from the Geometry
    %   object or options struct GEOM.
    %
    %   OBJ = Mesh3d(NAME) loads the files NAME_c4n.dat, NAME_n4e.dat and
    %   the existing files NAME_n4Db.dat, NAME_n4Nb.dat, NAME_n4Sb.dat and
    %   NAME_t4e.dat.
    %
    %   OBJ = Mesh3d(C4N,N4E,N4DB) creates the mesh from the node
    %   coordinates C4N (nrNodes-by-3), the node ids N4E (nrElems-by-4)
    %   of the tetrahedra and the node triples N4DB of the Dirichlet
    %   boundary faces. Every boundary face must be in exactly one
    %   boundary list; otherwise an error is raised.
    %
    %   OBJ = Mesh3d(C4N,N4E,N4DB,N4NB) also specifies the node triples
    %   N4NB of the Neumann boundary faces.
    %
    %   OBJ = Mesh3d(C4N,N4E,N4DB,N4NB,N4SB) also specifies the symmetry
    %   (gliding) boundary faces as rows [n1 n2 n3 axis]: axis 1, 2 or 3
    %   is the fixed displacement component, i.e. the direction of the
    %   face normal, and 0 a general normal.
    %
    %   OBJ = Mesh3d(...,'t4e',T4E) also specifies the bisection labels:
    %   the vertex order of N4E and the labels T4E (see below). Without
    %   labels (also for NAME and GEOM without them) every element is
    %   bisected first at its longest edge and every face marks its
    %   longest edge, which gives the vertex order and the types (Arnold,
    %   Mukherjee, Pouly); edges of equal computed length (equal in
    %   floating point) are ordered by their node ids.
    %
    %   OBJ = Mesh3d(...,'init','coloring') takes the labels of a mesh
    %   without labels from a greedy coloring of the vertices instead:
    %   every element gets its vertices sorted by color, the largest
    %   color first, and type 0 (Diening, Gehring, Storn). The default
    %   is 'init','longestedge'.
    %
    %   The tetrahedra are refined by recursive bisection (Maubach,
    %   Traxler): an element (x1,x2,x3,x4) of type t is bisected at its
    %   refinement edge x1-x4 into (x1,m,x2,x3) and (x4,m,x2,x3), resp.
    %   (x4,m,x3,x2) for t = 0, both of type mod(t+1,3). The longest edge
    %   labels can also give the initial types M and O, elements
    %   (x1,x2,x3,x4) whose faces opposite x4 and x1 mark the edges x1-x3
    %   and x2-x3 (M) resp. x2-x3 and x2-x3 (O); their first bisection
    %   gives (x1,m,x2,x3) for M resp. (x2,m,x1,x3) for O and
    %   (x2,m,x4,x3), both of type 1. The stored n4e is this vertex
    %   order, with columns 2 and 3 swapped where the element would be
    %   negatively oriented, and t4e records the swap, t4e = type + 3*flip
    %   resp. 6 + 2*flip for M and 7 + 2*flip for O; so every element is
    %   positively oriented and the refinement edge is n4e(:,[1 4]).
    %   Every boundary face belongs to exactly one of n4Db, n4Nb and
    %   n4Sb.
    %
    %   Row y of every x4y lists the ids of the x of y. The incidences
    %   e4n and f4n (nrNodes x nrElems resp. nrFaces) and e4ed and f4ed
    %   (nrEdges x nrElems resp. nrFaces) hold the ids as values and are
    %   computed on access; ed4n is the symmetric node-pair lookup
    %   ed4n(i,j) of the edge between the nodes i and j.
    %
    %   Mesh3d functions:
    %       getBoundary     - Extract boundary faces of tetrahedra
    %       deleteNodes     - Delete nodes and their elements
    %       refine          - Refine by recursive bisection
    %       length4ed       - Compute the lengths of edges
    %       point4e         - Compute points in elements
    %       point4f         - Compute points on faces
    %       point4ed        - Compute points on edges
    %       normal4e        - Compute the outer face normals of elements
    %       normal4f        - Compute the outer normals of faces
    %       tangent4ed      - Compute the unit edge tangents
    %       area4e          - Compute the volumes of tetrahedra
    %       area4f          - Compute the areas of faces
    %       diam4e          - Compute the diameters of tetrahedra
    %       mid4e           - Compute the centroids of tetrahedra
    %       transformation  - Compute inverse mapping matrix
    %       jacobian        - Compute the mapping Jacobian
    %       color4n         - Color nodes joined by edges
    %       plot            - Plot the tetrahedral mesh
    %       plotBd          - Plot the boundary faces by type
    %       plotEnumeration - Label nodes, elements and faces with ids
    %       save            - Write the mesh data to .dat files
    %       disp            - Display a summary of the mesh
    %       bisection       - Removed, use refine
    %       redRefine       - Removed, use refine(obj,'uniform')
    %
    %   Mesh3d properties:
    %       c4n          - Coordinates of nodes
    %       n4e          - Node ids of the tetrahedra (bisection order)
    %       n4Db         - Node triples of Dirichlet boundary
    %       n4Nb         - Node triples of Neumann boundary
    %       n4Sb         - Node triples and axis of symmetry boundary
    %       t4e          - Bisection labels (type + 3*flip, M/O 6..9)
    %       nrNodes      - Number of nodes
    %       nrEdges      - Number of edges
    %       nrFaces      - Number of faces
    %       nrElems      - Number of elements
    %       f4e          - Face ids of the four faces of each element
    %       n4f          - Sorted node ids of the faces
    %       ed4f         - Edge ids of the three edges of a face
    %       e4f          - Adjacent element ids for each face (0 if none)
    %       e4ed         - Sparse matrix of element ids for edges
    %       ed4e         - Edge ids of the six edges of an element
    %       e4n          - Sparse matrix of element ids for nodes
    %       ed4n         - Sparse matrix of edge ids for node pairs
    %       f4n          - Sparse matrix of face ids for nodes
    %       f4ed         - Sparse matrix of face ids for edges
    %       n4ed         - Sorted node pairs of unique edges
    %       ed4Db        - Edge ids of Dirichlet boundary faces
    %       ed4Nb        - Edge ids of Neumann boundary faces
    %       ed4Sb        - Edge ids of symmetry boundary faces
    %       f4Db         - Face ids of Dirichlet boundary
    %       f4Nb         - Face ids of Neumann boundary
    %       f4Sb         - Face ids of symmetry boundary
    %       fNr4f        - Local face numbers of faces in elements
    %       parent4e     - Ancestor element before the last refine call
    %       prolongation - P1 prolongation from the mesh before refine
    %
    %   See also Mesh, Geometry, FE

    % Copyright (c) 2026 by the fem4m authors

    properties
        c4n                % coordinates of nodes
        n4e                % node ids of the tetrahedra (bisection order)
        n4Db               % node triples of Dirichlet boundary
        n4Nb               % node triples of Neumann boundary
        n4Sb               % node triples and axis of symmetry boundary
        t4e                % bisection labels (type + 3*flip, M/O 6..9)
        %
        nrNodes            % number of nodes
        nrEdges            % number of edges
        nrFaces            % number of faces
        nrElems            % number of elements
        %
        f4e                % face ids of the four faces of each element
        n4f                % sorted node ids of the faces
        ed4f               % edge ids of the three edges of a face
        e4f                % adjacent element ids for each face (0 if none)
        ed4e               % edge ids of the six edges of an element
        n4ed               % sorted node pairs of unique edges
        ed4Db              % edge ids of Dirichlet boundary faces
        ed4Nb              % edge ids of Neumann boundary faces
        ed4Sb              % edge ids of symmetry boundary faces
        f4Db               % face ids of Dirichlet boundary
        f4Nb               % face ids of Neumann boundary
        f4Sb               % face ids of symmetry boundary
        fNr4f              % local face numbers of faces in elements
        %
        parent4e           % ancestor element before the last refine call
        prolongation       % P1 prolongation from the mesh before refine
    end

    properties (Dependent)
        e4ed               % sparse matrix of element ids for edges
        e4n                % sparse matrix of element ids for nodes
        ed4n               % sparse matrix of edge ids for node pairs
        f4n                % sparse matrix of face ids for nodes
        f4ed               % sparse matrix of face ids for edges
    end

    properties (Dependent, Hidden)
        pkInterpolation    % removed: use FE.prolongation
    end

    methods(Static)

        function [n4f,e4f] = getBoundary(n4e)
            %getBoundary - Extract boundary faces of tetrahedra
            %   N4F = getBoundary(N4E) returns the node triples of all
            %   faces contained in exactly one tetrahedron of N4E, each
            %   sorted ascending (not oriented), as uint32.
            %
            %   [N4F,E4F] = getBoundary(N4E) also returns for each
            %   boundary face the id of its adjacent element.
            n4e = uint32(n4e);
            faces = Mesh3d.localSides(n4e,[1 2 3;1 2 4;4 2 3;1 3 4]);
            [p,st] = Mesh3d.sortedGroups(faces,max(n4e(:)));
            count = diff([find(st);numel(st)+1]);
            first = p(st);
            first = first(count == 1);
            n4f = faces(first,:);
            e4f = uint32(ceil(double(first)/4));
        end

        function [c4n,n4e,t4e,n4bd] = deleteNodes(c4n,n4e,t4e,nodes)
            %deleteNodes - Delete nodes and their elements
            %   C4N = deleteNodes(C4N,N4E,T4E,NODES) removes the nodes
            %   NODES and all tetrahedra containing them and returns the
            %   coordinates C4N of the remaining nodes.
            %
            %   [C4N,N4E,T4E] = deleteNodes(...) also returns the node
            %   ids N4E of the remaining elements, renumbered to the
            %   remaining nodes, and their bisection labels T4E.
            %
            %   [C4N,N4E,T4E,N4BD] = deleteNodes(...) also returns the
            %   sorted node triples N4BD of all boundary faces of the
            %   remaining tetrahedra (getBoundary).
            c4n(nodes,:) = [];
            nodes = sort(nodes,'descend');
            for curNode = nodes(:)'
                for k=1:4
                    indx = n4e(:,k)==curNode;
                    n4e(indx,:)=[];
                    t4e(indx) = [];
                end
                n4e(n4e(:)>curNode) = n4e(n4e(:)>curNode) - 1;
            end
            n4bd = Mesh3d.getBoundary(n4e);
        end

    end

    methods
        function obj = Mesh3d(c4n,varargin)
            %Mesh3d - Create a tetrahedral mesh object
            t4e = [];
            n4Sb = [];
            init = 'longestedge';
            rest = varargin;
            if ischar(c4n) || isstring(c4n) || isstruct(c4n) || isa(c4n,'Geometry')
                if ischar(c4n) || isstring(c4n)
                    name = char(c4n);
                    obj.c4n = load([name,'_c4n.dat']);
                    obj.n4e = uint32(load([name,'_n4e.dat']));
                    obj.n4Db = Mesh3d.loadList([name,'_n4Db.dat']);
                    obj.n4Nb = Mesh3d.loadList([name,'_n4Nb.dat']);
                    n4Sb = Mesh3d.loadList([name,'_n4Sb.dat']);
                    if exist([name,'_t4e.dat'],'file')
                        t4e = load([name,'_t4e.dat']);
                    end
                else
                    geom = c4n;
                    if isstruct(geom)
                        geom = Geometry(geom);
                    end
                    obj.c4n = geom.c4n;
                    obj.n4e = uint32(geom.n4e);
                    obj.n4Db = uint32(geom.n4Db);
                    obj.n4Nb = uint32(geom.n4Nb);
                    n4Sb = geom.n4Sb;
                    t4e = geom.t4e;
                end
            else
                % n4e, n4Db and n4Nb up to the first option name
                isName = cellfun(@(a)ischar(a) || isstring(a),[rest,{''}]);
                nLists = min(find(isName,1) - 1,3);
                if nLists == 0
                    error('fem4m:Mesh3d:missingN4e', ...
                        'Mesh3d(c4n,n4e,n4Db,...) needs the node ids n4e of the elements.');
                end
                lists = [rest(1:nLists),cell(1,3-nLists)];
                rest = rest(nLists+1:end);
                if ~isempty(rest) && ~(ischar(rest{1}) || isstring(rest{1}))
                    n4Sb = rest{1};
                    rest = rest(2:end);
                end
                obj.c4n = c4n;
                obj.n4e = uint32(lists{1});
                obj.n4Db = uint32(lists{2});
                obj.n4Nb = uint32(lists{3});
            end
            while ~isempty(rest)
                if ~(ischar(rest{1}) || isstring(rest{1}))
                    Mesh3d.pkInterpolationRemoved();
                elseif strcmpi(rest{1},'t4e') && numel(rest) > 1
                    t4e = rest{2};
                elseif strcmpi(rest{1},'init') && numel(rest) > 1
                    init = rest{2};
                    if ~((ischar(init) || isstring(init)) && any(strcmpi(init,{'longestedge','coloring'})))
                        error('fem4m:Mesh3d:unknownOption', ...
                            'The option ''init'' takes ''longestedge'' or ''coloring''.');
                    end
                    init = lower(char(init));
                else
                    error('fem4m:Mesh3d:unknownOption', ...
                        ['Unknown option %s; Mesh3d takes the name-value pairs ''t4e'' ', ...
                        'and ''init'' (''longestedge'' or ''coloring'').'],char(rest{1}));
                end
                rest = rest(3:end);
            end
            obj.n4Sb = Mesh3d.checkN4Sb(n4Sb);
            if isempty(t4e) && strcmp(init,'coloring')
                [obj.n4e,obj.t4e] = Mesh3d.colorLabels(obj.n4e,size(obj.c4n,1));
            elseif isempty(t4e)
                [obj.n4e,obj.t4e] = Mesh3d.longestEdgeLabels(obj.c4n,obj.n4e);
            else
                if numel(t4e) ~= size(obj.n4e,1) || any(t4e(:) ~= round(t4e(:))) ...
                        || any(t4e(:) < 0 | t4e(:) > 9)
                    error('fem4m:Mesh3d:invalidT4e', ...
                        't4e must hold one bisection label 0..9 per element.');
                end
                obj.t4e = int8(t4e(:));
            end
            obj = normalizeOrientation(obj);
            obj = enumerate(obj);
            obj = fixBoundaryOrientation(obj);
        end

        function val = get.e4n(obj)
            %e4n - Sparse matrix of element ids for nodes
            J = (1:size(obj.n4e,1))'*ones(1,size(obj.n4e,2));
            val = sparse(double(obj.n4e(:)),J(:),J(:),obj.nrNodes,obj.nrElems);
        end

        function val = get.f4n(obj)
            %f4n - Sparse matrix of face ids for nodes
            J = (1:size(obj.n4f,1))'*ones(1,size(obj.n4f,2));
            val = sparse(double(obj.n4f(:)),J(:),J(:),obj.nrNodes,obj.nrFaces);
        end

        function val = get.ed4n(obj)
            %ed4n - Sparse matrix of edge ids for node pairs
            val = sparse(double(obj.n4ed(:,[1 2])),double(obj.n4ed(:,[2 1])),[1:obj.nrEdges,1:obj.nrEdges]', ...
                obj.nrNodes,obj.nrNodes);
        end

        function val = get.e4ed(obj)
            %e4ed - Sparse matrix of element ids for edges
            J = (1:size(obj.ed4e,1))'*ones(1,size(obj.ed4e,2));
            val = sparse(double(obj.ed4e(:)),J(:),J(:),obj.nrEdges,obj.nrElems);
        end

        function val = get.f4ed(obj)
            %f4ed - Sparse matrix of face ids for edges
            J = (1:size(obj.ed4f,1))'*ones(1,size(obj.ed4f,2));
            val = sparse(double(obj.ed4f(:)),J(:),J(:),obj.nrEdges,obj.nrFaces);
        end

        function length4ed = length4ed(obj)
            %length4ed - Compute the lengths of edges
            %   LENGTH4ED = length4ed(OBJ) computes the Euclidean length
            %   of every edge listed in OBJ.n4ed (nrEdges-by-1).
            coord1 = obj.c4n(obj.n4ed(:,1),:);
            coord2 = obj.c4n(obj.n4ed(:,2),:);
            length4ed = sqrt(sum((coord2-coord1).^2,2));
        end

        function point4e = point4e(obj,refp)
            %point4e - Compute points in elements
            %   POINT4E = point4e(OBJ,REFP) computes for every
            %   tetrahedron the point with the reference coordinates REFP
            %   (1-by-3) with respect to the first node and the three edge
            %   directions emanating from it (nrElems-by-3).
            P1 = obj.c4n(obj.n4e(:,1),:);
            P2 = obj.c4n(obj.n4e(:,2),:);
            P3 = obj.c4n(obj.n4e(:,3),:);
            P4 = obj.c4n(obj.n4e(:,4),:);
            point4e = P1 + refp(1)*(P2-P1) + refp(2)*(P3-P1) + refp(3)*(P4-P1);
        end

        function point4f = point4f(obj,refp)
            %point4f - Compute points on faces
            %   POINT4F = point4f(OBJ,REFP) computes for every face the
            %   point with the reference coordinates REFP (1-by-2) with
            %   respect to its first node n4f(:,1) and the two edge
            %   directions to n4f(:,2) and n4f(:,3) (nrFaces-by-3).
            P1 = obj.c4n(obj.n4f(:,1),:);
            P2 = obj.c4n(obj.n4f(:,2),:);
            P3 = obj.c4n(obj.n4f(:,3),:);
            point4f = P1 + refp(1)*(P2-P1) + refp(2)*(P3-P1);
        end

        function point4ed = point4ed(obj,refp)
            %point4ed - Compute points on edges
            %   POINT4ED = point4ed(OBJ,REFP) computes for every edge
            %   the point with the scalar reference coordinate REFP,
            %   where 0 and 1 are its end nodes n4ed(:,1) and n4ed(:,2)
            %   (nrEdges-by-3).
            P1 = obj.c4n(obj.n4ed(:,1),:);
            P2 = obj.c4n(obj.n4ed(:,2),:);
            point4ed = P1 + refp(1)*(P2-P1);
        end

        function normal4e = normal4e(obj)
            %normal4e - Compute the outer face normals of elements
            %   NORMAL4E = normal4e(OBJ) computes for every tetrahedron
            %   the outward unit normal vectors of its four faces
            %   (nrElems-by-4-by-3); NORMAL4E(:,j,:) belongs to the face
            %   f4e(:,j).
            %outer normals
            perm = [1 2 3;1 4 2;2 4 3;1 3 4];
            normal4e = zeros(obj.nrElems,4,3);
            for j = 1 : 4
                curNormal4e = cross(obj.c4n(obj.n4e(:,perm(j,3)),:)-obj.c4n(obj.n4e(:,perm(j,1)),:),...
                    obj.c4n(obj.n4e(:,perm(j,2)),:)-obj.c4n(obj.n4e(:,perm(j,1)),:));
                normal4e(:,j,:) = curNormal4e./((sqrt(sum(curNormal4e.^2,2)))*ones(1,3));
            end
        end

        function normal4f = normal4f(obj)
            %normal4f - Compute the outer normals of faces
            %   NORMAL4F = normal4f(OBJ) computes the outward unit normal
            %   vector of every face with respect to its first adjacent
            %   element e4f(:,1) (nrFaces-by-3).
            %outer normals with respect to e4f(:,1)
            perm = [1 2 3;1 4 2;2 4 3;1 3 4];
            nodes4f = [rowaddr(obj.n4e,obj.e4f(:,1),perm(obj.fNr4f(:,1),1)), ...
                rowaddr(obj.n4e,obj.e4f(:,1),perm(obj.fNr4f(:,1),2)), ...
                rowaddr(obj.n4e,obj.e4f(:,1),perm(obj.fNr4f(:,1),3))];
            normal4f = cross(obj.c4n(nodes4f(:,3),:)-obj.c4n(nodes4f(:,1),:),...
                obj.c4n(nodes4f(:,2),:)-obj.c4n(nodes4f(:,1),:));
            normal4f = normal4f./((sqrt(sum(normal4f.^2,2)))*ones(1,3));
        end

        function tangent4ed = tangent4ed(obj)
            %tangent4ed - Compute the unit edge tangents
            %   TANGENT4ED = tangent4ed(OBJ) computes the unit tangent
            %   vector of every edge pointing from its first to its
            %   second node in OBJ.n4ed (nrEdges-by-3).
            c4start = obj.c4n(obj.n4ed(:,1),:);
            c4end = obj.c4n(obj.n4ed(:,2),:);
            l4ed = obj.length4ed;
            tangent4ed = (c4end-c4start)./[l4ed l4ed l4ed];
        end

        function area4e = area4e(obj)
            %area4e - Compute the volumes of tetrahedra
            %   AREA4E = area4e(OBJ) computes for every tetrahedron the
            %   signed volume from the coordinates of its four nodes; it
            %   is positive, as every element is positively oriented
            %   (nrElems-by-1).
            % area = zeros(size(n4e,1),1);
            %
            % for curElem = 1 : size(n4e,1)
            %     area(curElem) = 1/6*det([ 1 1 1 1; c4n(n4e(curElem,:),:)']);
            % end
            % syms x1 x2 x3 x4 y1 y2 y3 y4 z1 z2 z3 z4
            %
            % det([1 1 1 1; x1 x2 x3 x4; y1 y2 y3 y4; z1 z2 z3 z4])
            x = obj.c4n(:,1); x = x(obj.n4e);
            y = obj.c4n(:,2); y = y(obj.n4e);
            z = obj.c4n(:,3); z = z(obj.n4e);
            if size(obj.n4e,1) == 1
                x=x';y=y';z=z';
            end
            area4e = x(:,2).*y(:,3).*z(:,4)-x(:,2).*y(:,4).*z(:,3)-y(:,2).*x(:,3).*z(:,4)+y(:,2).*x(:,4).*z(:,3)+z(:,2).*x(:,3).*y(:,4)-z(:,2).*x(:,4).*y(:,3)-x(:,1).*y(:,3).*z(:,4)...
                +x(:,1).*y(:,4).*z(:,3)+x(:,1).*y(:,2).*z(:,4)-x(:,1).*y(:,2).*z(:,3)-x(:,1).*z(:,2).*y(:,4)+x(:,1).*z(:,2).*y(:,3)+y(:,1).*x(:,3).*z(:,4)...
                -y(:,1).*x(:,4).*z(:,3)-y(:,1).*x(:,2).*z(:,4)+y(:,1).*x(:,2).*z(:,3)+y(:,1).*z(:,2).*x(:,4)-y(:,1).*z(:,2).*x(:,3)-z(:,1).*x(:,3).*y(:,4)...
                +z(:,1).*x(:,4).*y(:,3)+z(:,1).*x(:,2).*y(:,4)-z(:,1).*x(:,2).*y(:,3)-z(:,1).*y(:,2).*x(:,4)+z(:,1).*y(:,2).*x(:,3);
            area4e = area4e/6;
        end

        function area4f = area4f(obj)
            %area4f - Compute the areas of faces
            %   AREA4F = area4f(OBJ) computes for every face the area of
            %   the triangle spanned by its three nodes (nrFaces-by-1).
            area4f = 1/2 * sqrt(sum((cross(obj.c4n(obj.n4f(:,3),:)-obj.c4n(obj.n4f(:,1),:),...
                obj.c4n(obj.n4f(:,2),:)-obj.c4n(obj.n4f(:,1),:))).^2,2));
        end

        function diam4e = diam4e(obj)
            %diam4e - Compute the diameters of tetrahedra
            %   DIAM4E = diam4e(OBJ) computes the diameter of every
            %   tetrahedron as the length of its longest edge
            %   (nrElems-by-1).
            length4ed = obj.length4ed;
            diam4e = max(length4ed(obj.ed4e),[],2);
        end

        function mid4e = mid4e(obj)
            %mid4e - Compute the centroids of tetrahedra
            %   MID4E = mid4e(OBJ) computes the centroid of every
            %   tetrahedron as the mean of its four nodes (nrElems-by-3).
            mid4e = obj.point4e([1/4,1/4,1/4]);
        end

        function val = transformation(obj)
            %TRANSFORMATION - Compute inverse mapping matrix
            %   VAL = TRANSFORMATION(OBJ) computes for every tetrahedron
            %   the inverse C = inv(B) of the Jacobian B = [x2 x3 x4;
            %   y2 y3 y4; z2 z3 z4] (see jacobian) of the affine map from
            %   the reference tetrahedron. Row e of VAL (nrElems-by-9)
            %   holds C of element e column by column,
            %   [c11,c21,c31,c12,c22,c32,c13,c23,c33], so
            %   reshape(VAL,[],3,3) gives C(i,j) at (:,i,j).
            P1 = obj.c4n(obj.n4e(:,1),:);
            P2 = obj.c4n(obj.n4e(:,2),:);
            P3 = obj.c4n(obj.n4e(:,3),:);
            P4 = obj.c4n(obj.n4e(:,4),:);
            x2 = P2(:,1)-P1(:,1); y2 = P2(:,2)-P1(:,2); z2 = P2(:,3)-P1(:,3);
            x3 = P3(:,1)-P1(:,1); y3 = P3(:,2)-P1(:,2); z3 = P3(:,3)-P1(:,3);
            x4 = P4(:,1)-P1(:,1); y4 = P4(:,2)-P1(:,2); z4 = P4(:,3)-P1(:,3);

            T = 6*obj.area4e;
            % [c11,c21,c31,c12,c22,c32,c13,c23,c33]
            val = [(y3.*z4-y4.*z3)./T,(-y2.*z4+y4.*z2)./T,(y2.*z3-y3.*z2)./T,...
                (-x3.*z4+x4.*z3)./T,(x2.*z4-x4.*z2)./T,(-x2.*z3+x3.*z2)./T,...
                (x3.*y4-x4.*y3)./T,(-x2.*y4+x4.*y2)./T,(x2.*y3-x3.*y2)./T];
        end

        function val = jacobian(obj)
            %JACOBIAN - Compute the mapping Jacobian
            %   VAL = JACOBIAN(OBJ) computes for every tetrahedron the
            %   entries of the Jacobian B = [P2-P1,P3-P1,P4-P1] of the
            %   affine map from the reference tetrahedron, P1 to P4 the
            %   coordinates (as columns) of the nodes n4e(e,:). Row e of
            %   VAL (nrElems-by-9) holds B row by row as
            %   [x2,x3,x4,y2,y3,y4,z2,z3,z4], xk, yk, zk the components
            %   of Pk-P1.
            P1 = obj.c4n(obj.n4e(:,1),:);
            P2 = obj.c4n(obj.n4e(:,2),:);
            P3 = obj.c4n(obj.n4e(:,3),:);
            P4 = obj.c4n(obj.n4e(:,4),:);
            x2 = P2(:,1)-P1(:,1); y2 = P2(:,2)-P1(:,2); z2 = P2(:,3)-P1(:,3);
            x3 = P3(:,1)-P1(:,1); y3 = P3(:,2)-P1(:,2); z3 = P3(:,3)-P1(:,3);
            x4 = P4(:,1)-P1(:,1); y4 = P4(:,2)-P1(:,2); z4 = P4(:,3)-P1(:,3);

            % B row by row, i.e. the entries of B' column by column
            val = [x2,x3,x4,y2,y3,y4,z2,z3,z4];
        end

        function color4n = color4n(obj,levels)
            %color4n - Color nodes joined by edges
            %   COLOR4N = color4n(OBJ) assigns colors 1, 2, ... to the
            %   nodes (nrNodes-by-1) such that nodes joined by an edge
            %   receive different colors. The nodes are colored greedily
            %   in the order of their ids with the smallest free color.
            %
            %   COLOR4N = color4n(OBJ,LEVELS) treats nodes connected by
            %   paths of up to LEVELS edges as neighbors (default 1, also
            %   for an empty LEVELS).
            if nargin < 2 || isempty(levels)
                levels = 1;
            end
            color4n = zeros(obj.nrNodes,1);
            ed4nodes = obj.ed4n;
            for curNode = 1 : obj.nrNodes
                neighbours = curNode;
                for cuLevel = 1 : levels
                    [row,~] = find(ed4nodes(:,neighbours));
                    neighbours = [neighbours;row]; %#ok<AGROW>
                end
                neighbour_colors = color4n(neighbours);
                color4n(curNode) = min(setdiff(1:max(neighbour_colors)+1,neighbour_colors));
            end
        end

        function val = get.pkInterpolation(~) %#ok<STOUT>
            %pkInterpolation - Removed, use FE.prolongation
            Mesh3d.pkInterpolationRemoved();
        end

        function obj = refine(obj,varargin)
            %REFINE - Refine by recursive bisection
            %   OBJ = REFINE(OBJ,MARKED4E) bisects every element marked by
            %   the logical or 0/1 vector MARKED4E of length nrElems once
            %   at its refinement edge n4e(:,[1 4]) and bisects further
            %   elements until the mesh is conforming. Unlike in 2d, a
            %   marked element is bisected only once.
            %
            %   OBJ = REFINE(OBJ,MARKED4ED) with a vector of length
            %   nrEdges halves the marked edges, as in 2d. Element marks
            %   win if nrElems equals nrEdges.
            %
            %   OBJ = REFINE(OBJ,'uniform') bisects every element three
            %   times and then completes the mesh. On the initial or a
            %   uniformly refined mesh every element gets exactly 8
            %   children; after adaptive steps the completion can bisect
            %   elements next to earlier refinement further.
            %
            %   OBJ = REFINE(OBJ,'uniform',N) repeats this N times.
            %
            %   parent4e and the P1 prolongation always refer to the mesh
            %   before the call. Marks of another length raise an error.
            if isempty(varargin)
                error('fem4m:Mesh3d:unknownRefineInput', ...
                    'refine needs element or edge marks or ''uniform''.');
            end
            marks = varargin{1};
            if ischar(marks) || isstring(marks)
                if ~strcmpi(marks,'uniform')
                    error('fem4m:Mesh3d:unknownRefineInput', ...
                        'Unknown refinement %s; use marks or ''uniform''.',char(marks));
                end
                rounds = 1;
                if numel(varargin) > 1
                    rounds = varargin{2};
                end
                parent4e = (1:obj.nrElems)';
                P = speye(obj.nrNodes);
                for k = 1 : rounds
                    obj = bisect(obj,[],[],3);
                    parent4e = parent4e(obj.parent4e);
                    P = obj.prolongation*P;
                end
                obj.parent4e = parent4e;
                obj.prolongation = P;
            elseif numel(marks) == obj.nrElems
                obj = bisect(obj,marks(:) ~= 0,[],0);
            elseif numel(marks) == obj.nrEdges
                obj = bisect(obj,[],marks(:) ~= 0,0);
            else
                error('fem4m:Mesh3d:unknownRefineInput', ...
                    'Marks must be a vector of length nrElems (%d) or nrEdges (%d).', ...
                    obj.nrElems,obj.nrEdges);
            end
        end

        function obj = bisection(obj,varargin) %#ok<INUSD>
            %BISECTION - Removed, use refine
            %   OBJ = BISECTION(OBJ,...) raises an error; use
            %   OBJ = REFINE(OBJ,MARKED4E) instead.
            error('fem4m:Mesh3d:bisectionRemoved', ...
                'bisection was replaced by refine: use mesh.refine(marked4e).');
        end

        function obj = redRefine(obj) %#ok<MANU>
            %redRefine - Removed, use refine(obj,'uniform')
            %   OBJ = redRefine(OBJ) raises an error; use
            %   OBJ = REFINE(OBJ,'uniform') instead.
            error('fem4m:Mesh3d:redRefineRemoved', ...
                'redRefine was replaced by uniform bisection: use mesh.refine(''uniform'').');
        end

        function plot(obj,opt)
            %PLOT - Plot the tetrahedral mesh
            %   PLOT(OBJ) clears the current figure and draws all
            %   tetrahedra with tetramesh if the mesh has at most 500
            %   elements, and the boundary surface (option "2d")
            %   otherwise.
            %
            %   PLOT(OBJ,OPT) selects the view. OPT must contain one of
            %   the following, checked in this order (the first match
            %   wins):
            %       "3d" - Plot all tetrahedra with tetramesh.
            %       "2d" - Plot the boundary faces as white patches.
            %       "bd" - Plot the boundary faces by type (plotBd).
            %   Other values give a warning.
            if nargin < 2
                opt = '2d';
            end
            if (nargin < 2 && obj.nrElems <= 500) || contains(opt,'3d')
                clf;
                tetramesh(double(obj.n4e),obj.c4n);
            elseif nargin < 2 || contains(opt,'2d')
                n4fBd = obj.n4f(obj.e4f(:,2) == 0,:);
                coordX = reshape(obj.c4n(n4fBd',1),3,[]);
                coordY = reshape(obj.c4n(n4fBd',2),3,[]);
                coordZ = reshape(obj.c4n(n4fBd',3),3,[]);
                clf;
                patch(coordX,coordY,coordZ,'w');
                axis equal;
                axis([min(coordX(:)),max(coordX(:)),min(coordY(:)),max(coordY(:)),min(coordZ(:)),max(coordZ(:))]);
            elseif contains(opt,'bd')
                clf;
                obj.plotBd();
            else
                warning('fem4m:Mesh3d:unknownPlotOption','option: %s not recognized',opt);
            end
        end

        function plotBd(obj)
            %plotBd - Plot the boundary faces by type
            %   plotBd(OBJ) draws the Dirichlet faces red, the Neumann
            %   faces blue and the symmetry faces green into the current
            %   axes.
            lists = {obj.n4Db,[1 0 0];obj.n4Nb,[0 0 1];obj.n4Sb,[0 0.6 0]};
            hold on
            for k = 1 : size(lists,1)
                n4b = lists{k,1};
                if isempty(n4b)
                    continue
                end
                n4b = double(n4b(:,1:3))';
                patch(reshape(obj.c4n(n4b,1),3,[]),reshape(obj.c4n(n4b,2),3,[]), ...
                    reshape(obj.c4n(n4b,3),3,[]),lists{k,2},'FaceAlpha',0.8);
            end
            hold off
            axis equal
            view(3)
        end

        function plotEnumeration(obj)
            %plotEnumeration - Label nodes, elements and faces with ids
            %   plotEnumeration(OBJ) writes the ids of all nodes,
            %   elements and faces as text at the nodes and at the
            %   centroids of the elements and faces into the current
            %   axes; it does not draw the mesh (see plot).
            mid4e = obj.mid4e;
            mid4f = obj.point4f([1/3,1/3]);
            text(obj.c4n(:,1),obj.c4n(:,2),obj.c4n(:,3),num2str((1:size(obj.c4n,1))'),...
                'color',[0.26 0.49 0.19],'FontWeight','bold','FontSize',24,'HorizontalAlignment','center');
            text(mid4e(:,1),mid4e(:,2),mid4e(:,3),num2str((1:size(obj.n4e,1))'),...
                'color',[0.95 0.91 0],'FontWeight','bold','FontSize',24,'HorizontalAlignment','center');
            text(mid4f(:,1),mid4f(:,2),mid4f(:,3),num2str((1:obj.nrFaces)'),...
                'color',[0.38 0.02 0.68],'FontWeight','bold','FontSize',24,'HorizontalAlignment','center');
        end

        function obj = save(obj,filename)
            %SAVE - Write the mesh data to .dat files
            %   SAVE(OBJ,FILENAME) writes for the char vector FILENAME
            %   the files FILENAME_c4n.dat (17 significant digits),
            %   FILENAME_n4e.dat, FILENAME_t4e.dat and the boundary lists
            %   FILENAME_n4Db.dat, FILENAME_n4Nb.dat and FILENAME_n4Sb.dat,
            %   which Mesh3d(FILENAME) reads again. An empty boundary list
            %   writes no file and deletes the file of an earlier save.
            % all 17 digits, writematrix alone writes 15
            writematrix(compose("%.17g",obj.c4n),[filename,'_c4n.dat'], ...
                'Delimiter',' ','QuoteStrings','none');
            writematrix(obj.n4e,[filename,'_n4e.dat'],'Delimiter',' ');
            writematrix(obj.t4e,[filename,'_t4e.dat'],'Delimiter',' ');
            names = {'n4Db','n4Nb','n4Sb'};
            for k = 1 : numel(names)
                file = [filename,'_',names{k},'.dat'];
                if ~isempty(obj.(names{k}))
                    writematrix(obj.(names{k}),file,'Delimiter',' ');
                elseif isfile(file)
                    % no older list is read back
                    delete(file);
                end
            end
        end

        function disp(obj)
            %DISP - Display a summary of the mesh
            %   DISP(OBJ) prints the numbers of nodes, edges, faces,
            %   elements and boundary faces.
            fprintf('  Mesh3d with %d nodes, %d edges, %d faces and %d elements\n', ...
                obj.nrNodes,obj.nrEdges,obj.nrFaces,obj.nrElems);
            fprintf('  boundary faces: %d Dirichlet, %d Neumann, %d symmetry\n', ...
                size(obj.n4Db,1),size(obj.n4Nb,1),size(obj.n4Sb,1));
        end

    end

    methods (Access=private)

        function obj = normalizeOrientation(obj)
            %normalizeOrientation - Positive orientation with the flip bit
            %   Where an element is negatively oriented, its columns 2 and
            %   3 are swapped and its flip bit is toggled.
            negative = Mesh3d.signedVolumes(obj.c4n,obj.n4e) < 0;
            obj.n4e(negative,[2 3]) = obj.n4e(negative,[3 2]);
            [kind,flip] = Mesh3d.decodeLabels(obj.t4e);
            obj.t4e = Mesh3d.encodeLabels(kind,xor(flip,negative));
        end

        function [nodes4e,type,label] = maubachForm(obj)
            %maubachForm - Bisection order, types and boundary labels
            %   NODES4E is the vertex order of the bisection (refinement
            %   edge nodes4e(:,[1 4])), TYPE the bisection type and
            %   LABEL(e,k) the boundary code of the face opposite vertex
            %   k: 0 interior, 1 Dirichlet, 2 Neumann, 3+axis symmetry.
            code4f = zeros(obj.nrFaces,1);
            code4f(obj.f4Db) = 1;
            code4f(obj.f4Nb) = 2;
            if ~isempty(obj.n4Sb)
                code4f(obj.f4Sb) = 3 + double(obj.n4Sb(:,4));
            end
            label = reshape(code4f(obj.f4e(:,[3 4 2 1])),[],4);
            nodes4e = double(obj.n4e);
            [type,flip] = Mesh3d.decodeLabels(obj.t4e);
            nodes4e(flip,[2 3]) = nodes4e(flip,[3 2]);
            label(flip,[2 3]) = label(flip,[3 2]);
        end

        function obj = bisect(obj,marked4e,req4ed,fullSweeps)
            %BISECT - Recursive bisection with closure
            %   OBJ = BISECT(OBJ,MARKED4E,REQ4ED,FULLSWEEPS) bisects the
            %   elements MARKED4E (logical over the elements) once, the
            %   elements of the edges REQ4ED (logical over the edges)
            %   until these edges are halved, and in the first FULLSWEEPS
            %   sweeps all elements; in every sweep the elements with an
            %   edge that is halved are bisected as well, until the mesh
            %   is conforming. Edges are identified by the keys
            %   min*2^32+max of their node ids.
            K = uint64(2^32);
            localEdges = [1 2;2 3;3 1;1 4;2 4;3 4];
            c4n = obj.c4n;
            nN0 = size(c4n,1);
            [nodes4e,type,label] = obj.maubachForm();
            parent = (1:obj.nrElems)';
            if isempty(marked4e)
                marked4e = false(obj.nrElems,1);
            end
            reqKeys = zeros(0,1,'uint64');
            if ~isempty(req4ed)
                reqKeys = Mesh3d.edgeKeys(double(obj.n4ed(req4ed,:)),K);
            end
            midKeys = zeros(0,1,'uint64');
            Pt = speye(nN0);
            for sweep = 1 : 50 + 3*fullSweeps
                keys = zeros(size(nodes4e,1),6,'uint64');
                for j = 1 : 6
                    keys(:,j) = Mesh3d.edgeKeys(nodes4e(:,localEdges(j,:)),K);
                end
                reqKeys = reqKeys(~ismember(reqKeys,midKeys));
                if sweep <= fullSweeps
                    mark = true(size(nodes4e,1),1);
                else
                    mark = any(ismember(keys,[midKeys;reqKeys]),2);
                    if sweep == 1
                        mark = mark | marked4e;
                    end
                end
                if ~any(mark)
                    break
                end
                % closure within the sweep: elements with a halved edge
                halved = [midKeys;keys(mark,4)];
                while true
                    more = ~mark & any(ismember(keys,halved),2);
                    if ~any(more)
                        break
                    end
                    mark = mark | more;
                    halved = [halved;keys(more,4)]; %#ok<AGROW>
                end
                % midpoints of the refinement edges, new or from an earlier sweep
                [refKeys,~,refId] = unique(keys(mark,4));
                newKeys = refKeys(~ismember(refKeys,midKeys));
                lo = double(idivide(newKeys,K,'floor'));
                hi = double(newKeys - uint64(lo)*K);
                c4n = [c4n;(c4n(lo,:) + c4n(hi,:))/2]; %#ok<AGROW>
                Pt = [Pt,(Pt(:,lo) + Pt(:,hi))/2]; %#ok<AGROW>
                midKeys = [midKeys;newKeys]; %#ok<AGROW>
                [~,pos] = ismember(refKeys,midKeys);
                mid = nN0 + pos;
                m = mid(refId);
                % children (Traxler): (x1,m,x2,x3) and (x4,m,x2,x3), for
                % type 0 (x4,m,x3,x2), both of type mod(t+1,3); labels of
                % the opposite faces follow the vertices, 0 for the cut
                X = nodes4e(mark,:);
                L = label(mark,:);
                t = type(mark);
                t0 = t == 0;
                nM = size(X,1);
                child1 = [X(:,1),m,X(:,2),X(:,3)];
                label1 = [zeros(nM,1),L(:,4),L(:,2),L(:,3)];
                child2 = [X(:,4),m,X(:,2),X(:,3)];
                label2 = [zeros(nM,1),L(:,1),L(:,2),L(:,3)];
                child2(t0,3:4) = X(t0,[3 2]);
                label2(t0,3:4) = L(t0,[3 2]);
                % initial types M (3) and O (4) of the longest edge labels:
                % (x1,m,x2,x3) resp. (x2,m,x1,x3) for O and (x2,m,x4,x3),
                % both of type 1
                tMO = t >= 3;
                tO = t == 4;
                child1(tO,:) = [X(tO,2),m(tO),X(tO,1),X(tO,3)];
                label1(tO,:) = [L(tO,2),L(tO,4),zeros(nnz(tO),1),L(tO,3)];
                child2(tMO,:) = [X(tMO,2),m(tMO),X(tMO,4),X(tMO,3)];
                label2(tMO,:) = [L(tMO,2),L(tMO,1),zeros(nnz(tMO),1),L(tMO,3)];
                tChild = mod(t+1,3);
                tChild(tMO) = 1;
                nodes4e = [nodes4e(~mark,:);child1;child2];
                label = [label(~mark,:);label1;label2];
                type = [type(~mark);tChild;tChild];
                parent = [parent(~mark);parent(mark);parent(mark)];
            end
            if any(mark)
                error('fem4m:Mesh3d:closureDiverged', ...
                    'The bisection closure does not terminate: inconsistent t4e labels.');
            end
            % positive orientation, boundary lists from the labels
            negative = Mesh3d.signedVolumes(c4n,nodes4e) < 0;
            nodes4e(negative,[2 3]) = nodes4e(negative,[3 2]);
            label(negative,[2 3]) = label(negative,[3 2]);
            others = [2 3 4;1 3 4;1 2 4;1 2 3];
            n4b = zeros(0,3);
            code = zeros(0,1);
            for k = 1 : 4
                onBd = label(:,k) > 0;
                n4b = [n4b;nodes4e(onBd,others(k,:))]; %#ok<AGROW>
                code = [code;label(onBd,k)]; %#ok<AGROW>
            end
            nodes4Sb = [n4b(code >= 3,:),code(code >= 3) - 3];
            obj = Mesh3d(c4n,nodes4e,n4b(code == 1,:),n4b(code == 2,:),nodes4Sb, ...
                't4e',Mesh3d.encodeLabels(type,negative));
            obj.parent4e = parent;
            obj.prolongation = Pt';
        end

        function [obj,indDb,indNb,indSb] = fixBoundaryOrientation(obj)
            %fixBoundaryOrientation - Orient boundary faces outward
            %   Flips the boundary triples whose normal points into the
            %   adjacent element, together with their edge ids.
            [obj.n4Db,obj.ed4Db,indDb] = orientOutward(obj,obj.n4Db,obj.f4Db,obj.ed4Db);
            [obj.n4Nb,obj.ed4Nb,indNb] = orientOutward(obj,obj.n4Nb,obj.f4Nb,obj.ed4Nb);
            [obj.n4Sb,obj.ed4Sb,indSb] = orientOutward(obj,obj.n4Sb,obj.f4Sb,obj.ed4Sb);
        end

        function [n4b,ed4b,ind] = orientOutward(obj,n4b,f4b,ed4b)
            %orientOutward - Flip boundary triples that point inward
            ind = false(size(n4b,1),1);
            if isempty(n4b)
                return
            end
            elems = obj.e4f(f4b,1);
            % the node of the adjacent element opposite to the face
            opposite = sum(double(obj.n4e(elems,:)),2) - sum(double(n4b(:,1:3)),2);
            c1 = obj.c4n(n4b(:,1),:);
            normal = cross(obj.c4n(n4b(:,2),:)-c1,obj.c4n(n4b(:,3),:)-c1,2);
            ind = dot(normal,obj.c4n(opposite,:)-c1,2) > 0;
            n4b(ind,[2 3]) = n4b(ind,[3 2]);
            ed4b(ind,:) = ed4b(ind,[3 2 1]);
        end

        function obj = enumerate(obj)
            %ENUMERATE - Build faces, edges and topological numbers
            %   Faces and edges are numbered in lexicographic order of
            %   their sorted node ids. Every face belongs to two elements
            %   or to one element and exactly one boundary triple, and
            %   ed4f(:,j) is the edge from n4f(:,j) to n4f(:,mod(j,3)+1).
            localFaces = [1 2 3;1 2 4;4 2 3;1 3 4];
            localEdges = [1 2;2 3;3 1;1 4;2 4;3 4];
            edges4face = [1 2 3;1 5 4;2 6 5;3 4 6];   % local edges of the local faces
            opposite4face = [3 1 2;4 1 2;4 2 3;4 3 1]; % local node opposite to these edges
            nE = size(obj.n4e,1);
            nN = size(obj.c4n,1);
            nodes4e = uint32(obj.n4e);
            n4b = {obj.n4Db,obj.n4Nb,obj.n4Sb};
            nB = cellfun(@(x)size(x,1),n4b);
            n4Bd = zeros(sum(nB),3,'uint32');
            offset = 0;
            for k = 1 : numel(n4b)
                if nB(k) > 0
                    n4Bd(offset+(1:nB(k)),:) = n4b{k}(:,1:3);
                end
                offset = offset + nB(k);
            end
            % faces, element-major: position 4*(e-1)+k is local face k of
            % element e, so the first occurrence is the smaller element
            faces = [Mesh3d.localSides(nodes4e,localFaces);sort(n4Bd,2)];
            [p,st] = Mesh3d.sortedGroups(faces,nN);
            count = diff([find(st);numel(st)+1]);
            if any(count > 2)
                error('fem4m:Mesh3d:duplicateFace', ...
                    'A face belongs to more than two elements and boundary triples.');
            end
            [f4i,first,last] = Mesh3d.numberGroups(p,st);
            if any(first > 4*nE)
                error('fem4m:Mesh3d:unknownBoundaryFace', ...
                    'A boundary triple is not a face of the mesh.');
            end
            if any(count == 1)
                error('fem4m:Mesh3d:unmatchedFaces', ...
                    ['Faces with one element are in no boundary list: incomplete boundary ', ...
                    'lists, or a non-conforming refinement (inconsistent t4e labels).']);
            end
            obj.n4f = faces(first,:);
            obj.f4e = reshape(f4i(1:4*nE),4,nE)';
            offset = 4*nE;
            obj.f4Db = f4i(offset+(1:nB(1)));
            obj.f4Nb = f4i(offset+nB(1)+(1:nB(2)));
            obj.f4Sb = f4i(offset+nB(1)+nB(2)+(1:nB(3)));
            last(last == first | last > 4*nE) = 0;
            elem1 = ceil(first/4);
            face1 = first - 4*(elem1-1);
            elem2 = ceil(last/4);
            face2 = last - 4*(elem2-1);
            face2(last == 0) = 0;
            obj.e4f = uint32([elem1,elem2]);
            obj.fNr4f = uint32([face1,face2]);
            clear faces p st count f4i last
            % edges, element-major as the faces
            edges = [Mesh3d.localSides(nodes4e,localEdges);
                Mesh3d.localSides(n4Bd,[1 2;2 3;3 1])];
            [p,st] = Mesh3d.sortedGroups(edges,nN);
            [ed4i,first] = Mesh3d.numberGroups(p,st);
            obj.n4ed = edges(first,:);
            obj.ed4e = reshape(ed4i(1:6*nE),6,nE)';
            ed4Bd = reshape(ed4i(6*nE+1:end),3,[])';
            obj.ed4Db = ed4Bd(1:nB(1),:);
            obj.ed4Nb = ed4Bd(nB(1)+(1:nB(2)),:);
            obj.ed4Sb = ed4Bd(nB(1)+nB(2)+(1:nB(3)),:);
            clear edges p st ed4i
            % ed4f relative to the sorted node triples of n4f
            nF = size(obj.n4f,1);
            obj.ed4f = zeros(nF,3,'uint32');
            for j = 1 : 3
                oppNode = reshape(nodes4e(sub2ind([nE,4],elem1,opposite4face(face1,j))),[],1);
                pos = 1*(oppNode == obj.n4f(:,1)) + 2*(oppNode == obj.n4f(:,2)) ...
                    + 3*(oppNode == obj.n4f(:,3));
                slot = mod(pos,3) + 1;
                obj.ed4f(sub2ind([nF,3],(1:nF)',slot)) = ...
                    reshape(obj.ed4e(sub2ind([nE,6],elem1,edges4face(face1,j))),[],1);
            end
            % numbers
            obj.nrNodes = nN;
            obj.nrElems = nE;
            obj.nrEdges = size(obj.n4ed,1);
            obj.nrFaces = nF;
        end

    end

    methods (Static, Access=private)

        function sides = localSides(n4x,local)
            %localSides - Sorted node tuples of local sides, element-major
            %   SIDES = localSides(N4X,LOCAL) returns the node ids of the
            %   local sides LOCAL (one row per side) of every row of N4X,
            %   sorted within each tuple; row K*(e-1)+k is local side k of
            %   row e, K = size(LOCAL,1).
            [nX,~] = size(n4x);
            [K,c] = size(local);
            sides = zeros(K*nX,c,'uint32');
            for k = 1 : K
                sides(k:K:end,:) = n4x(:,local(k,:));
            end
            sides = sort(sides,2);
        end

        function [x4i,first,last] = numberGroups(p,st)
            %numberGroups - Number groups of equal tuples in sorted order
            %   [X4I,FIRST,LAST] = numberGroups(P,ST) numbers the groups of
            %   the stably sorted tuples (sortedGroups) in lexicographic
            %   order. X4I is the group id of every original row, FIRST
            %   and LAST the first and last row of every group.
            first = double(p(st));
            last = double(p([st(2:end);true]));
            x4i = zeros(numel(p),1,'uint32');
            x4i(p) = cumsum(uint32(st));
        end

        function [p,st] = sortedGroups(rows,nrNodes)
            %sortedGroups - Stable sort of node tuples into groups
            %   [P,ST] = sortedGroups(ROWS,NRNODES) sorts the rows of the
            %   uint32 tuples ROWS stably in lexicographic order, ROWS(P,:),
            %   and marks with ST the first row of each group of equal
            %   tuples. The rows are combined into one uint64 key unless
            %   that overflows for NRNODES nodes.
            n = size(rows,1);
            if n == 0
                p = zeros(0,1);
                st = false(0,1);
                return
            end
            N1 = uint64(nrNodes) + 1;
            if double(N1)^size(rows,2) < 2^64
                key = uint64(rows(:,1));
                for j = 2 : size(rows,2)
                    key = key*N1 + uint64(rows(:,j));
                end
                [key,p] = sort(key);
                st = [true;key(2:end) ~= key(1:end-1)];
            else
                [rows,p] = sortrows(rows);
                st = [true;any(rows(2:end,:) ~= rows(1:end-1,:),2)];
            end
        end

        function keys = edgeKeys(pairs,K)
            %edgeKeys - Keys min*K+max of node pairs
            keys = uint64(min(pairs,[],2))*K + uint64(max(pairs,[],2));
        end

        function vol = signedVolumes(c4n,n4e)
            %signedVolumes - Six times the signed volumes of tetrahedra
            n4e = double(n4e);
            x1 = c4n(n4e(:,1),:);
            vol = dot(cross(c4n(n4e(:,2),:)-x1,c4n(n4e(:,3),:)-x1,2),c4n(n4e(:,4),:)-x1,2);
        end

        function [n4e,t4e] = colorLabels(n4e,nrNodes)
            %colorLabels - Bisection labels from a greedy vertex coloring
            %   The vertices are colored greedily in the order of their
            %   ids (as color4n); every element gets its vertices sorted
            %   by the key mod(color,N), N the largest color (the largest
            %   color first), and type 0 (Diening, Gehring, Storn).
            nE = size(n4e,1);
            edges = unique(Mesh3d.localSides(uint32(n4e),[1 2;1 3;1 4;2 3;2 4;3 4]),'rows');
            edges = double(edges);
            adjacency = sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],1,nrNodes,nrNodes);
            color = zeros(nrNodes,1);
            for v = 1 : nrNodes
                used = color(adjacency(:,v) ~= 0);
                color(v) = find(~ismember(1:numel(used)+1,used),1);
            end
            key = mod(color,max(color));         % colors 1..N: largest -> 0
            [~,order] = sort(reshape(key(n4e),nE,4),2);
            n4e = n4e(sub2ind(size(n4e),repmat((1:nE)',1,4),order));
            t4e = zeros(nE,1,'int8');
        end

        function [n4e,t4e] = longestEdgeLabels(c4n,n4e)
            %longestEdgeLabels - Bisection labels from the longest edges
            %   The edges are ranked by their computed length, equal values
            %   by their sorted node ids; every element is bisected at its
            %   edge of largest rank ab and every face marks its edge of
            %   largest rank (Arnold, Mukherjee, Pouly). With the marked
            %   edges of the faces acd and bcd the element is of type P, A,
            %   M or O: P (both marks meet at p) gives (a,q,p,b) of type 1,
            %   A (marks a-p and b-q) (a,q,p,b) of type 0, M (acd marks
            %   a-p, bcd marks cd) (a,q,p,b) of type M, with a and b
            %   swapped where bcd marks b-p and acd marks cd, and O (both
            %   mark cd) (a,c,d,b) of type O, with {p,q} = {c,d}.
            nE = size(n4e,1);
            n4e = double(n4e);
            local = [1 2;1 3;1 4;2 3;2 4;3 4];
            edges = Mesh3d.localSides(uint32(n4e),local);
            [p,st] = Mesh3d.sortedGroups(edges,size(c4n,1));
            ed4i = Mesh3d.numberGroups(p,st);
            n4ed = double(edges(p(st),:));
            len = vecnorm(c4n(n4ed(:,2),:) - c4n(n4ed(:,1),:),2,2);
            [~,order] = sortrows([len,n4ed]);
            rank4ed = zeros(numel(order),1);
            rank4ed(order) = 1 : numel(order);
            rank4e = reshape(rank4ed(ed4i),6,nE)';
            % local refinement edge ab and opposite edge cd
            [~,jRef] = max(rank4e,[],2);
            others = [3 4;2 4;2 3;1 4;1 3;1 2];
            a = local(jRef,1);
            b = local(jRef,2);
            c = others(jRef,1);
            d = others(jRef,2);
            pairId = [0 1 2 3;1 0 4 5;2 4 0 6;3 5 6 0];
            rk = @(u,v)rank4e(sub2ind([nE,6],(1:nE)',pairId(sub2ind([4 4],u,v))));
            pa = Mesh3d.markedVertex(rk(a,c),rk(a,d),rk(c,d),c,d);
            pb = Mesh3d.markedVertex(rk(b,c),rk(b,d),rk(c,d),c,d);
            q = c + d - pa;
            order4e = [a,q,pa,b];
            kind = ones(nE,1);
            % A: the marks a-pa and b-pb do not meet
            kind(pa > 0 & pb > 0 & pa ~= pb) = 0;
            % M: one face marks cd, the other is the face (a,q,p)
            I = pa == 0 & pb > 0;
            order4e(I,:) = [b(I),c(I)+d(I)-pb(I),pb(I),a(I)];
            kind(xor(pa == 0,pb == 0)) = 3;
            % O: both faces mark cd
            I = pa == 0 & pb == 0;
            order4e(I,:) = [a(I),c(I),d(I),b(I)];
            kind(I) = 4;
            n4e = uint32(n4e(sub2ind([nE,4],repmat((1:nE)',1,4),order4e)));
            t4e = Mesh3d.encodeLabels(kind,false(nE,1));
        end

        function p = markedVertex(rUc,rUd,rCd,c,d)
            %markedVertex - Vertex of the marked edge of the face (u,c,d)
            %   P = markedVertex(RUC,RUD,RCD,C,D) returns C or D where the
            %   edge u-c resp. u-d has the largest rank and 0 for c-d.
            [~,j] = max([rUc,rUd,rCd],[],2);
            p = zeros(size(c));
            p(j == 1) = c(j == 1);
            p(j == 2) = d(j == 2);
        end

        function t4e = encodeLabels(kind,flip)
            %encodeLabels - Labels from the bisection type and the flip
            %   T4E = encodeLabels(KIND,FLIP) is KIND + 3*FLIP for the
            %   types 0..2 and 6 + 2*FLIP resp. 7 + 2*FLIP for the
            %   initial types M (3) and O (4).
            kind = double(kind(:));
            flip = double(flip(:) ~= 0);
            t4e = kind + 3*flip;
            mo = kind >= 3;
            t4e(mo) = kind(mo) + 3 + 2*flip(mo);
            t4e = int8(t4e);
        end

        function [kind,flip] = decodeLabels(t4e)
            %decodeLabels - Bisection type and flip of the labels
            %   [KIND,FLIP] = decodeLabels(T4E) inverts encodeLabels.
            t4e = double(t4e(:));
            mo = t4e >= 6;
            kind = mod(t4e,3);
            kind(mo) = mod(t4e(mo),2) + 3;
            flip = (t4e >= 3 & ~mo) | t4e >= 8;
        end

        function n4Sb = checkN4Sb(n4Sb)
            %checkN4Sb - Symmetry boundary triples with their axis
            if isempty(n4Sb)
                n4Sb = zeros(0,4,'uint32');
                return
            end
            switch size(n4Sb,2)
                case 4
                case 1
                    error('fem4m:Mesh3d:invalidN4Sb', ...
                        'The 5th argument is n4Sb; pass bisection labels as Mesh3d(...,''t4e'',T4E).');
                otherwise
                    error('fem4m:Mesh3d:invalidN4Sb', ...
                        'n4Sb must have the rows [n1 n2 n3 axis].');
            end
            n4Sb = uint32(n4Sb);
        end

        function list = loadList(file)
            %loadList - Load a boundary list file if it exists
            if exist(file,'file')
                list = uint32(load(file));
            else
                list = [];
            end
        end

        function pkInterpolationRemoved()
            %pkInterpolationRemoved - Error for the removed prolongation object
            error('fem4m:Mesh3d:pkInterpolationRemoved', ...
                ['Meshes take no prolongation object anymore: use ', ...
                'P = fe.prolongation(feOld) after the refinement.']);
        end

    end

end