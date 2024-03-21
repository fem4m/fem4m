classdef Mesh
    %MESH - Triangular mesh with adaptive refinement and coarsening
    %   OBJ = MESH(NAME) loads the mesh from the files NAME_c4n.dat and
    %   NAME_n4e.dat and the boundary files NAME_n4Db.dat, NAME_n4Nb.dat
    %   and NAME_n4Sb.dat where they exist (save writes these files).
    %   NAME is a char vector.
    %
    %   OBJ = MESH(GEOM) creates the mesh of the Geometry object GEOM, or
    %   of Geometry(GEOM) for an options struct GEOM, including the
    %   projection of new boundary nodes onto a curved boundary where
    %   the geometry defines one (GEOM.projectBoundary).
    %
    %   OBJ = MESH(C4N,N4E) creates the mesh with the node coordinates
    %   C4N (one row per node) and the node ids N4E of the triangles (one
    %   row per element); all boundary edges are Dirichlet edges.
    %
    %   OBJ = MESH(C4N,N4E,n4Db) also specifies the node pairs n4Db of
    %   the Dirichlet boundary edges, which must cover the boundary; there
    %   is no Neumann and no symmetry boundary.
    %
    %   OBJ = MESH(C4N,N4E,n4Db,n4Nb) also specifies the node pairs n4Nb
    %   of the Neumann boundary edges. n4Db, n4Nb and n4Sb together must
    %   cover the boundary, every boundary edge with the orientation it
    %   has in its element.
    %
    %   OBJ = MESH(...,n4Sb) also specifies the symmetry boundary: node
    %   pairs plus the fixed displacement component (1 or 2, 0 for a
    %   general normal) in the third column.
    %
    %   Row y of every x4y lists the ids of the x of y: e4n is the sparse
    %   nrNodes-by-nrElems incidence with e4n(n,e) = e for the nodes n of
    %   the element e, e4ed(ed,:) are the two elements of the edge ed (0
    %   on the boundary), and ed4n is the symmetric node-pair lookup
    %   ed4n(i,j) of the edge between the nodes i and j.
    %
    %   Every refinement stores the finite element independent data of
    %   the prolongation from the mesh before the refine call: parent4e
    %   and the P1 prolongation matrix. Finite element prolongations are
    %   computed from them by FE.prolongation. Unless coarse returns the
    %   mesh unchanged, it sets both to [] and stores children4e instead,
    %   with one row per element before the coarsening: [e,f] if the
    %   element e absorbed its sibling f, [e,0] if it was kept, and zeros
    %   for the absorbed sibling; the nonzero rows belong to the new
    %   elements in order.
    %
    %   MESH functions:
    %       area4e          - Compute the areas of elements
    %       area4n          - Compute nodal patch areas
    %       area4ed         - Compute edge patch areas
    %       length4ed       - Compute the lengths of edges
    %       point4e         - Compute points in elements
    %       point4ed        - Compute points on edges
    %       normal4e        - Compute element edge normals
    %       normal4ed       - Compute the unit edge normals
    %       tangent4e       - Compute element edge tangents
    %       tangent4ed      - Compute the unit edge tangents
    %       color4n         - Color nodes of adjacent elements
    %       color4e         - Color adjacent elements
    %       patch4e         - Collect adjacent elements
    %       patch4n         - Collect elements per node
    %       transformation  - Compute inverse mapping matrix
    %       jacobian        - Compute the mapping Jacobian
    %       refine          - Refine the mesh
    %       coarse          - Coarsen marked element pairs
    %       disp            - Display summary or full data of the mesh
    %       plot            - Plot the triangular mesh
    %       plotBd          - Plot the boundary of the mesh
    %       plotEnumeration - Label nodes, elements, and edges with ids
    %       plotColor4n     - Plot mesh with node colors
    %       plotColor4e     - Plot mesh with element colors
    %       removeNode      - Remove nodes and their elements
    %       removeEdge      - Remove an edge and its elements
    %       addTriangle     - Add a triangle to the mesh
    %       save            - Write mesh data to .dat files
    %       closure         - Complete marked edges for RGB refinement
    %
    %   MESH properties:
    %       c4n          - Coordinates of nodes
    %       n4e          - Node ids of the elements
    %       n4Db         - Node pairs of Dirichlet boundary
    %       n4Nb         - Node pairs of Neumann boundary
    %       n4Sb         - Node pairs of symmetry boundary plus fixed axis
    %       n4ed         - Node pairs of unique edges
    %       ed4n         - Sparse matrix of edge ids for node pairs
    %       ed4e         - Edge ids of the three edges of each element
    %       e4ed         - Adjacent element ids for each edge
    %       e4n          - Sparse matrix of element ids for nodes
    %       ed4Db        - Edge ids of Dirichlet boundary
    %       ed4Nb        - Edge ids of Neumann boundary
    %       ed4Sb        - Edge ids of symmetry boundary
    %       edNr4ed      - Local edge numbers of edges within elements
    %       nNr4n        - Sparse matrix of local node numbers
    %       NbEd4n       - Neumann edge ids incident to each node
    %       SbEd4n       - Symmetry edge ids incident to each node
    %       DbEd4n       - Dirichlet edge ids incident to each node
    %       nrNodes      - Number of nodes
    %       nrEdges      - Number of edges
    %       nrElems      - Number of elements
    %       markedEdges  - 0/1 marks of edges to be refined
    %       refineEdges  - Closure hull of the marked edges
    %       prolongation - P1 prolongation from the mesh before refine
    %       parent4e     - Ancestor element before the last refine call
    %       children4e   - Element pairs merged by the last coarse call
    %
    %   See also MeshSurface, Mesh3d, MeshPoly, Geometry

    % Copyright (c) 2026 by the fem4m authors

    properties
        c4n         % coordinates of nodes
        n4e         % node ids of the elements
        n4Db        % node pairs of Dirichlet boundary
        n4Nb        % node pairs of Neumann boundary
        n4Sb        % node pairs of symmetry boundary plus fixed axis
        n4ed        % node pairs of unique edges
        ed4n        % sparse matrix of edge ids for node pairs
        ed4e        % edge ids of the three edges of each element
        e4ed        % adjacent element ids for each edge
        e4n         % sparse matrix of element ids for nodes
        ed4Db       % edge ids of Dirichlet boundary
        ed4Nb       % edge ids of Neumann boundary
        ed4Sb       % edge ids of symmetry boundary
        edNr4ed     % local edge numbers of edges within elements
        nNr4n       % sparse matrix of local node numbers
        NbEd4n      % Neumann edge ids incident to each node
        SbEd4n      % symmetry edge ids incident to each node
        DbEd4n      % Dirichlet edge ids incident to each node
        nrNodes     % number of nodes
        nrEdges     % number of edges
        nrElems     % number of elements
        markedEdges % 0/1 marks of edges to be refined
        refineEdges % closure hull of the marked edges
        prolongation % P1 prolongation from the mesh before refine
        parent4e    % ancestor element before the last refine call
        children4e  % element pairs merged by the last coarse call
    end

    properties (Dependent, Hidden)
        pkInterpolation % removed: use FE.prolongation
    end

    properties (Access = private)
        N0
        uniformC4n
        projectBoundary
    end

    methods
        function obj = Mesh(c4n,n4e,n4Db,n4Nb,n4Sb,varargin)
            %MESH - Create a triangular mesh object
            if ~isempty(varargin) || (nargin == 2 && (ischar(c4n) || isstruct(c4n) ...
                    || isa(c4n,'Geometry')))
                error('fem4m:Mesh:pkInterpolationRemoved', ...
                    ['Meshes take no prolongation object anymore: use ', ...
                    'P = fe.prolongation(feOld) after the refinement.']);
            end
            if ischar(c4n)
                name = c4n;
                obj.c4n   = load([name,'_c4n.dat']);
                obj.n4e   = load([name,'_n4e.dat']);
                if exist([name,'_n4Db.dat'],'file')
                    obj.n4Db = load([name,'_n4Db.dat']);
                else
                    obj.n4Db = [];
                end
                if exist([name,'_n4Nb.dat'],'file')
                    obj.n4Nb  = load([name,'_n4Nb.dat']);
                else
                    obj.n4Nb = [];
                end
                if exist([name,'_n4Sb.dat'],'file')
                    obj.n4Sb  = load([name,'_n4Sb.dat']);
                else
                    obj.n4Sb = [];
                end
            elseif isstruct(c4n) || strcmpi(class(c4n),'Geometry')
                geom = c4n;
                if isstruct(c4n)
                    geom = Geometry(geom);
                end
                obj.c4n = geom.c4n;
                obj.n4e = geom.n4e;
                obj.n4Db = geom.n4Db;
                obj.n4Nb = geom.n4Nb;
                obj.n4Sb = geom.n4Sb;
                obj.projectBoundary = geom.projectBoundary;
            elseif nargin > 1
                obj.c4n = c4n;
                obj.n4e = n4e;
                if nargin < 3
                    % all boundary edges (edges of one element) are Dirichlet
                    edges = [n4e(:,[1 2]);n4e(:,[2 3]);n4e(:,[3 1])];
                    [~,~,j] = unique(sort(edges,2),'rows');
                    once = accumarray(j,1) == 1;
                    n4Db = edges(once(j),:);
                end
                obj.n4Db = n4Db;
                if nargin > 3
                    obj.n4Nb = n4Nb;
                else
                    obj.n4Nb = [];
                end
                if nargin > 4
                    obj.n4Sb = n4Sb;
                else
                    obj.n4Sb = [];
                end
            end
            if ~isempty(obj.n4Sb) && size(obj.n4Sb,2) ~= 3
                error('fem4m:Mesh:invalidN4Sb', ...
                    'n4Sb must have the rows [n1 n2 axis].');
            end
            obj.N0 = size(obj.c4n,1);
            obj = getDataStructures(obj);
        end

        function val = get.pkInterpolation(~) %#ok<STOUT>
            %pkInterpolation - Removed, use FE.prolongation
            error('fem4m:Mesh:pkInterpolationRemoved', ...
                ['Meshes hold no prolongation object anymore: use ', ...
                'P = fe.prolongation(feOld) after the refinement.']);
        end

        function obj = set.pkInterpolation(obj,~) %#ok<INUSL>
            %pkInterpolation - Removed, use FE.prolongation
            error('fem4m:Mesh:pkInterpolationRemoved', ...
                ['Meshes hold no prolongation object anymore: use ', ...
                'P = fe.prolongation(feOld) after the refinement.']);
        end

        function area4e = area4e(obj)
            %area4e - Compute the areas of elements
            %   AREA4E = area4e(OBJ) computes the signed area of every
            %   triangle from its node coordinates, positive for a
            %   counterclockwise node order (nrElems-by-1).
            x1 = obj.c4n(obj.n4e(:,1),1);
            x2 = obj.c4n(obj.n4e(:,2),1);
            x3 = obj.c4n(obj.n4e(:,3),1);
            y1 = obj.c4n(obj.n4e(:,1),2);
            y2 = obj.c4n(obj.n4e(:,2),2);
            y3 = obj.c4n(obj.n4e(:,3),2);
            area4e = ( x1.*(y2 - y3) + x2.*(y3 - y1) + x3.*(y1 - y2) )/2;
        end

        function area4n = area4n(obj)
            %area4n - Compute nodal patch areas
            %   AREA4N = area4n(OBJ) computes for every node the area of
            %   its patch, the sum of the areas of all elements containing
            %   it (nrNodes-by-1).
            area = obj.area4e*ones(1,3);
            area4n = accumarray(obj.n4e(:),area(:),[size(obj.c4n,1),1]);
        end

        function area4n = area4ed(obj)
            %area4ed - Compute edge patch areas
            %   AREA4ED = area4ed(OBJ) computes for every edge the area of
            %   its patch, the sum of the areas of its one or two adjacent
            %   elements (nrEdges-by-1).
            area = obj.area4e*ones(1,3);
            area4n = accumarray(obj.ed4e(:),area(:),[obj.nrEdges,1]);
        end

        function length4ed = length4ed(obj)
            %length4ed - Compute the lengths of edges
            %   LENGTH4ED = length4ed(OBJ) computes the Euclidean length
            %   of every edge listed in OBJ.n4ed (nrEdges-by-1).
            length4ed = sqrt(sum((obj.c4n(obj.n4ed(:,2),:)-obj.c4n(obj.n4ed(:,1),:)).^2,2));
        end

        function val = point4e(obj,refp)
            %point4e - Compute points in elements
            %   VAL = point4e(OBJ,REFP) computes for every element the
            %   image P1 + (P2-P1)*REFP(1) + (P3-P1)*REFP(2) of the
            %   reference point REFP (1-by-2), where P1, P2, P3 are the
            %   nodes of the triangle; VAL has one row per element.
            P1 = obj.c4n(obj.n4e(:,1),:);
            P2 = obj.c4n(obj.n4e(:,2),:);
            P3 = obj.c4n(obj.n4e(:,3),:);
            val = P1 + (P2-P1)*refp(1) + (P3-P1)*refp(2);
        end

        function val = point4ed(obj,refp)
            %point4ed - Compute points on edges
            %   VAL = point4ed(OBJ,REFP) computes for every edge the
            %   point with the scalar reference coordinate REFP, where 0
            %   and 1 are its first and second node in OBJ.n4ed; VAL has
            %   one row per edge.
            P1 = obj.c4n(obj.n4ed(:,1),:);
            P2 = obj.c4n(obj.n4ed(:,2),:);
            val = P1 + (P2-P1)*refp(1);
        end

        function normal4e = normal4e(obj)
            %normal4e - Compute element edge normals
            %   NORMAL4E = normal4e(OBJ) computes for every element the
            %   unit normal vectors of its three edges by rotating their
            %   unit tangents (see tangent4e) clockwise by a quarter turn,
            %   which gives the outer normals for a counterclockwise node
            %   order. NORMAL4E is nrElems-by-3-by-2.
            tangents = obj.tangent4e;
            normal4e = zeros(obj.nrElems,3,2);
            normal4e(:,1,:) = [tangents(:,1,2), -tangents(:,1,1)];
            normal4e(:,2,:) = [tangents(:,2,2), -tangents(:,2,1)];
            normal4e(:,3,:) = [tangents(:,3,2), -tangents(:,3,1)];
        end

        function normal4ed = normal4ed(obj)
            %normal4ed - Compute the unit edge normals
            %   NORMAL4ED = normal4ed(OBJ) computes the unit normal vector
            %   of every edge by rotating its unit tangent (see
            %   tangent4ed) clockwise by a quarter turn (nrEdges-by-2).
            %   Boundary edges keep the orientation they have in their
            %   element, so for counterclockwise elements their normals
            %   are the outer unit normals.
            tangents = obj.tangent4ed;
            normal4ed = [tangents(:,2), -tangents(:,1)];
        end

        function tangent4e = tangent4e(obj)
            %tangent4e - Compute element edge tangents
            %   TANGENT4E = tangent4e(OBJ) computes for every element the
            %   unit tangent vectors of its three edges, oriented along
            %   the node order of the triangle: TANGENT4E(e,k,:) points
            %   from the node k to the node mod(k,3)+1 of the element e
            %   (nrElems-by-3-by-d).
            allEdges = [obj.n4e(:,[1 2]); obj.n4e(:,[2 3]); obj.n4e(:,[3 1])];
            c4start  = obj.c4n(allEdges(:,1),:);
            c4end    = obj.c4n(allEdges(:,2),:);
            lengths  = sqrt(sum((c4end-c4start).^2,2));
            tangents = (c4end - c4start)./lengths;
            tangent4e = zeros(obj.nrElems,3,size(tangents,2));
            tangent4e(:,1,:) = tangents(1:size(obj.n4e,1),:);
            tangent4e(:,2,:) = tangents(size(obj.n4e,1)+1:2*size(obj.n4e,1),:);
            tangent4e(:,3,:) = tangents(2*size(obj.n4e,1)+1:3*size(obj.n4e,1),:);
        end

        function tangent4ed = tangent4ed(obj)
            %tangent4ed - Compute the unit edge tangents
            %   TANGENT4ED = tangent4ed(OBJ) computes the unit tangent
            %   vector of every edge pointing from its first to its
            %   second node in OBJ.n4ed (nrEdges-by-d).
            c4start = obj.c4n(obj.n4ed(:,1),:);
            c4end = obj.c4n(obj.n4ed(:,2),:);
            tangent4ed = (c4end-c4start)./obj.length4ed;
        end

        function color4n = color4n(obj,levels)
            %color4n - Color nodes of adjacent elements
            %   COLOR4N = color4n(OBJ) colors every node greedily with a
            %   small positive integer such that the two nodes of every
            %   edge receive different colors (nrNodes-by-1).
            %
            %   COLOR4N = color4n(OBJ,LEVELS) gives different colors to
            %   any two nodes connected by a path of at most LEVELS edges
            %   (default 1).
            if nargin < 2 || isempty(levels)
                levels = 1;
            end
            color4n = zeros(obj.nrNodes,1);
            for curNode = 1 : obj.nrNodes
                neighbours = curNode;
                for cuLevel = 1 : levels
                    [row,~] = find(obj.ed4n(:,neighbours));
                    neighbours = [neighbours;row]; %#ok<AGROW>
                end
                neighbour_colors = color4n(neighbours);
                color4n(curNode) = min(setdiff(1:max(neighbour_colors)+1,neighbour_colors));
            end
        end

        function color4e = color4e(obj,levels)
            %color4e - Color adjacent elements
            %   COLOR4E = color4e(OBJ) colors every element greedily with
            %   a small positive integer such that two elements receive
            %   different colors if they share a node or both share a
            %   node with a third element (nrElems-by-1).
            %
            %   COLOR4E = color4e(OBJ,LEVELS) gives different colors to
            %   any two elements connected by a chain of at most LEVELS+2
            %   elements, each sharing a node with the next (default
            %   LEVELS = 1; LEVELS = 0 separates only elements that share
            %   a node).
            if nargin < 2 || isempty(levels)
                levels = 1;
            end
            color4e = zeros(obj.nrElems,1);
            e4nT = obj.e4n';
            for curElem = 1 : obj.nrElems
                nodes = obj.n4e(curElem,:);
                [~,~,neighbours] = find(e4nT(:,nodes));
                neighbours = unique(neighbours);
                for cuLevel = 1 : levels
                    nodes = unique(obj.n4e(neighbours,:));
                    [~,~,v] = find(e4nT(:,nodes));
                    neighbours = [neighbours;unique(v(:))]; %#ok<AGROW>
                end
                neighbour_colors = color4e(neighbours);
                color4e(curElem) = min(setdiff(1:max(neighbour_colors)+1,neighbour_colors));
            end
        end

        function patch4e = patch4e(obj)
            %patch4e - Collect adjacent elements
            %   PATCH4E = patch4e(OBJ) returns for every element the ids
            %   of all elements sharing at least one node with it, the
            %   element itself included, as an nrElems-by-1 cell array of
            %   sorted column vectors.
            patch4e = cell(obj.nrElems,1);
            e4nT = obj.e4n';
            for curElem = 1 : obj.nrElems
                [~,~,v] = find(e4nT(:,obj.n4e(curElem,:)));
                patch4e{curElem} = unique(v);
            end
        end

        function patch4n = patch4n(obj)
            %patch4n - Collect elements per node
            %   PATCH4N = patch4n(OBJ) returns for every node the ids of
            %   all elements containing it, as an nrNodes-by-1 cell array
            %   of sorted column vectors.
            patch4n = cell(obj.nrNodes,1);
            e4nT = obj.e4n';
            for curNode = 1 : obj.nrNodes
                [~,~,v] = find(e4nT(:,curNode));
                patch4n{curNode} = unique(v);
            end
        end

        function val = transformation(obj,~)
            %TRANSFORMATION - Compute inverse mapping matrix
            %   VAL = TRANSFORMATION(OBJ) computes for every element the
            %   transposed inverse of the Jacobian B = [P2-P1, P3-P1] of
            %   the affine map from the reference triangle to the element
            %   with the nodes P1, P2, P3. VAL is nrElems-by-2-by-2 with
            %   VAL(e,:,:) = inv(B)' of the element e.
            %
            %   VAL = TRANSFORMATION(OBJ,REFP) returns the same; the
            %   reference point REFP is ignored since the map is affine
            %   (the elements pass it for curved meshes like MeshSurface).
            % x = F(x_ref) = Bx_ref + d
            % \nabla_x phi(F(x_ref)) = (B^(-1))' nabla_{x_ref}(phi(F(x_ref)))
            %                      = (B^(-1))' nabla_{x_ref} phi_ref(x_ref)
            % val = (B^(-1))' with B = [n2-n1, n3-n1]
            n1 = obj.c4n(obj.n4e(:,1),:);
            n2 = obj.c4n(obj.n4e(:,2),:);
            n3 = obj.c4n(obj.n4e(:,3),:);
            c21 = (0.5./obj.area4e).*(n2-n1);
            c31 = (0.5./obj.area4e).*(n3-n1);
            val = zeros(obj.nrElems,2,2);
            val(:,1,1) = c31(:,2);
            val(:,1,2) = -c21(:,2);
            val(:,2,1) = -c31(:,1);
            val(:,2,2) = c21(:,1);
        end

        function val = jacobian(obj,~)
            %JACOBIAN - Compute the mapping Jacobian
            %   VAL = JACOBIAN(OBJ) computes for every element the
            %   transposed Jacobian B' of the affine map from the
            %   reference triangle to the element, B = [P2-P1, P3-P1]
            %   with the nodes P1, P2, P3. VAL is nrElems-by-2-by-2 with
            %   the rows VAL(e,1,:) = P2-P1 and VAL(e,2,:) = P3-P1.
            %
            %   VAL = JACOBIAN(OBJ,REFP) returns the same; the reference
            %   point REFP is ignored since the map is affine.
            % x = F(x_ref) = Bx_ref + d
            % val = B' with B = [n2-n1, n3-n1]
            n1 = obj.c4n(obj.n4e(:,1),:);
            n2 = obj.c4n(obj.n4e(:,2),:);
            n3 = obj.c4n(obj.n4e(:,3),:);
            c21 = (n2-n1);
            c31 = (n3-n1);
            val = zeros(obj.nrElems,2,2);
            val(:,1,1) = c21(:,1);
            val(:,1,2) = c21(:,2);
            val(:,2,1) = c31(:,1);
            val(:,2,2) = c31(:,2);
        end

        function obj = refine(obj,varargin)
            %REFINE - Refine the mesh
            %   OBJ = REFINE(OBJ) refines the mesh along the marked edges
            %   OBJ.markedEdges by red-green-blue refinement, or uniformly
            %   with the warning fem4m:Mesh:noMarkedEdges when no edge is
            %   marked. The first edge n4e(:,[1 2]) of every element is
            %   its refinement edge.
            %
            %   OBJ = REFINE(OBJ,MARKED4ED) marks the edges given by the
            %   0/1 or logical vector MARKED4ED (one entry per edge; it
            %   replaces OBJ.markedEdges) before refining; a vector with
            %   one entry per element marks the three edges of every
            %   selected element instead.
            %
            %   OBJ = REFINE(OBJ,'uniform') refines all edges once.
            %
            %   OBJ = REFINE(OBJ,'uniform',N) repeats the uniform
            %   refinement N times; N = 0 keeps the mesh and sets
            %   parent4e and prolongation to the identity.
            %
            %   OBJ = REFINE(OBJ,'RGB') uses red-green-blue refinement
            %   (default), and OBJ = REFINE(OBJ,'NVB') uses newest vertex
            %   bisection instead. Both combine with MARKED4ED and
            %   'uniform', e.g. OBJ = REFINE(OBJ,'NVB','uniform').
            %
            %   OBJ = REFINE(OBJ,'graded',MU,H,GRADPOINTS,R) refines the
            %   first edges of the elements that are too large by
            %   red-green-blue refinement until every element T satisfies
            %   hT <= H and hT <= H*(rT/R(k))^(1-MU(k)) for every grading
            %   point GRADPOINTS(k,:), where hT is the longest edge of T
            %   and rT the distance of its centroid to the point.
            %   GRADPOINTS has one point per row (default [0,0]), R
            %   defaults to 1, and scalar MU and R hold for all points.
            %
            %   OBJ = REFINE(OBJ,'gradedbeta',MU,GRADPOINTS,R) builds a
            %   beta-graded mesh: it refines the mesh uniformly once and
            %   moves every node at a distance r < R(k) from the point
            %   GRADPOINTS(k,:) to the distance R(k)*(r/R(k))^(1/MU(k)) on
            %   the same ray. The unshifted coordinates are stored, so repeated
            %   'gradedbeta' calls without other refinements in between grade
            %   successive uniform refinements. GRADPOINTS defaults to
            %   [0,0] and R to 1, and scalar MU and R hold for all
            %   points. The meshes are not nested: parent4e and
            %   prolongation are those of the uniform refinement.
            %
            %   OBJ = REFINE(OBJ,'gradedred',KAPPA,GRADPOINTS) red-refines
            %   the mesh once; the new node of every edge at a grading
            %   point GRADPOINTS(k,:) (a node of the mesh, one per row)
            %   lies at the fraction KAPPA(k) of the edge from that point
            %   instead of at its midpoint (KAPPA = 1/2). A scalar KAPPA
            %   holds for all points. The error fem4m:Mesh:gradedRedArgs
            %   is raised without KAPPA or GRADPOINTS, and
            %   fem4m:Mesh:gradPointNotNode for a point that is no node.
            %
            %   The options are case-insensitive char vectors (other types
            %   are ignored); an unknown option raises the warning
            %   fem4m:Mesh:unknownRefineOption and is ignored, the removed
            %   option 'longestedge' the error fem4m:Mesh:longestEdgeRemoved.
            %   parent4e and the P1 prolongation always refer to the mesh
            %   before the call, also over several rounds.
            type = [];
            num = [];
            numVector = [];
            uniform = false;
            for k=1:length(varargin)
                if islogical(varargin{k})
                    % logical masks of marked edges or elements
                    numVector{end+1} = double(varargin{k}); %#ok<AGROW>
                elseif isa(varargin{k},'numeric')
                    if isscalar(varargin{k})
                        num{end+1} = varargin{k}; %#ok<AGROW>
                    else
                        numVector{end+1} = varargin{k}; %#ok<AGROW>
                    end
                elseif isa(varargin{k},'char')
                    if strcmpi(varargin{k},'UNIFORM')
                        uniform = true;
                    elseif strcmpi(varargin{k},'LONGESTEDGE')
                        error('fem4m:Mesh:longestEdgeRemoved', ...
                            ['The refine option ''longestedge'' was removed: order the ', ...
                            'vertices in n4e before creating the mesh (the refinement ', ...
                            'edge is n4e(:,[1 2])).']);
                    elseif strcmpi(varargin{k},'RGB') || strcmpi(varargin{k},'NVB') || strcmpi(varargin{k},'GRADED') || strcmpi(varargin{k},'GRADEDRED') || strcmpi(varargin{k},'GRADEDBETA')
                        type = varargin{k};
                    else
                        warning('fem4m:Mesh:unknownRefineOption','%s is no suitable refinement option (''UNIFORM'',''RGB'', ''NVB'', ''GRADED'', (int) number of uniform refinements)',varargin{k});
                    end
                end
            end
            if isempty(type); type = 'RGB'; end
            graded = any(strcmpi(type,{'GRADED','GRADEDBETA','GRADEDRED'}));
            if ~isempty(numVector) && ~graded
                if numel(numVector{1})==obj.nrElems
                    marked4e = numVector{1}(:);
                    marked4e = marked4e(:,[1 1 1]);
                    marked4ed = accumarray(obj.ed4e(:),marked4e(:));
                    obj.markedEdges = marked4ed>0;
                else
                    obj.markedEdges = numVector{1};
                end
            end
            if isempty(num) || ~uniform
                anz = 1;
            else
                anz = num{1};
            end
            if ~uniform && ~graded && nnz(obj.markedEdges)==0
                warning('fem4m:Mesh:noMarkedEdges','No marked edges. Use UNIFORM refinement.');
                uniform = true;
            end
            if graded
                % the numeric arguments after the type, by position
                first = find(cellfun(@(a)ischar(a) && strcmpi(a,type),varargin),1);
                args = varargin(first+1:end);
                args = [args(~cellfun(@ischar,args)),cell(1,4)];
            end
            if strcmpi(type,'GRADED')
                mu = args{1}; H = args{2}; gradPoints = args{3}; R = args{4};
                anz = 1;
            end
            if strcmpi(type,'GRADEDBETA')
                mu = args{1}; gradPoints = args{2}; R = args{3};
            end
            if strcmpi(type,'GRADEDRED')
                mu = args{1}; gradPoints = args{2};
                uniform = true;
            end

            if anz < 1
                % no round: the identity refinement
                obj.parent4e = (1:obj.nrElems)';
                obj.prolongation = speye(obj.nrNodes);
            end
            for k = 1 : anz
                if k > 1
                    parent4e = obj.parent4e;
                    P = obj.prolongation;
                end
                if uniform
                    obj.markedEdges = ones(size(obj.n4ed,1),1);
                end
                if strcmpi(type,'RGB')
                    obj = closure(obj); % closure
                    obj = refineRGB(obj); % RGB refinement
                elseif strcmpi(type,'NVB')
                    obj = closure(obj); % closure
                    obj = refineNVB(obj); % NVB refinement
                elseif strcmpi(type,'GRADED')
                    obj = refineGraded(obj,mu,H,gradPoints,R);
                elseif strcmpi(type,'GRADEDBETA')
                    obj = refineGradedBeta(obj,mu,gradPoints,R);
                elseif strcmpi(type,'GRADEDRED')
                    obj = refineGradedRed(obj,mu,gradPoints);
                end
                obj.markedEdges = zeros(size(obj.n4ed,1),1);
                obj.refineEdges = [];
                obj = finishRound(obj);
                if k > 1
                    % compose to the mesh before the refine call
                    obj.parent4e = parent4e(obj.parent4e);
                    obj.prolongation = obj.prolongation*P;
                end
            end
        end

        function obj = coarse(obj,elems)
            %COARSE - Coarsen marked element pairs
            %   OBJ = COARSE(OBJ,ELEMS) undoes the last newest vertex
            %   bisection at the newest nodes (third nodes) of the
            %   elements ELEMS, a logical vector with one entry per
            %   element or a vector of element ids: such a node is
            %   removed if it was created by a refinement and belongs to
            %   2 (boundary) or 4 (interior) elements, and the pairs of
            %   sibling elements at it are merged. parent4e and
            %   prolongation are set to [], and children4e records the
            %   merged pairs. If no element of ELEMS has a newest node
            %   created by a refinement, OBJ is returned unchanged.
            if isa(elems,'logical')
                obj = coarseNVB(obj,elems);
            else
                e = zeros(size(obj.n4e,1),1);
                e(elems) = 1;
                obj = coarseNVB(obj,e);
            end
        end

        function disp(obj)
            %DISP - Display summary or full data of the mesh
            %   DISP(OBJ) displays the node coordinates, the element node
            %   lists and the nonempty boundary lists for meshes with at
            %   most 100 nodes, and only the numbers of nodes, elements
            %   and edges otherwise.
            if obj.nrNodes <= 100
                disp('coordinates for nodes:')
                disp(obj.c4n)
                disp('nodes for elements:')
                disp(obj.n4e)
                if ~isempty(obj.n4Db)
                    disp('nodes for Dirichlet boundary:')
                    disp(obj.n4Db)
                end
                if ~isempty(obj.n4Nb)
                    disp('nodes for Neumann boundary:')
                    disp(obj.n4Nb)
                end
                if ~isempty(obj.n4Sb)
                    disp('nodes for symmetric boundary:')
                    disp(obj.n4Sb)
                end
            else
                disp(['no. nodes: ',num2str(obj.nrNodes)]);
                disp(['no. elems: ',num2str(obj.nrElems)]);
                disp(['no. edges: ',num2str(obj.nrEdges)]);
            end
        end

        function plot(obj,eind,clear)
            %PLOT - Plot the triangular mesh
            %   PLOT(OBJ) clears the current figure and plots the edges of
            %   all elements in black, or the elements as a white surface
            %   if the nodes have three coordinates (surface meshes).
            %
            %   PLOT(OBJ,EIND) plots only the elements selected by the
            %   index or logical vector EIND; an empty EIND selects all.
            %
            %   PLOT(OBJ,EIND,CLEAR) clears the figure only if CLEAR is
            %   true.
            if nargin == 1 || isempty(eind)
                eind = 1:obj.nrElems;
            end
            if nargin < 3 || clear
                clf;
            end
            if size(obj.c4n,2)==3
                trisurf(obj.n4e(eind,:),obj.c4n(:,1),obj.c4n(:,2),obj.c4n(:,3),zeros(obj.nrNodes,1),'FaceColor',[1 1 1])
            else
                triplot(obj.n4e(eind,:),obj.c4n(:,1),obj.c4n(:,2),'black')
            end
        end

        function plotBd(obj,color)
            %plotBd - Plot the boundary of the mesh
            %   plotBd(OBJ) draws the Dirichlet boundary as solid red
            %   lines, the Neumann boundary as dashed blue lines, and the
            %   symmetry boundary as dotted green lines into the current
            %   axes, which it leaves on hold.
            %
            %   plotBd(OBJ,COLOR) also uses the color COLOR, a char vector
            %   color name or short color spec such as 'k', for all three
            %   boundaries.
            hold on;
            if nargin < 2
                color = 'red';
            end
            if ~isempty(obj.n4Db)
                plot([obj.c4n(obj.n4Db(:,1),1)';obj.c4n(obj.n4Db(:,2),1)'],...
                    [obj.c4n(obj.n4Db(:,1),2)';obj.c4n(obj.n4Db(:,2),2)'],color,'LineWidth',2);
            end
            if nargin < 2
                color = 'blue';
            end
            if ~isempty(obj.n4Nb)
                plot([obj.c4n(obj.n4Nb(:,1),1)';obj.c4n(obj.n4Nb(:,2),1)'],...
                    [obj.c4n(obj.n4Nb(:,1),2)';obj.c4n(obj.n4Nb(:,2),2)'],['--',color],'LineWidth',2);
            end
            if nargin < 2
                color = 'green';
            end
            if ~isempty(obj.n4Sb)
                plot([obj.c4n(obj.n4Sb(:,1),1)';obj.c4n(obj.n4Sb(:,2),1)'],...
                    [obj.c4n(obj.n4Sb(:,1),2)';obj.c4n(obj.n4Sb(:,2),2)'],[':',color],'LineWidth',2);
            end
        end

        function plotEnumeration(obj)
            %plotEnumeration - Label nodes, elements, and edges with ids
            %   plotEnumeration(OBJ) writes the ids of the nodes (green),
            %   of the elements at their centroids (yellow) and of the
            %   edges at their midpoints (purple) as text into the current
            %   axes; it does not draw the mesh (see plot).
            mid4e = obj.point4e([1/3,1/3]);
            mid4ed = obj.point4ed(1/2);
            text(obj.c4n(:,1),obj.c4n(:,2),num2str((1:size(obj.c4n,1))'),...
                'color',[0.26 0.49 0.19],'FontWeight','bold','FontSize',24,'HorizontalAlignment','center');
            text(mid4e(:,1),mid4e(:,2),num2str((1:size(obj.n4e,1))'),...
                'color',[0.95 0.91 0],'FontWeight','bold','FontSize',24,'HorizontalAlignment','center');
            text(mid4ed(:,1),mid4ed(:,2),num2str((1:obj.nrEdges)'),...
                'color',[0.38 0.02 0.68],'FontWeight','bold','FontSize',24,'HorizontalAlignment','center');
        end

        function plotColor4n(obj,color)
            %plotColor4n - Plot mesh with node colors
            %   plotColor4n(OBJ) plots the mesh and marks every node with
            %   a color from the greedy coloring of color4n.
            %
            %   plotColor4n(OBJ,COLOR) also uses the given color vector
            %   COLOR instead of computing one.
            plot(obj);hold on;
            if nargin < 2 || isempty(color)
                color = obj.color4n;
            end
            scatter(obj.c4n(:,1),obj.c4n(:,2),200,color,'filled');
            colorbar;
        end

        function plotColor4e(obj,color)
            %plotColor4e - Plot mesh with element colors
            %   plotColor4e(OBJ) plots the mesh and marks every element
            %   center with a color from the greedy coloring of color4e.
            %
            %   plotColor4e(OBJ,COLOR) also uses the given color vector
            %   COLOR instead of computing one.
            plot(obj);hold on;
            if nargin < 2 || isempty(color)
                color = obj.color4e;
            end
            mid4e = obj.point4e([1/3,1/3]);
            scatter(mid4e(:,1),mid4e(:,2),200,color,'filled');
            colorbar;
        end

        function obj = removeNode(obj,nodes)
            %removeNode - Remove nodes and their elements
            %   OBJ = removeNode(OBJ,NODES) removes the nodes NODES (a
            %   row vector of node ids) together with all their elements.
            %   The exposed edges are appended to the Neumann boundary if
            %   the node lies on it, else to the symmetry boundary if it
            %   lies on that, else to the Dirichlet boundary. The result
            %   is a new Mesh built from the remaining lists, without
            %   refinement data and boundary projection.
            % node should not be connected to a node that does not share
            % another triangle, if so remove that node first
            nodes = sort(nodes,'descend');
            for curNode = nodes
                patch = nonzeros(obj.e4n(curNode,:));
                patchBd = zeros(numel(patch),2);
                for curElem = 1:numel(patch)
                    indx = find(obj.n4e(patch(curElem),:)==curNode);
                    patchBd(curElem,:) = [obj.n4e(patch(curElem),mod(indx+2-1,3)+1),...
                        obj.n4e(patch(curElem),mod(indx+1-1,3)+1)];
                end
                if ~isempty(find(obj.n4Nb==curNode, 1))
                    obj.n4Nb = [obj.n4Nb;patchBd];
                elseif ~isempty(obj.n4Sb) && ~isempty(find(obj.n4Sb(:,1:2)==curNode, 1))
                    % exposed edges inherit the symmetry axis of the node
                    sbAxis = obj.n4Sb(find(any(obj.n4Sb(:,1:2)==curNode,2),1),3);
                    obj.n4Sb = [obj.n4Sb;patchBd,sbAxis*ones(size(patchBd,1),1)];
                else
                    obj.n4Db = [obj.n4Db;patchBd];
                end
                obj.n4e(patch,:) = [];
                obj.n4e(obj.n4e(:)>curNode) = obj.n4e(obj.n4e(:)>curNode) - 1;
                if ~isempty(obj.n4Db)
                    obj.n4Db(obj.n4Db(:,1)==curNode | obj.n4Db(:,2)==curNode,:) = [];
                    obj.n4Db(obj.n4Db(:)>curNode) = obj.n4Db(obj.n4Db(:)>curNode) - 1;
                end
                if ~isempty(obj.n4Nb)
                    obj.n4Nb(obj.n4Nb(:,1)==curNode | obj.n4Nb(:,2)==curNode,:) = [];
                    obj.n4Nb(obj.n4Nb(:)>curNode) = obj.n4Nb(obj.n4Nb(:)>curNode) - 1;
                end
                if ~isempty(obj.n4Sb)
                    obj.n4Sb(obj.n4Sb(:,1)==curNode | obj.n4Sb(:,2)==curNode,:) = [];
                    sbNodes = obj.n4Sb(:,1:2);
                    sbNodes(sbNodes>curNode) = sbNodes(sbNodes>curNode) - 1;
                    obj.n4Sb(:,1:2) = sbNodes;
                end
                obj.c4n(curNode,:) = [];
                obj = Mesh(obj.c4n,obj.n4e,obj.n4Db,obj.n4Nb,obj.n4Sb);
            end
        end

        function obj = removeEdge(obj,curEdge)
            %removeEdge - Remove an edge and its elements
            %   OBJ = removeEdge(OBJ,CUREDGE) removes the edge CUREDGE (a
            %   single edge id; otherwise OBJ is returned unchanged)
            %   together with its adjacent elements. Their other edges
            %   are appended to the Neumann or symmetry boundary if
            %   CUREDGE belongs to it, else to the Dirichlet boundary. The
            %   result is a new Mesh built from the remaining lists.
            %remove bd edges first!!!
            % only one edge at a time!!!
            if isscalar(curEdge)
                patch = obj.e4ed(curEdge,:);
                for ind = 1 : 2
                    if patch(ind) ~=0
                        curElem = patch(ind);
                        if obj.edNr4ed(curEdge,ind) == 1
                            bd = [obj.n4e(curElem,[3 2]); obj.n4e(curElem,[1 3])];
                        elseif obj.edNr4ed(curEdge,ind) == 2
                            bd = [obj.n4e(curElem,[1 3]); obj.n4e(curElem,[2 1])];
                        else
                            bd = [obj.n4e(curElem,[2 1]); obj.n4e(curElem,[3 2])];
                        end
                        if ~isempty(find(obj.ed4Nb==curEdge, 1))
                            obj.n4Nb = [obj.n4Nb;bd];
                        elseif ~isempty(find(obj.ed4Sb==curEdge, 1))
                            % exposed sides inherit the symmetry axis of the edge
                            sbAxis = obj.n4Sb(obj.ed4Sb==curEdge,3);
                            obj.n4Sb = [obj.n4Sb;bd,sbAxis*ones(size(bd,1),1)];
                        else
                            obj.n4Db = [obj.n4Db;bd];
                        end
                    end
                end
                obj.n4e(patch(patch~=0),:) = [];
                obj.n4Db(obj.ed4Db==curEdge,:) = [];
                obj.n4Nb(obj.ed4Nb==curEdge,:) = [];
                obj.n4Sb(obj.ed4Sb==curEdge,:) = [];
                obj = Mesh(obj.c4n,obj.n4e,obj.n4Db,obj.n4Nb,obj.n4Sb);
            end
        end

        function obj = addTriangle(obj,nodes)
            %addTriangle - Add a triangle to the mesh
            %   OBJ = addTriangle(OBJ,NODES) adds the triangle with the
            %   existing node ids NODES (1-by-3) and updates the boundary
            %   lists: edges that become interior are removed from them,
            %   and its new edges are appended to the Dirichlet boundary.
            %   The result is a new Mesh built from the lists.
            obj.n4e(end+1,:) = nodes;

            indDb = [];
            indNb = [];
            indSb = [];

            for indEdge = 1 : 3
                curEdge = obj.ed4n(nodes(mod(indEdge-1,3)+1),nodes(mod(indEdge+1-1,3)+1));
                if curEdge == 0
                    obj.n4Db(end+1,:) = [nodes(mod(indEdge-1,3)+1),nodes(mod(indEdge+1-1,3)+1)];
                else
                    indDb = [indDb,find(obj.ed4Db==curEdge)]; %#ok<AGROW>
                    indNb = [indNb,find(obj.ed4Nb==curEdge)]; %#ok<AGROW>
                    indSb = [indSb,find(obj.ed4Sb==curEdge)]; %#ok<AGROW>
                end
            end

            obj.n4Db(indDb,:) = [];
            obj.n4Nb(indNb,:) = [];
            obj.n4Sb(indSb,:) = [];

            obj = Mesh(obj.c4n,obj.n4e,obj.n4Db,obj.n4Nb,obj.n4Sb);
        end

        function obj = save(obj,filename)
            %SAVE - Write mesh data to .dat files
            %   SAVE(OBJ) writes the node coordinates (17 significant
            %   digits), the element node lists and the boundary lists to
            %   the files mesh_c4n.dat, mesh_n4e.dat, mesh_n4Db.dat,
            %   mesh_n4Nb.dat and mesh_n4Sb.dat, which Mesh('mesh') reads
            %   back. An empty boundary list writes no file and deletes
            %   the file of an earlier save.
            %
            %   SAVE(OBJ,FILENAME) uses the char vector FILENAME as the
            %   prefix of the file names instead of 'mesh'.
            if nargin < 2 || isempty(filename)
                filename = 'mesh';
            end
            % all 17 digits, writematrix alone writes 15
            writematrix(compose("%.17g",obj.c4n),[filename,'_c4n.dat'], ...
                'Delimiter',' ','QuoteStrings','none');
            writematrix(obj.n4e,[filename,'_n4e.dat'],'Delimiter',' ');
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

        function obj = closure(obj)
            %CLOSURE - Complete marked edges for RGB refinement
            %   OBJ = CLOSURE(OBJ) sets OBJ.refineEdges to the marked
            %   edges OBJ.markedEdges and marks repeatedly the first edge
            %   of every triangle with a marked second or third edge,
            %   until the closure condition of red-green-blue and newest
            %   vertex bisection refinement holds.
            obj.refineEdges = obj.markedEdges;
            indx = find(obj.refineEdges(obj.ed4e(:,1))==0 & (obj.refineEdges(obj.ed4e(:,2))==1 | obj.refineEdges(obj.ed4e(:,3))==1));
            while ~isempty(indx)
                obj.refineEdges(obj.ed4e(indx,1)) = true;
                indx = find(obj.refineEdges(obj.ed4e(:,1))==0 & (obj.refineEdges(obj.ed4e(:,2))==1 | obj.refineEdges(obj.ed4e(:,3))==1));
            end
        end


    end

    methods (Access = protected)
        function obj = finishRound(obj)
            %finishRound - Complete a change of the mesh
            %   OBJ = finishRound(OBJ) projects the new boundary nodes onto
            %   the curved boundary of the geometry, if there is one,
            %   after every refinement round and coarsening. Subclasses
            %   override it to finish their own data.
            if ~isempty(obj.projectBoundary)
                obj = obj.projectBoundary(obj);
            end
        end
    end

    methods (Access = private)
        % Extended version of "EFFICIENT IMPLEMENTATION OF ADAPTIVE P1-FEM
        % IN MATLAB"
        function obj = getDataStructures(obj)
            %getDataStructures - Compute all derived mesh data structures
            obj.nrNodes = size(obj.c4n,1);
            obj.nrElems = size(obj.n4e,1);
            if isempty(obj.n4Db) && isempty(obj.n4Nb) && isempty(obj.n4Sb)
                twoN4ed = [obj.n4e(:,[1 2]);obj.n4e(:,[2 3]);obj.n4e(:,[3 1])];
            elseif isempty(obj.n4Nb) && isempty(obj.n4Sb)
                twoN4ed = [obj.n4e(:,[1 2]);obj.n4e(:,[2 3]);obj.n4e(:,[3 1]);obj.n4Db(:,[2 1])];
            elseif isempty(obj.n4Db) && isempty(obj.n4Sb)
                twoN4ed = [obj.n4e(:,[1 2]);obj.n4e(:,[2 3]);obj.n4e(:,[3 1]);obj.n4Nb(:,[2 1])];
            elseif isempty(obj.n4Db) && isempty(obj.n4Nb)
                twoN4ed = [obj.n4e(:,[1 2]);obj.n4e(:,[2 3]);obj.n4e(:,[3 1]);obj.n4Sb(:,[2 1])];
            elseif isempty(obj.n4Nb)
                twoN4ed = [obj.n4e(:,[1 2]);obj.n4e(:,[2 3]);obj.n4e(:,[3 1]);obj.n4Db(:,[2 1]);obj.n4Sb(:,[2 1])];
            elseif isempty(obj.n4Db)
                twoN4ed = [obj.n4e(:,[1 2]);obj.n4e(:,[2 3]);obj.n4e(:,[3 1]);obj.n4Nb(:,[2 1]);obj.n4Sb(:,[2 1])];
            elseif isempty(obj.n4Sb)
                twoN4ed = [obj.n4e(:,[1 2]);obj.n4e(:,[2 3]);obj.n4e(:,[3 1]);obj.n4Db(:,[2 1]);obj.n4Nb(:,[2 1])];
            else
                twoN4ed = [obj.n4e(:,[1 2]);obj.n4e(:,[2 3]);obj.n4e(:,[3 1]);obj.n4Db(:,[2 1]);obj.n4Nb(:,[2 1]);obj.n4Sb(:,[2 1])];
            end
            indx = find(twoN4ed(:,1) < twoN4ed(:,2));
            indy = find(twoN4ed(:,1) > twoN4ed(:,2));
            obj.n4ed = twoN4ed(indx,:);
            obj.nrEdges = size(obj.n4ed,1);
            twoEd = zeros(size(twoN4ed,1),1);
            twoEd(indx) = (1:size(obj.n4ed,1));
            obj.ed4n = sparse(obj.n4ed(:,1),obj.n4ed(:,2),twoEd(indx)',size(obj.c4n,1),size(obj.c4n,1));
            [~,~,permx] = find(obj.ed4n);
            [~,~,permy] = find(sparse(twoN4ed(indy,2),twoN4ed(indy,1),indy));
            twoEd(permy) = permx;
            obj.ed4n = obj.ed4n + obj.ed4n';
            obj.ed4e = reshape(twoEd(1:numel(obj.n4e)),size(obj.n4e));
            obj.ed4Db = twoEd(numel(obj.n4e)+1:numel(obj.n4e)+size(obj.n4Db,1));
            obj.ed4Nb = twoEd(numel(obj.n4e)+size(obj.n4Db,1)+1:numel(obj.n4e)+size(obj.n4Db,1)+size(obj.n4Nb,1));
            obj.ed4Sb = twoEd(numel(obj.n4e)+size(obj.n4Db,1)+size(obj.n4Nb,1)+1:end);
            obj.e4ed = zeros(size(obj.n4ed,1),2);
            if ~isempty(obj.n4Db)
                obj.n4ed(obj.ed4Db,:) = obj.n4Db;
            end
            if ~isempty(obj.n4Nb)
                obj.n4ed(obj.ed4Nb,:) = obj.n4Nb;
            end
            if ~isempty(obj.n4Sb)
                obj.n4ed(obj.ed4Sb,:) = obj.n4Sb(:,[1 2]);
            end
            elems = (1:size(obj.n4e,1))'*ones(1,size(obj.n4e,2));
            elems = elems(:);
            indx = indx(indx <= numel(obj.n4e));
            indy = indy(indy <= numel(obj.n4e));
            obj.e4ed(twoEd(indx),1) = elems(indx);
            obj.e4ed(twoEd(indy),2) = elems(indy);
            indz = find(obj.e4ed(:,1)==0);
            obj.e4ed(indz,1) = obj.e4ed(indz,2);
            obj.e4ed(indz,2) = 0;
            obj.e4n = sparse(obj.n4e(:),elems,elems,obj.nrNodes,obj.nrElems);
            obj.edNr4ed = zeros(size(obj.e4ed));
            ind = (1:size(obj.n4ed,1))';
            indz = find(obj.e4ed(:,2)>0);
            obj.edNr4ed(:,1) = mod(find(obj.ed4e(obj.e4ed(:,1),:)'==(ind*ones(1,3))')-1,3)+1;
            obj.edNr4ed(indz,2) = mod(find(obj.ed4e(obj.e4ed(indz,2),:)'==(ind(indz)*ones(1,3))')-1,3)+1;
            obj.DbEd4n = zeros(obj.nrNodes,2);
            obj.NbEd4n = zeros(obj.nrNodes,2);
            obj.SbEd4n = zeros(obj.nrNodes,2);
            if ~isempty(obj.n4Db)
                obj.DbEd4n(obj.n4Db(:,1),1) = obj.ed4Db;
                obj.DbEd4n(obj.n4Db(:,2),2) = obj.ed4Db;
            end
            if ~isempty(obj.n4Nb)
                obj.NbEd4n(obj.n4Nb(:,1),1) = obj.ed4Nb;
                obj.NbEd4n(obj.n4Nb(:,2),2) = obj.ed4Nb;
            end
            if ~isempty(obj.n4Sb)
                obj.SbEd4n(obj.n4Sb(:,1),1) = obj.ed4Sb;
                obj.SbEd4n(obj.n4Sb(:,2),2) = obj.ed4Sb;
            end
            nNr = ((1:3)'*ones(1,size(obj.n4e,1)))';
            obj.nNr4n = sparse(twoN4ed(1:numel(obj.n4e),1),elems,nNr(:));
        end

        % red green blue
        function obj = refineRGB(obj)
            %refineRGB - Refine along marked edges by red-green-blue
            mid4ed = obj.point4ed(1/2);
            c4nNew = [obj.c4n;mid4ed(obj.refineEdges==1,:)];
            newNode4ed = zeros(size(obj.n4ed,1),1);
            newNode4ed(obj.refineEdges==1) = size(obj.c4n,1)+ (1 : nnz(obj.refineEdges));
            newNode4ed4e = reshape(newNode4ed(obj.ed4e),[],3);
            %no refinement
            indx = find(~newNode4ed4e(:,1) & ~newNode4ed4e(:,2) & ~newNode4ed4e(:,3));
            n4eNew = obj.n4e(indx,:);
            parents = indx;
            %red
            indx = find(newNode4ed4e(:,1) & newNode4ed4e(:,2) & newNode4ed4e(:,3));
            n4eNew = [n4eNew;...
                obj.n4e(indx,1),newNode4ed4e(indx,1),newNode4ed4e(indx,3);...
                newNode4ed4e(indx,1),obj.n4e(indx,2), newNode4ed4e(indx,2);...
                newNode4ed4e(indx,3),newNode4ed4e(indx,2),obj.n4e(indx,3);...
                newNode4ed4e(indx,2),newNode4ed4e(indx,3),newNode4ed4e(indx,1)];
            parents = [parents;indx;indx;indx;indx];
            %green
            indx = find(newNode4ed4e(:,1) & ~newNode4ed4e(:,2) & ~newNode4ed4e(:,3));
            n4eNew = [n4eNew;...
                obj.n4e(indx,3),obj.n4e(indx,1),newNode4ed4e(indx,1);...
                obj.n4e(indx,2),obj.n4e(indx,3),newNode4ed4e(indx,1)];
            parents = [parents;indx;indx];
            %blue left
            indx = find(newNode4ed4e(:,1) & ~newNode4ed4e(:,2) & newNode4ed4e(:,3));
            n4eNew = [n4eNew;...
                newNode4ed4e(indx,1),obj.n4e(indx,3),newNode4ed4e(indx,3);...
                obj.n4e(indx,1),newNode4ed4e(indx,1),newNode4ed4e(indx,3);...
                obj.n4e(indx,2),obj.n4e(indx,3),newNode4ed4e(indx,1)];
            parents = [parents;indx;indx;indx];
            %blue right
            indx = find(newNode4ed4e(:,1) & newNode4ed4e(:,2) & ~newNode4ed4e(:,3));
            n4eNew = [n4eNew;...
                obj.n4e(indx,3),obj.n4e(indx,1),newNode4ed4e(indx,1);...
                newNode4ed4e(indx,1),obj.n4e(indx,2), newNode4ed4e(indx,2);...
                obj.n4e(indx,3),newNode4ed4e(indx,1),newNode4ed4e(indx,2)];
            parents = [parents;indx;indx;indx];
            % Boundary
            if ~isempty(obj.n4Db)
                newNode4Db = newNode4ed(obj.ed4Db);
                indx = find(newNode4Db~=0);
                n4DbNew = [obj.n4Db(~newNode4Db,:);...
                    obj.n4Db(indx,1),newNode4Db(indx);...
                    newNode4Db(indx),obj.n4Db(indx,2)];
            else
                n4DbNew = [];
            end
            if ~isempty(obj.n4Nb)
                newNode4Nb = newNode4ed(obj.ed4Nb);
                indx = find(newNode4Nb~=0);
                n4NbNew = [obj.n4Nb(~newNode4Nb,:);...
                    obj.n4Nb(indx,1),newNode4Nb(indx);...
                    newNode4Nb(indx),obj.n4Nb(indx,2)];
            else
                n4NbNew = [];
            end
            if ~isempty(obj.n4Sb)
                newNode4Sb = newNode4ed(obj.ed4Sb);
                indx = find(newNode4Sb~=0);
                n4SbNew = [obj.n4Sb(~newNode4Sb,:);...
                    obj.n4Sb(indx,1),newNode4Sb(indx),obj.n4Sb(indx,3);...
                    newNode4Sb(indx),obj.n4Sb(indx,2),obj.n4Sb(indx,3)];
            else
                n4SbNew = [];
            end
            I = ([1:size(obj.c4n,1),newNode4ed(obj.refineEdges==1)']'*ones(1,2))';
            J = [(1:size(obj.c4n,1))'*ones(1,2);obj.n4ed(obj.refineEdges==1,:)]';
            S = 1/2 * ones(numel(I),1);
            P = sparse(I(:),J(:),S,size(c4nNew,1),size(obj.c4n,1));
            obj.prolongation = P;
            obj.parent4e = parents;
            obj.c4n = c4nNew;
            obj.n4e = n4eNew;
            obj.n4Db = n4DbNew;
            obj.n4Nb = n4NbNew;
            obj.n4Sb = n4SbNew;
            obj = getDataStructures(obj);
        end

        %  NVB refinement
        function obj = refineNVB(obj)
            %refineNVB - Refine along marked edges by newest vertex bisection
            mid4ed = obj.point4ed(1/2);
            c4nNew = [obj.c4n;mid4ed(obj.refineEdges==1,:)];
            newNode4ed = zeros(size(obj.n4ed,1),1);
            newNode4ed(obj.refineEdges==1) = size(obj.c4n,1)+ (1 : nnz(obj.refineEdges));
            newNode4ed4e = newNode4ed(obj.ed4e);
            if obj.nrElems == 1
                newNode4ed4e = newNode4ed4e';
            end
            %type
            %indxNo = find(~newNode4ed4e(:,1) & ~newNode4ed4e(:,2) & ...
            %    ~newNode4ed4e(:,3));
            indxGreen =  find(newNode4ed4e(:,1) & ~newNode4ed4e(:,2) & ~newNode4ed4e(:,3));
            indxBlueLeft =  find(newNode4ed4e(:,1) & ~newNode4ed4e(:,2) & newNode4ed4e(:,3));
            indxBlueRight =  find(newNode4ed4e(:,1) & newNode4ed4e(:,2) & ~newNode4ed4e(:,3));
            indxbisec3 =  find(newNode4ed4e(:,1) & newNode4ed4e(:,2) & newNode4ed4e(:,3));
            % the first (left) new element inherits the element number
            n4eNew = obj.n4e;
            n4eNew(indxGreen,:) = [obj.n4e(indxGreen,3),obj.n4e(indxGreen,1),newNode4ed4e(indxGreen,1)];
            n4eNew(indxBlueLeft,:) = [newNode4ed4e(indxBlueLeft,1),obj.n4e(indxBlueLeft,3),newNode4ed4e(indxBlueLeft,3)];
            n4eNew(indxBlueRight,:) = [obj.n4e(indxBlueRight,3),obj.n4e(indxBlueRight,1),newNode4ed4e(indxBlueRight,1)];
            n4eNew(indxbisec3,:) = [newNode4ed4e(indxbisec3,1),obj.n4e(indxbisec3,3),newNode4ed4e(indxbisec3,3)];
            parents = (1:obj.nrElems)';
            % bisec3
            indx = find(newNode4ed4e(:,1) & newNode4ed4e(:,2) & newNode4ed4e(:,3));
            n4eNew = [n4eNew;...
                obj.n4e(indx,1),newNode4ed4e(indx,1),newNode4ed4e(indx,3);...
                newNode4ed4e(indx,1),obj.n4e(indx,2), newNode4ed4e(indx,2);...
                obj.n4e(indx,3),newNode4ed4e(indx,1),newNode4ed4e(indx,2)];
            parents = [parents;indx;indx;indx];
            %green
            indx = find(newNode4ed4e(:,1) & ~newNode4ed4e(:,2) & ~newNode4ed4e(:,3));
            n4eNew = [n4eNew;...
                obj.n4e(indx,2),obj.n4e(indx,3),newNode4ed4e(indx,1)];
            parents = [parents;indx];
            %blue left
            indx = find(newNode4ed4e(:,1) & ~newNode4ed4e(:,2) & newNode4ed4e(:,3));
            n4eNew = [n4eNew;...
                obj.n4e(indx,1),newNode4ed4e(indx,1),newNode4ed4e(indx,3);...
                obj.n4e(indx,2),obj.n4e(indx,3),newNode4ed4e(indx,1)];
            parents = [parents;indx;indx];
            %blue right
            indx = find(newNode4ed4e(:,1) & newNode4ed4e(:,2) & ~newNode4ed4e(:,3));
            n4eNew = [n4eNew;...
                newNode4ed4e(indx,1),obj.n4e(indx,2), newNode4ed4e(indx,2);...
                obj.n4e(indx,3),newNode4ed4e(indx,1),newNode4ed4e(indx,2)];
            parents = [parents;indx;indx];
            % Boundary
            if ~isempty(obj.n4Db)
                newNode4Db = newNode4ed(obj.ed4Db);
                indx = find(newNode4Db~=0);
                n4DbNew = [obj.n4Db(~newNode4Db,:);...
                    obj.n4Db(indx,1),newNode4Db(indx);...
                    newNode4Db(indx),obj.n4Db(indx,2)];
            else
                n4DbNew = [];
            end
            if ~isempty(obj.n4Nb)
                newNode4Nb = newNode4ed(obj.ed4Nb);
                indx = find(newNode4Nb~=0);
                n4NbNew = [obj.n4Nb(~newNode4Nb,:);...
                    obj.n4Nb(indx,1),newNode4Nb(indx);...
                    newNode4Nb(indx),obj.n4Nb(indx,2)];
            else
                n4NbNew = [];
            end
            if ~isempty(obj.n4Sb)
                newNode4Sb = newNode4ed(obj.ed4Sb);
                indx = find(newNode4Sb~=0);
                n4SbNew = [obj.n4Sb(~newNode4Sb,:);...
                    obj.n4Sb(indx,1),newNode4Sb(indx),obj.n4Sb(indx,3);...
                    newNode4Sb(indx),obj.n4Sb(indx,2),obj.n4Sb(indx,3)];
            else
                n4SbNew = [];
            end
            I = ([1:size(obj.c4n,1),newNode4ed(obj.refineEdges==1)']'*ones(1,2))';
            J = [(1:size(obj.c4n,1))'*ones(1,2);obj.n4ed(obj.refineEdges==1,:)]';
            S = 1/2 * ones(numel(I),1);
            P = sparse(I(:),J(:),S,size(c4nNew,1),size(obj.c4n,1));
            obj.parent4e = parents;
            obj.prolongation = P;
            obj.c4n = c4nNew;
            obj.n4e = n4eNew;
            obj.n4Db = n4DbNew;
            obj.n4Nb = n4NbNew;
            obj.n4Sb = n4SbNew;
            obj = getDataStructures(obj);
        end

        function obj = coarseNVB(obj,markedElems)
            %coarseNVB - Undo newest vertex bisection for marked elements
            % Find marked Nodes
            newestNode = unique(obj.n4e(markedElems & obj.n4e(:,3)>obj.N0,3));
            if isempty(newestNode); return; end
            valence = accumarray(obj.n4e(:),1,[size(obj.c4n,1) 1]);
            markedNodes = zeros(size(obj.c4n,1),1);
            markedNodes(newestNode((valence(newestNode) == 2 | valence(newestNode) == 4))) = 1;
            % Find elements
            brother = obj.e4ed(obj.ed4e(:,3),:);
            brother(brother(:,1)==(1:size(obj.n4e,1))',1) = brother(brother(:,1)==(1:size(obj.n4e,1))',2);
            brother = brother(:,1);
            elems = find(markedNodes(obj.n4e(:,3))  & (brother > (1:size(obj.n4e,1))') );
            markedElements = zeros(size(obj.n4e,1),1);
            markedElements(elems) = 1;
            for element = elems'
                if markedElements(element)
                    markedElements(brother(element)) = 0;
                end
            end
            elems = find(markedElements);
            % Coarse
            n4eNew = obj.n4e;
            n4eNew(elems,[1 3 2]) = [n4eNew(elems,[2 1]) n4eNew(brother(elems),1)];
            n4eNew(brother(elems),:) = [];
            obj.children4e = (1:obj.nrElems)';
            obj.children4e(elems,2) = brother(elems);
            obj.children4e(brother(elems),:) = 0;
            % Update nodes
            activeNodes = find(markedNodes==0);
            c4nNew = obj.c4n(activeNodes,:);
            n4nNew = zeros(size(obj.c4n,1),1);
            n4nNew(activeNodes) = 1:length(activeNodes);
            n4eNew = reshape(n4nNew(n4eNew),[],3);
            %Db
            if ~isempty(obj.n4Db)
                n4b = zeros(size(obj.c4n,1),2);
                n4DbNew = obj.n4Db;
                n4b(n4DbNew(:,1),1) = 1:size(n4DbNew,1);
                n4b(n4DbNew(:,2),2) = 1:size(n4DbNew,1);
                idx = ( markedNodes & n4b(:,2) );
                n4DbNew(n4b(idx,2),2) = n4DbNew(n4b(idx,1),2);
                n4DbNew(n4b(idx,1),2) = 0;
                n4DbNew = reshape(n4nNew(n4DbNew(n4DbNew(:,2)>0,:)),[],2);
            else
                n4DbNew = [];
            end
            %Nb
            if ~isempty(obj.n4Nb)
                n4b = zeros(size(obj.c4n,1),2);
                n4NbNew = obj.n4Nb;
                n4b(n4NbNew(:,1),1) = 1:size(n4NbNew,1);
                n4b(n4NbNew(:,2),2) = 1:size(n4NbNew,1);
                idx = ( markedNodes & n4b(:,2) );
                n4NbNew(n4b(idx,2),2) = n4NbNew(n4b(idx,1),2);
                n4NbNew(n4b(idx,1),2) = 0;
                n4NbNew = reshape(n4nNew(n4NbNew(n4NbNew(:,2)>0,:)),[],2);
            else
                n4NbNew = [];
            end
            %Sb
            if ~isempty(obj.n4Sb)
                n4b = zeros(size(obj.c4n,1),2);
                n4SbNew = obj.n4Sb;
                n4b(n4SbNew(:,1),1) = 1:size(n4SbNew,1);
                n4b(n4SbNew(:,2),2) = 1:size(n4SbNew,1);
                idx = ( markedNodes & n4b(:,2) );
                n4SbNew(n4b(idx,2),2) = n4SbNew(n4b(idx,1),2);
                n4SbNew(n4b(idx,1),2) = 0;
                n4SbNew = n4SbNew(n4SbNew(:,2)>0,:);
                % renumber the nodes, not the axis in the third column
                n4SbNew(:,1:2) = reshape(n4nNew(n4SbNew(:,1:2)),[],2);
            else
                n4SbNew = [];
            end
            % coarsening is no refinement: no parent elements, no prolongation
            obj.parent4e = [];
            obj.prolongation = [];
            obj.c4n = c4nNew;
            obj.n4e = n4eNew;
            obj.n4Db = n4DbNew;
            obj.n4Nb = n4NbNew;
            obj.n4Sb = n4SbNew;
            obj = getDataStructures(obj);
            if ~isempty(elems)
                obj = finishRound(obj);
            end
        end

        function obj = refineGraded(obj,mu,H,gradPoints,R)
            %refineGraded - Refine with a mesh size graded around points
            if nargin < 4 || isempty(gradPoints)
                gradPoints = [0,0];
            end
            if nargin < 5 || isempty(R)
                R = 1;
            end
            % one exponent and radius per grading point
            mu = mu(:).*ones(size(gradPoints,1),1);
            R = R(:).*ones(size(gradPoints,1),1);
            parent4e = (1:obj.nrElems)';
            P = speye(obj.nrNodes);
            count = 1;
            while count ~= 0
                midPoint4e = obj.point4e([1/3,1/3]);
                obj.markedEdges = false(obj.nrEdges,1);
                refineElems = false(obj.nrElems,1);
                length4ed = obj.length4ed;
                hT = max(length4ed(obj.ed4e),[],2);
                refineElems(hT > H) = true;
                for k = 1 : size(gradPoints,1)
                    rT = sqrt( (midPoint4e(:,1)-gradPoints(k,1)).^2+...
                        (midPoint4e(:,2)-gradPoints(k,2)).^2  );
                    refineElems(hT > H.*(rT/R(k)).^(1-mu(k))) = true;
                end
                count = nnz(refineElems);
                obj.markedEdges(obj.ed4e(refineElems,1)) = true;
                obj = closure(obj);
                obj = refineRGB(obj);
                obj = finishRound(obj);
                parent4e = parent4e(obj.parent4e);
                P = obj.prolongation*P;
            end
            obj.parent4e = parent4e;
            obj.prolongation = P;
        end

        function obj = refineGradedBeta(obj,mu,gradPoints,R)
            %refineGradedBeta - Build a beta-graded mesh by shifting nodes
            if nargin < 3 || isempty(gradPoints)
                gradPoints = [0,0];
            end
            if nargin < 4 || isempty(R)
                R = 1;
            end
            % one exponent and radius per grading point
            mu = mu(:).*ones(size(gradPoints,1),1);
            R = R(:).*ones(size(gradPoints,1),1);

            % refine uniform mesh
            if ~isempty(obj.uniformC4n)
                obj.c4n = obj.uniformC4n;
            end
            obj = obj.refine('uniform','rgb',1);
            obj.uniformC4n = obj.c4n;

            % Shift coordinates
            for k = 1 : size(gradPoints,1)
                rT = sqrt( (obj.c4n(:,1)-gradPoints(k,1)).^2+...
                    (obj.c4n(:,2)-gradPoints(k,2)).^2  );
                d(:,1) = -(gradPoints(k,1) - obj.c4n(:,1))./rT*R(k);
                d(:,2) = -(gradPoints(k,2) - obj.c4n(:,2))./rT*R(k);
                %move mesh points
                indx = rT < R(k) & rT>0;
                obj.c4n(indx,:) = gradPoints(k*ones(nnz(indx),1),:)...
                    +(rT(indx)/R(k)).^(1/mu(k)).*d(indx,:);
            end
        end

        % graded red
        function obj = refineGradedRed(obj,kappa,gradPoints)
            %refineGradedRed - Red-refine with graded midpoints
            %kappa = 2.^(-1./mu);
            % strictly less than optimal value, i.e. round down to one
            % decimal
            %kappa = max(floor(kappa*10)/10,0.1);
            if nargin < 3 || isempty(gradPoints) || isempty(kappa)
                error('fem4m:Mesh:gradedRedArgs', ...
                    '''gradedred'' needs the factor kappa and at least one grading point.');
            end
            % one factor and node per grading point
            kappa = kappa(:).*ones(size(gradPoints,1),1);
            node4p = zeros(size(gradPoints,1),1);
            for k = 1 : size(gradPoints,1)
                curNode = find(obj.c4n(:,1)==gradPoints(k,1) & obj.c4n(:,2)==gradPoints(k,2),1);
                if isempty(curNode)
                    error('fem4m:Mesh:gradPointNotNode', ...
                        'The grading point (%g,%g) of ''gradedred'' must be a node of the mesh.',gradPoints(k,1),gradPoints(k,2));
                end
                node4p(k) = curNode;
            end
            c4nNew = obj.point4ed(1/2);
            %grade singular nodes
            for k = 1 : size(gradPoints,1)
                curNode = node4p(k);
                gradedA4ed = obj.c4n(obj.n4ed(:,1),:) + kappa(k) * ( obj.c4n(obj.n4ed(:,2),:) - obj.c4n(obj.n4ed(:,1),:) );
                gradedB4ed = obj.c4n(obj.n4ed(:,2),:) + kappa(k) * ( obj.c4n(obj.n4ed(:,1),:) - obj.c4n(obj.n4ed(:,2),:) );
                c4nNew(obj.n4ed(:,1)==curNode,:) = gradedA4ed(obj.n4ed(:,1)==curNode,:);
                c4nNew(obj.n4ed(:,2)==curNode,:) = gradedB4ed(obj.n4ed(:,2)==curNode,:);
            end
            c4nNew = [obj.c4n;c4nNew];
            %refine red
            newNode4ed = size(obj.c4n,1) + (1 : obj.nrEdges)';
            newNode4ed4e = reshape(newNode4ed(obj.ed4e),[],3);
            n4eNew = [obj.n4e(:,1),newNode4ed4e(:,1),newNode4ed4e(:,3);...
                newNode4ed4e(:,1),obj.n4e(:,2), newNode4ed4e(:,2);...
                newNode4ed4e(:,3),newNode4ed4e(:,2),obj.n4e(:,3);...
                newNode4ed4e(:,2),newNode4ed4e(:,3),newNode4ed4e(:,1)];
            parents = [(1:obj.nrElems)';(1:obj.nrElems)';(1:obj.nrElems)';(1:obj.nrElems)'];
            % Boundary
            if ~isempty(obj.n4Db)
                newNode4Db = newNode4ed(obj.ed4Db);
                n4DbNew = [obj.n4Db(~newNode4Db,:);...
                    obj.n4Db(:,1),newNode4Db(:);...
                    newNode4Db(:),obj.n4Db(:,2)];
            else
                n4DbNew = [];
            end
            if ~isempty(obj.n4Nb)
                newNode4Nb = newNode4ed(obj.ed4Nb);
                n4NbNew = [obj.n4Nb(~newNode4Nb,:);...
                    obj.n4Nb(:,1),newNode4Nb(:);...
                    newNode4Nb(:),obj.n4Nb(:,2)];
            else
                n4NbNew = [];
            end
            if ~isempty(obj.n4Sb)
                newNode4Sb = newNode4ed(obj.ed4Sb);
                n4SbNew = [obj.n4Sb(~newNode4Sb,:);...
                    obj.n4Sb(:,1),newNode4Sb(:),obj.n4Sb(:,3);...
                    newNode4Sb(:),obj.n4Sb(:,2),obj.n4Sb(:,3)];
            else
                n4SbNew = [];
            end
            % P1 weights of the new nodes, moved towards the grading points
            weight = 1/2*ones(obj.nrEdges,2);
            for k = 1 : size(gradPoints,1)
                curNode = node4p(k);
                weight(obj.n4ed(:,1)==curNode,:) = ones(nnz(obj.n4ed(:,1)==curNode),1)*[1-kappa(k),kappa(k)];
                weight(obj.n4ed(:,2)==curNode,:) = ones(nnz(obj.n4ed(:,2)==curNode),1)*[kappa(k),1-kappa(k)];
            end
            P = [speye(size(obj.c4n,1));
                sparse([1:obj.nrEdges,1:obj.nrEdges]',double(obj.n4ed(:)),weight(:), ...
                obj.nrEdges,size(obj.c4n,1))];
            obj.prolongation = P;
            obj.parent4e = parents;
            obj.c4n = c4nNew;
            obj.n4e = n4eNew;
            obj.n4Db = n4DbNew;
            obj.n4Nb = n4NbNew;
            obj.n4Sb = n4SbNew;
            obj = getDataStructures(obj);
        end

    end
end