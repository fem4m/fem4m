classdef MeshPoly < handle
    %MeshPoly - Polygon mesh with arbitrary polygonal elements
    %   OBJ = MeshPoly(GEOM) creates the polygon mesh of the Geometry
    %   object GEOM with its boundary node pairs; an options struct GEOM
    %   is passed to Geometry first.
    %
    %   OBJ = MeshPoly(NAME) loads the files NAME_c4n.dat and
    %   NAME_n4e.dat (one line of node ids per element) and the existing
    %   boundary files NAME_n4Db.dat, NAME_n4Nb.dat and NAME_n4Sb.dat,
    %   e.g. written by save.
    %
    %   OBJ = MeshPoly(c4n,n4e) creates the mesh with the node
    %   coordinates c4n and the node ids n4e of the elements (a cell array
    %   or a matrix); all boundary edges are Dirichlet edges.
    %
    %   OBJ = MeshPoly(...,n4Db) also specifies the node pairs of the
    %   Dirichlet boundary.
    %
    %   OBJ = MeshPoly(...,n4Nb) also specifies the node pairs of the
    %   Neumann boundary.
    %
    %   OBJ = MeshPoly(...,n4Sb) also specifies the node pairs of the
    %   symmetry boundary, whose third column is the fixed axis (1 or 2,
    %   0 for a general normal).
    %
    %   In every form, the node pairs of n4Db, n4Nb and n4Sb (given, read
    %   from the files or taken from GEOM) must together cover the
    %   boundary, each pair in the node order of its element; otherwise
    %   the edge numbering is incomplete or wrong.
    %
    %   Row y of every x4y lists the ids of the x of y (n4e and ed4e are
    %   cell arrays): e4n is the sparse nrNodes-by-nrElems incidence with
    %   e4n(n,e) = e for the nodes n of the element e, e4ed(ed,:) are the
    %   two elements of the edge ed (0 on the boundary), and ed4n is the
    %   symmetric node-pair lookup ed4n(i,j) of the edge between the
    %   nodes i and j.
    %
    %   MeshPoly itself cannot be refined: refine raises an error, and
    %   the subclasses MeshPolyT, MeshPolyQ and MeshPolyH implement it.
    %   Unlike Mesh, the polygon meshes store no refinement data such as
    %   parent4e or a prolongation matrix.
    %
    %   The property convex is a user setting (default true), not a
    %   test of the elements: integrate splits every element into a fan
    %   of triangles from its point mid4e (true) or from its first node
    %   (false), each weighted with its signed area; the quadrature
    %   points lie inside the element only if it is star-shaped with
    %   respect to that point.
    %
    %   MeshPoly functions:
    %       length4ed         - Compute the lengths of all edges
    %       point4ed          - Compute points on all edges
    %       normal4ed         - Compute the unit normals of all edges
    %       tangent4ed        - Compute the unit tangents of all edges
    %       tangent4e         - Compute the unit tangents of element edges
    %       normal4e          - Compute the unit normals of element edges
    %       edNr4ed           - Compute local edge numbers in elements
    %       computeColor4n    - Color nodes of neighboring elements
    %       area4e            - Compute the areas of all elements
    %       diam4e            - Compute the diameters of all elements
    %       mid4e             - Compute the centroids of all elements
    %       disp              - Display summary or full data of the mesh
    %       plot              - Plot the polygonal mesh
    %       plotBd            - Plot the boundary of the mesh
    %       plotEnumeration   - Label nodes, elements, and edges with ids
    %       plotColor4n       - Plot the mesh colored by node colors
    %       plus              - Union of two polygon meshes
    %       refine            - Refine marked elements
    %       refineGraded      - Refine the mesh with a graded distribution
    %       subtriangulateMid - Split polygons into triangles around mid4e
    %       save              - Write the mesh data to .dat files
    %       removeGhostNodes  - Remove unused nodes from mesh data
    %       signedArea        - Twice the signed areas of triangles
    %
    %   MeshPoly properties:
    %       c4n     - Coordinates of nodes
    %       n4e     - Cell array of node ids for elements
    %       n4Db    - Node pairs of Dirichlet boundary
    %       n4Nb    - Node pairs of Neumann boundary
    %       n4Sb    - Node pairs of symmetry boundary plus fixed axis
    %       n4ed    - Node pairs of unique edges
    %       ed4n    - Sparse matrix of edge ids for node pairs
    %       ed4e    - Cell array of edge ids for elements
    %       e4ed    - Adjacent element ids for each edge
    %       e4n     - Sparse matrix of element ids for nodes
    %       ed4Db   - Edge ids of Dirichlet boundary
    %       ed4Nb   - Edge ids of Neumann boundary
    %       ed4Sb   - Edge ids of symmetry boundary
    %       NbEd4n  - Neumann edge ids incident to each node
    %       SbEd4n  - Symmetry edge ids incident to each node
    %       DbEd4n  - Dirichlet edge ids incident to each node
    %       nrNodes - Number of nodes
    %       nrEdges - Number of edges
    %       nrElems - Number of elements
    %       nrN4e   - Number of nodes per element
    %       convex  - Integrate on fans from mid4e (true) or node 1
    %
    %   See also Mesh, Mesh3d, MeshPolyH, MeshPolyT, MeshPolyQ

    % Copyright (c) 2026 by the fem4m authors

    properties
        c4n     % coordinates of nodes
        n4e     % cell array of node ids for elements
        n4Db    % node pairs of Dirichlet boundary
        n4Nb    % node pairs of Neumann boundary
        n4Sb    % node pairs of symmetry boundary plus fixed axis
        n4ed    % node pairs of unique edges
        ed4n    % sparse matrix of edge ids for node pairs
        ed4e    % cell array of edge ids for elements
        e4ed    % adjacent element ids for each edge
        e4n     % sparse matrix of element ids for nodes
        ed4Db   % edge ids of Dirichlet boundary
        ed4Nb   % edge ids of Neumann boundary
        ed4Sb   % edge ids of symmetry boundary
        NbEd4n  % Neumann edge ids incident to each node
        SbEd4n  % symmetry edge ids incident to each node
        DbEd4n  % Dirichlet edge ids incident to each node
        nrNodes % number of nodes
        nrEdges % number of edges
        nrElems % number of elements
        nrN4e   % number of nodes per element
        convex = true % integrate on fans from mid4e (true) or node 1
    end

    properties(Access=protected)
        area4eSave
        diam4eSave
        mid4eSave
        tangent4eSave
        normal4eSave
        edNr4edSave
    end

    methods
        function obj = MeshPoly(c4n,n4e,n4Db,n4Nb,n4Sb)
            %MeshPoly - Create a polygon mesh object
            if nargin < 5
                n4Sb = [];
            end
            explicit = isnumeric(c4n);
            if ischar(c4n)
                name = c4n;
                c4n   = load([name,'_c4n.dat']);
                n4e   = readmatrix([name,'_n4e.dat'],'NumHeaderLines',0);
                n4e = mat2cell(n4e,ones(size(n4e,1),1));
                n4e = cellfun(@(z)(z(~isnan(z))),n4e,'UniformOutput',false);
                if exist([name,'_n4Db.dat'],'file')
                    n4Db = load([name,'_n4Db.dat']);
                else
                    n4Db = [];
                end
                if exist([name,'_n4Nb.dat'],'file')
                    n4Nb  = load([name,'_n4Nb.dat']);
                else
                    n4Nb = [];
                end
                if exist([name,'_n4Sb.dat'],'file')
                    n4Sb  = load([name,'_n4Sb.dat']);
                else
                    n4Sb = [];
                end
            elseif isstruct(c4n) || strcmpi(class(c4n),'Geometry')
                geom = c4n;
                if isstruct(c4n)
                    geom = Geometry(geom);
                end
                c4n = geom.c4n;
                n4e = geom.n4e;
                n4Db = geom.n4Db;
                n4Nb = geom.n4Nb;
                n4Sb = geom.n4Sb;
            end
            if ~iscell(n4e)
                n4e = mat2cell(n4e,ones(size(n4e,1),1));
            end
            if explicit && nargin < 3
                % all boundary edges (edges of one element) are Dirichlet
                edges = cell2mat(cellfun(@(z)[z(:),circshift(z(:),-1)],n4e(:),'UniformOutput',false));
                [~,~,j] = unique(sort(edges,2),'rows');
                once = accumarray(j,1) == 1;
                n4Db = edges(once(j),:);
            end
            if explicit && nargin < 4
                n4Nb = [];
            end
            if ~isempty(n4Sb) && size(n4Sb,2) ~= 3
                error('fem4m:MeshPoly:invalidN4Sb', ...
                    'n4Sb must have the rows [n1 n2 axis].');
            end
            obj.c4n = c4n;
            obj.n4e = n4e;
            obj.n4Db = n4Db;
            obj.n4Nb = n4Nb;
            obj.n4Sb = n4Sb;
            obj = getDataStructures(obj);
        end

        function length4ed = length4ed(obj)
            %length4ed - Compute the lengths of all edges
            %   LENGTH4ED = length4ed(OBJ) computes the Euclidean length of
            %   every edge listed in the node pairs OBJ.n4ed.
            length4ed = sqrt(sum((obj.c4n(obj.n4ed(:,2),:)-obj.c4n(obj.n4ed(:,1),:)).^2,2));
        end

        function point4ed = point4ed(obj,refp)
            %point4ed - Compute points on all edges
            %   POINT4ED = point4ed(OBJ,REFP) computes for every edge in
            %   OBJ.n4ed the point given by the reference coordinate REFP
            %   (a scalar, only REFP(1) is used), where 0 and 1 are its
            %   first and second node. POINT4ED is nrEdges-by-2.
            P1 = obj.c4n(obj.n4ed(:,1),:);
            P2 = obj.c4n(obj.n4ed(:,2),:);
            point4ed = P1 + (P2-P1)*refp(1);
        end

        function normal4ed = normal4ed(obj)
            %normal4ed - Compute the unit normals of all edges
            %   NORMAL4ED = normal4ed(OBJ) computes the unit normal vector
            %   [t2,-t1] of every edge by rotating its unit tangent t from
            %   OBJ.tangent4ed clockwise by a quarter turn.
            tangents = obj.tangent4ed;
            normal4ed = [tangents(:,2), -tangents(:,1)];
        end

        function tangent4ed = tangent4ed(obj)
            %tangent4ed - Compute the unit tangents of all edges
            %   TANGENT4ED = tangent4ed(OBJ) computes the unit tangent
            %   vector of every edge pointing from its first to its second
            %   node in OBJ.n4ed.
            c4start = obj.c4n(obj.n4ed(:,1),:);
            c4end = obj.c4n(obj.n4ed(:,2),:);
            tangent4ed = (c4end-c4start)./[obj.length4ed obj.length4ed];
        end

        function tangent4e = tangent4e(obj)
            %tangent4e - Compute the unit tangents of element edges
            %   TANGENT4E = tangent4e(OBJ) computes the unit tangent
            %   vectors of the edges of every element, oriented along the
            %   node order of the element, and caches the result.
            %   TANGENT4E is an nrElems-by-1 cell array whose entry for an
            %   element with J nodes is a 1-by-J-by-2 array; its K-th
            %   tangent belongs to the edge from node K to node K+1.
            tangent4e = obj.tangent4eSave;
            if isempty(tangent4e)
                tangent4e = cell(obj.nrElems,1);
                for j = unique(obj.nrN4e)'
                    parts = find(obj.nrN4e==j);
                    curN4e = cell2mat(obj.n4e(parts));
                    tangents = zeros(numel(parts),j,2);
                    for k = 1 : j
                        allEdges = curN4e(:,[k,mod(k,j)+1]);
                        c4start  = obj.c4n(allEdges(:,1),:);
                        c4end    = obj.c4n(allEdges(:,2),:);
                        lengths  = sqrt(sum((c4end-c4start).^2,2));
                        tangents(:,k,:) = (c4end - c4start)./[lengths lengths];
                    end
                    for curElem = 1:numel(parts)
                        tangent4e{parts(curElem)} = tangents(curElem,:,:);
                    end
                end
                obj.tangent4eSave = tangent4e;
            end
        end

        function normal4e = normal4e(obj)
            %normal4e - Compute the unit normals of element edges
            %   NORMAL4E = normal4e(OBJ) computes the unit normal vectors
            %   of the edges of every element by rotating the tangents of
            %   tangent4e clockwise by a quarter turn (outer normals for
            %   counterclockwise elements), and caches the result.
            %   NORMAL4E has the cell layout of tangent4e.
            normal4e = obj.normal4eSave;
            if isempty(normal4e)
                t4e = obj.tangent4e;
                normal4e = cell(obj.nrElems,1);
                for j = unique(obj.nrN4e)'
                    parts = find(obj.nrN4e==j);
                    tangents = cell2mat(t4e(parts));
                    normals = zeros(numel(parts),j,2);
                    for k = 1 : j
                        normals(:,k,:) = [tangents(:,k,2), -tangents(:,k,1)];
                    end
                    for curElem = 1:numel(parts)
                        normal4e{parts(curElem)} = normals(curElem,:,:);
                    end
                end
                obj.normal4eSave = normal4e;
            end
        end

        function edNr4ed = edNr4ed(obj)
            %edNr4ed - Compute local edge numbers in elements
            %   EDNR4ED = edNr4ed(OBJ) computes for every edge the number
            %   of its position within each of its adjacent elements, and
            %   caches the result. EDNR4ED is nrEdges-by-2; EDNR4ED(ed,k)
            %   is the local number of the edge ed in the element
            %   e4ed(ed,k), 0 for the missing second element of a
            %   boundary edge.
            edNr4ed = obj.edNr4edSave;
            if isempty(edNr4ed)
                edNr4ed = zeros(size(obj.e4ed));
                for j = unique(obj.nrN4e)'
                    parts1 = find(obj.nrN4e(obj.e4ed(:,1))==j);
                    parts2 = find(obj.e4ed(:,2)>0);
                    parts2 = parts2(obj.nrN4e(obj.e4ed(parts2,2))==j);
                    ind1 = obj.e4ed(parts1,1);
                    ind2 = obj.e4ed(parts2,2);
                    curEd4e1 = cell2mat(obj.ed4e(ind1));
                    curEd4e2 = cell2mat(obj.ed4e(ind2));
                    if ~isempty(parts1)
                        edNr4ed(parts1,1) = mod(find(curEd4e1'==(parts1*ones(1,j))')-1,j)+1;
                    end
                    if ~isempty(parts2)
                        edNr4ed(parts2,2) = mod(find(curEd4e2'==(parts2*ones(1,j))')-1,j)+1;
                    end
                end
                obj.edNr4edSave = edNr4ed;
            end
        end

        function color4n = computeColor4n(obj,levels)
            %computeColor4n - Color nodes of neighboring elements
            %   COLOR4N = computeColor4n(OBJ) colors every node with a
            %   small integer such that no two nodes sharing an element
            %   receive the same color.
            %
            %   COLOR4N = computeColor4n(OBJ,LEVELS) also widens the
            %   neighborhood to LEVELS layers of adjacent elements
            %   (default 1).
            if nargin < 2 || isempty(levels)
                levels = 1;
            end
            color4n = zeros(obj.nrNodes,1);
            e4nT = obj.e4n';
            for curNode = 1 : obj.nrNodes
                neighbours = curNode;
                for l = 1 : levels
                    elems = nonzeros(e4nT(:,neighbours));
                    indx = [0;cumsum(obj.nrN4e(elems))];
                    neighbours = zeros(indx(end),1);
                    for k = 1:numel(elems)
                        neighbours(indx(k)+1:indx(k+1)) = obj.n4e{elems(k)};
                    end
                end
                neighbour_colors = color4n(neighbours);
                color4n(curNode) = min(setdiff(1:max(neighbour_colors)+1,neighbour_colors));
            end
        end

        function area4e = area4e(obj)
            %area4e - Compute the areas of all elements
            %   AREA4E = area4e(OBJ) computes the area of every polygonal
            %   element with the shoelace formula and caches the result.
            area4e = obj.area4eSave;
            if isempty(area4e)
                xCoord = obj.c4n(:,1);
                yCoord = obj.c4n(:,2);
                area4e = zeros(size(obj.n4e,1),1);
                for j = unique(obj.nrN4e)'
                    curN4e = cell2mat(obj.n4e(obj.nrN4e==j));
                    x4e = xCoord(curN4e')';
                    y4e = yCoord(curN4e')';
                    area4e(obj.nrN4e==j) = abs(1/2*sum(x4e(:,1:end).*y4e(:,[2:end,1])-x4e(:,[2:end,1]).*y4e(:,1:end),2));
                end
                obj.area4eSave = area4e;
            end
        end

        function diam4e = diam4e(obj)
            %diam4e - Compute the diameters of all elements
            %   DIAM4E = diam4e(OBJ) computes the diameter of every
            %   element as the largest distance between two of its nodes,
            %   and caches the result.
            diam4e = obj.diam4eSave;
            if isempty(diam4e)
                diam4e = zeros(size(obj.n4e,1),1);
                for j = unique(obj.nrN4e)'
                    curN4e = cell2mat(obj.n4e(obj.nrN4e==j));
                    v = reshape(obj.c4n(curN4e,:),size(curN4e,1),j,1,2);
                    diam4e(obj.nrN4e==j) = max(vecnorm(v - permute(v,[1 3 2 4]),2,4),[],[2 3]);
                end
                obj.diam4eSave = diam4e;
            end
        end

        function mid4e = mid4e(obj)
            %mid4e - Compute the centroids of all elements
            %   MID4E = mid4e(OBJ) computes the centroid of every
            %   polygonal element as the area-weighted center obtained
            %   with the shoelace formula, and caches the result.
            mid4e = obj.mid4eSave;
            if isempty(mid4e)
                xCoord = obj.c4n(:,1);
                yCoord = obj.c4n(:,2);
                mid4e = zeros(size(obj.n4e,1),2);
                for j = unique(obj.nrN4e)'
                    curN4e = cell2mat(obj.n4e(obj.nrN4e==j));
                    x4e = xCoord(curN4e')';
                    y4e = yCoord(curN4e')';
                    A = 1/2*sum(x4e(:,1:end).*y4e(:,[2:end,1])-x4e(:,[2:end,1]).*y4e(:,1:end),2);
                    mid4e(obj.nrN4e==j,1) = sum( (x4e(:,1:end)+x4e(:,[2:end,1])).*(x4e(:,1:end).*y4e(:,[2:end,1])-x4e(:,[2:end,1]).*y4e(:,1:end)),2)./(6*A);
                    mid4e(obj.nrN4e==j,2) = sum( (y4e(:,1:end)+y4e(:,[2:end,1])).*(x4e(:,1:end).*y4e(:,[2:end,1])-x4e(:,[2:end,1]).*y4e(:,1:end)),2)./(6*A);
                end
                obj.mid4eSave = mid4e;
            end
        end

        function disp(obj)
            %DISP - Display summary or full data of the mesh
            %   DISP(OBJ) displays the node coordinates, element node
            %   lists, and boundary nodes for meshes with at most 100
            %   nodes and only the numbers of nodes, elements, and edges
            %   otherwise.
            if obj.nrNodes <= 100
                disp('coordinates for nodes:')
                disp(obj.c4n)
                disp('nodes for elements:')
                celldisp(obj.n4e,'n4e')
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

        function plot(obj)
            %PLOT - Plot the polygonal mesh
            %   PLOT(OBJ) draws every polygonal element as a white patch
            %   with black edges into the current axes and leaves hold
            %   on.
            X = cellfun(@(curElem)(obj.c4n(curElem,1)),obj.n4e,'UniformOutput',false)';
            Y = cellfun(@(curElem)(obj.c4n(curElem,2)),obj.n4e,'UniformOutput',false)';
            for j = min(obj.nrN4e) : max(obj.nrN4e)
                curX = cell2mat(X(1,obj.nrN4e==j));
                curY = cell2mat(Y(1,obj.nrN4e==j));
                patch(curX,curY,'w','Edgecolor','k');
                hold on;
            end
        end

        function plotBd(obj,color)
            %plotBd - Plot the boundary of the mesh
            %   plotBd(OBJ) draws the Dirichlet boundary edges in red, the
            %   Neumann boundary edges in green, and the symmetry boundary
            %   edges dotted in blue as patches.
            %
            %   plotBd(OBJ,COLOR) also uses the color COLOR for all three
            %   boundaries.
            X = obj.c4n(:,1);
            Y = obj.c4n(:,2);
            if nargin < 2
                color = 'red';
            end
            if ~isempty(obj.n4Db)
                patch(X(obj.n4Db)',Y(obj.n4Db)','w','Edgecolor',color);
            end
            if nargin < 2
                color = 'green';
            end
            if ~isempty(obj.n4Nb)
                patch(X(obj.n4Nb)',Y(obj.n4Nb)','w','Edgecolor',color);
            end
            if nargin < 2
                color = 'blue';
            end
            if ~isempty(obj.n4Sb)
                patch(X(obj.n4Sb(:,[1 2]))',Y(obj.n4Sb(:,[1 2]))','w','Edgecolor',color,'LineStyle',':');
            end
        end

        function plotEnumeration(obj,fontsize)
            %plotEnumeration - Label nodes, elements, and edges with ids
            %   plotEnumeration(OBJ) writes the ids of the nodes, of the
            %   elements at their points mid4e (centroids; vertex averages
            %   in MeshPolyQ), and of the edges at their midpoints in
            %   different colors into the current axes; it does not draw
            %   the mesh.
            %
            %   plotEnumeration(OBJ,FONTSIZE) also specifies the font size
            %   of the labels, which defaults to 24.
            if nargin < 2 || isempty(fontsize)
                fontsize = 24;
            end
            mid4e = obj.mid4e();
            mid4ed = obj.point4ed(1/2);
            text(obj.c4n(:,1),obj.c4n(:,2),num2str((1:size(obj.c4n,1))'),...
                'color',[0.26 0.49 0.19],'FontWeight','bold','FontSize',fontsize,'HorizontalAlignment','center');
            text(mid4e(:,1),mid4e(:,2),num2str((1:size(obj.n4e,1))'),...
                'color',[0.95 0.91 0],'FontWeight','bold','FontSize',fontsize,'HorizontalAlignment','center');
            text(mid4ed(:,1),mid4ed(:,2),num2str((1:obj.nrEdges)'),...
                'color',[0.38 0.02 0.68],'FontWeight','bold','FontSize',fontsize,'HorizontalAlignment','center');
        end

        function plotColor4n(obj,color)
            %plotColor4n - Plot the mesh colored by node colors
            %   plotColor4n(OBJ) plots the mesh and marks every node with
            %   a color from the greedy coloring of computeColor4n.
            %
            %   plotColor4n(OBJ,COLOR) also uses the given color vector
            %   COLOR (one value per node) instead of computing one.
            plot(obj);hold on;
            if nargin < 2 || isempty(color)
                color = obj.computeColor4n;
            end
            scatter(obj.c4n(:,1),obj.c4n(:,2),200,color,'filled');
            colorbar;
        end

        function obj = plus(obj1,obj2)
            %PLUS - Union of two polygon meshes
            %   OBJ3 = PLUS(OBJ1,OBJ2), or OBJ3 = OBJ1 + OBJ2, merges two
            %   polygon meshes into a single MeshPoly (also for subclass
            %   inputs) by identifying nodes with identical coordinates
            %   and removing all ghost nodes and the boundary edges shared
            %   by both meshes, also those with different boundary
            %   conditions in OBJ1 and OBJ2.
            c4nPlus = [obj1.c4n;obj2.c4n];
            n4ePlus = [obj1.n4e;obj2.n4e];
            for k = obj1.nrElems+(1:obj2.nrElems)
                n4ePlus{k} = obj1.nrNodes+n4ePlus{k};
            end
            n4DbPlus = [obj1.n4Db;obj1.nrNodes+obj2.n4Db];
            n4NbPlus = [obj1.n4Nb;obj1.nrNodes+obj2.n4Nb];
            n4SbPlus = obj2.n4Sb;
            if ~isempty(n4SbPlus)
                n4SbPlus(:,1:2) = obj1.nrNodes+n4SbPlus(:,1:2);
            end
            n4SbPlus = [obj1.n4Sb;n4SbPlus];
            % find double nodes
            [~,i1,i2] = intersect(obj1.c4n,obj2.c4n,'rows');
            for k = 1 : numel(i1)
                for j = 1 : size(n4ePlus,1)
                    curN4e = n4ePlus{j};
                    curN4e(curN4e==obj1.nrNodes+i2(k)) = i1(k);
                    n4ePlus{j} = curN4e;
                end
                n4DbPlus(n4DbPlus==obj1.nrNodes+i2(k)) = i1(k);
                n4NbPlus(n4NbPlus==obj1.nrNodes+i2(k)) = i1(k);
                if ~isempty(n4SbPlus)
                    sbNodes = n4SbPlus(:,1:2);
                    sbNodes(sbNodes==obj1.nrNodes+i2(k)) = i1(k);
                    n4SbPlus(:,1:2) = sbNodes;
                end
            end
            %double boundary edges, also from different boundary lists
            n4bd = [n4DbPlus;n4NbPlus];
            if ~isempty(n4SbPlus)
                n4bd = [n4bd;n4SbPlus(:,1:2)];
            end
            interior = ismember(n4bd,n4bd(:,[2 1]),'rows');
            nrDb = size(n4DbPlus,1);
            nrNb = size(n4NbPlus,1);
            n4DbPlus = n4DbPlus(~interior(1:nrDb),:);
            n4NbPlus = n4NbPlus(~interior(nrDb+(1:nrNb)),:);
            n4SbPlus = n4SbPlus(~interior(nrDb+nrNb+1:end),:);
            % remove additional nodes
            [c4nPlus,n4ePlus,n4DbPlus,n4NbPlus,n4SbPlus] = MeshPoly.removeGhostNodes(c4nPlus,n4ePlus,n4DbPlus,n4NbPlus,n4SbPlus);
            %
            obj = MeshPoly(c4nPlus,n4ePlus,n4DbPlus,n4NbPlus,n4SbPlus);
        end

        function obj = refine(obj,~) %#ok<INUSD>
            %REFINE - Refine marked elements
            %   OBJ = REFINE(OBJ,R4E) raises an error and ignores the
            %   marked elements R4E: MeshPoly does not implement
            %   refinement, which depends on the element shape of a
            %   subclass (MeshPolyT, MeshPolyQ, MeshPolyH).
            error('fem4m:MeshPoly:refineNotImplemented',' refinement routine not implemented');
        end

        function obj = refineGraded(obj,mu,H,gradPoints,R)
            %refineGraded - Refine the mesh with a graded distribution
            %   OBJ = refineGraded(OBJ,MU,H) repeatedly refines the
            %   elements with refine of the subclass until every element
            %   T has a diameter h_T of at most H and of at most
            %   H*(r_T/R(k))^(1-MU(k)) for every grading point k, where
            %   r_T is the distance of the point mid4e of T to the grading
            %   point (h_T/2 if it is zero).
            %
            %   OBJ = refineGraded(...,gradPoints) also specifies the
            %   grading points, one per row, which default to [0,0].
            %
            %   OBJ = refineGraded(...,R) also specifies the radii R(k) of
            %   the regions graded by the k-th grading point (default 1).
            %   Scalar MU and R hold for all grading points.
            if nargin < 4 || isempty(gradPoints)
                gradPoints = [0,0];
            end
            if nargin < 5 || isempty(R)
                R = 1;
            end
            % one exponent and radius per grading point
            mu = mu(:).*ones(size(gradPoints,1),1);
            R = R(:).*ones(size(gradPoints,1),1);

            count = 1;
            while count ~= 0
                midPoint4e = obj.mid4e;
                refineElems = false(obj.nrElems,1);
                hT = obj.diam4e;
                refineElems(hT > H) = true;
                for k = 1 : size(gradPoints,1)
                    rT = sqrt( (midPoint4e(:,1)-gradPoints(k,1)).^2+...
                        (midPoint4e(:,2)-gradPoints(k,2)).^2  );
                    rT(rT<eps) = hT(rT<eps)./2;
                    refineElems(hT > H.*(rT/R(k)).^(1-mu(k))) = true;
                end
                count = nnz(refineElems);
                obj = obj.refine(refineElems);
            end
        end

        function [mesh,parent4ed,children4ed] = subtriangulateMid(obj)
            %subtriangulateMid - Split polygons into triangles around mid4e
            %   MESH = subtriangulateMid(OBJ) splits every polygonal
            %   element into triangles spanned by two consecutive nodes
            %   and the point mid4e of the element (its centroid; the
            %   vertex average of the quad in MeshPolyQ), and returns them
            %   as a Mesh. Elements that are not star-shaped with respect
            %   to mid4e, such as the element at the tip of a slit, are
            %   split into the triangles spanned by a vertex they are
            %   star-shaped to and the remaining edges; if there is no
            %   such vertex, it raises an error. MESH.parent4e holds the
            %   element of OBJ of every triangle.
            %
            %   [MESH,PARENT4ED] = subtriangulateMid(OBJ) also returns the
            %   column vector PARENT4ED that holds for every edge of MESH
            %   the id of the edge of OBJ it lies on, 0 for the new inner
            %   edges.
            %
            %   [...,CHILDREN4ED] = subtriangulateMid(OBJ) also returns
            %   the ids of the edges of MESH that lie on edges of OBJ,
            %   ordered by their parent edge; every edge of OBJ is one
            %   edge of MESH, so CHILDREN4ED(k) is the edge of MESH on the
            %   edge k of OBJ.
            c4nNew = [obj.c4n;obj.mid4e];
            n4eNew = [];
            parent4e = [];
            for j = unique(obj.nrN4e)'
                curIndx = find(obj.nrN4e==j);
                curParts = cell2mat(obj.n4e(curIndx));
                for k = 1:j
                    n4eNew = [n4eNew;curParts(:,[k mod(k,j)+1]), obj.nrNodes+curIndx]; %#ok<AGROW>
                    parent4e = [parent4e;curIndx]; %#ok<AGROW>
                end
            end
            % polygons that are not star-shaped with respect to their
            % mid4e, e.g. at the tip of a slit, are split from a vertex
            % they are star-shaped to
            area = MeshPoly.signedArea(c4nNew(n4eNew(:,1),:),c4nNew(n4eNew(:,2),:),c4nNew(n4eNew(:,3),:));
            bad = find(accumarray(parent4e,area <= 0,[obj.nrElems,1]));
            if ~isempty(bad)
                remove = ismember(parent4e,bad);
                n4eNew(remove,:) = [];
                parent4e(remove) = [];
                for curElem = bad'
                    nodes = obj.n4e{curElem}(:)';
                    j = numel(nodes);
                    next = [2:j,1];
                    tri = [];
                    for k = 1 : j
                        edges = setdiff(1:j,[k,mod(k-2,j)+1]);
                        if all(MeshPoly.signedArea(obj.c4n(nodes(edges),:), ...
                                obj.c4n(nodes(next(edges)),:),obj.c4n(nodes(k),:)) > 0)
                            tri = [nodes(edges)',nodes(next(edges))',nodes(k)*ones(numel(edges),1)];
                            break
                        end
                    end
                    if isempty(tri)
                        error('fem4m:MeshPoly:notStarShaped', ...
                            'Polygon %d is not star-shaped with respect to its centroid or a vertex',curElem);
                    end
                    n4eNew = [n4eNew;tri]; %#ok<AGROW>
                    parent4e = [parent4e;curElem*ones(size(tri,1),1)]; %#ok<AGROW>
                end
                % remove the unused centroids
                keep = true(size(c4nNew,1),1);
                keep(obj.nrNodes+bad) = false;
                newNode = cumsum(keep);
                n4eNew = newNode(n4eNew);
                c4nNew = c4nNew(keep,:);
            end
            mesh = Mesh(c4nNew,n4eNew,obj.n4Db,obj.n4Nb,obj.n4Sb);
            mesh.parent4e = parent4e;
            sn4ed = sparse(obj.n4ed(:,1),obj.n4ed(:,2),1:obj.nrEdges,mesh.nrNodes,mesh.nrNodes);
            parent4ed = rowaddr(sn4ed,mesh.n4ed(:,1),mesh.n4ed(:,2));
            [children4ed,~,ed] = find(parent4ed);
            [~,perm] = sort(ed);
            children4ed = children4ed(perm);
        end

        function obj = save(obj,filename)
            %SAVE - Write the mesh data to .dat files
            %   SAVE(OBJ) writes the node coordinates (all digits), element
            %   node lists, and boundary node pairs to text files named
            %   mesh_*.dat, which MeshPoly('mesh') reads again. An empty
            %   boundary list writes no file and deletes the file of an
            %   earlier save.
            %
            %   SAVE(OBJ,FILENAME) also uses FILENAME as the prefix of the
            %   file names instead of 'mesh'.
            if nargin < 2 || isempty(filename)
                filename = 'mesh';
            end
            % all 17 digits, writematrix alone writes 15
            writematrix(compose("%.17g",obj.c4n),[filename,'_c4n.dat'], ...
                'Delimiter',' ','QuoteStrings','none');
            writecell(obj.n4e,[filename,'_n4e.dat'],'Delimiter','tab');
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

    end

    methods (Access = protected)

        function [sub,hnodes] = subMesh(obj,name,k,file)
            %subMesh - Underlying mesh of a subclass and its hanging nodes
            %   [SUB,HNODES] = subMesh(OBJ,NAME,K,FILE) reads the elements
            %   SUB (K vertices each) of the underlying mesh and the
            %   hanging nodes HNODES from the files NAME_FILE.dat and
            %   NAME_hnodes.dat written by save, if NAME is a file name
            %   prefix with such files for this mesh (one row of FILE per
            %   element). Otherwise every polygon must have K vertices and
            %   is an element of SUB, and HNODES is the all-zero sparse
            %   nrNodes-by-nrNodes matrix.
            nN = size(obj.c4n,1);
            if ischar(name) && exist([name,'_',file,'.dat'],'file')
                sub = load([name,'_',file,'.dat']);
                if size(sub,1) == numel(obj.n4e)
                    ijv = zeros(0,3);
                    if exist([name,'_hnodes.dat'],'file')
                        ijv = [ijv;load([name,'_hnodes.dat'])];
                    end
                    hnodes = sparse(ijv(:,1),ijv(:,2),ijv(:,3),nN,nN);
                    return
                end
            end
            if any(cellfun(@numel,obj.n4e) ~= k)
                error('fem4m:MeshPoly:vertexCount', ...
                    ['Every element of %s needs %d vertices; meshes with ', ...
                    'hanging nodes are read from the files of save.'],class(obj),k);
            end
            sub = cell2mat(obj.n4e);
            hnodes = sparse(nN,nN);
        end

    end

    methods (Access = private)

        function obj = getDataStructures(obj)
            %getDataStructures - Compute all derived mesh data structures
            obj.nrNodes = size(obj.c4n,1);
            obj.nrElems = size(obj.n4e,1);
            obj.nrN4e = cellfun(@numel,obj.n4e);
            twoN4ed = zeros(sum(obj.nrN4e),2);
            indpoly = [0;cumsum(obj.nrN4e)];
            for curElem = 1 : size(obj.n4e,1)
                elem = obj.n4e{curElem};
                twoN4ed(indpoly(curElem)+1:indpoly(curElem+1),:) = elem([1:numel(elem);2:numel(elem),1])';
            end
            if ~isempty(obj.n4Db)
                twoN4ed = [twoN4ed;obj.n4Db(:,[2 1])];
            end
            if ~isempty(obj.n4Nb)
                twoN4ed = [twoN4ed;obj.n4Nb(:,[2 1])];
            end
            if ~isempty(obj.n4Sb)
                twoN4ed = [twoN4ed;obj.n4Sb(:,[2 1])];
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
            obj.ed4e = mat2cell(twoEd(1:sum(obj.nrN4e))',1,obj.nrN4e)';
            obj.ed4Db = twoEd(sum(obj.nrN4e)+1:sum(obj.nrN4e)+size(obj.n4Db,1));
            obj.ed4Nb = twoEd(sum(obj.nrN4e)+size(obj.n4Db,1)+1:sum(obj.nrN4e)+size(obj.n4Db,1)+size(obj.n4Nb,1));
            obj.ed4Sb = twoEd(sum(obj.nrN4e)+size(obj.n4Db,1)+size(obj.n4Nb,1)+1:end);
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
            elems = zeros(sum(obj.nrN4e),1);
            for curElem = 1 : size(obj.n4e,1)
                elems(indpoly(curElem)+1:indpoly(curElem+1)) = curElem;
            end
            indx = indx(indx <= sum(obj.nrN4e));
            indy = indy(indy <= sum(obj.nrN4e));
            obj.e4ed(twoEd(indx),1) = elems(indx);
            obj.e4ed(twoEd(indy),2) = elems(indy);
            indz = find(obj.e4ed(:,1)==0);
            obj.e4ed(indz,1) = obj.e4ed(indz,2);
            obj.e4ed(indz,2) = 0;
            obj.e4n = sparse(twoN4ed(1:sum(obj.nrN4e),1),elems,elems,obj.nrNodes,obj.nrElems);
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
        end

    end

    methods(Static)
        function a = signedArea(P,Q,c)
            %signedArea - Twice the signed areas of triangles
            %   A = signedArea(P,Q,C) returns twice the signed areas of the
            %   triangles with the vertices given by the rows of P, Q and
            %   C, positive for counterclockwise triangles.
            a = (Q(:,1)-P(:,1)).*(c(:,2)-P(:,2)) - (Q(:,2)-P(:,2)).*(c(:,1)-P(:,1));
        end

        function [c4n,n4e,n4Db,n4Nb,n4Sb] = removeGhostNodes(c4n,n4e,n4Db,n4Nb,n4Sb)
            %removeGhostNodes - Remove unused nodes from mesh data
            %   [c4n,n4e] = removeGhostNodes(c4n,n4e) removes all nodes
            %   that are not used by any element of the cell array n4e and
            %   renumbers the remaining ones in n4e.
            %
            %   [c4n,n4e,n4Db] = removeGhostNodes(...,n4Db) also renumbers
            %   the Dirichlet boundary node pairs n4Db accordingly.
            %
            %   [...,n4Nb] = removeGhostNodes(...,n4Nb) also renumbers the
            %   Neumann boundary node pairs n4Nb accordingly.
            %
            %   [...,n4Sb] = removeGhostNodes(...,n4Sb) also renumbers the
            %   symmetry boundary node pairs in the first two columns of
            %   n4Sb accordingly.
            if nargin < 3
                n4Db = [];
            end
            if nargin < 4
                n4Nb = [];
            end
            if nargin < 5
                n4Sb = [];
            end
            rN = setdiff(1:size(c4n,1),cell2mat(cellfun(@transpose,n4e,'UniformOutput',false)));
            c4n(rN,:) = [];
            k = 1;
            while  k <= numel(rN)
                m = 1;
                while( k+m <= numel(rN) && rN(k+m)==rN(k)+m)
                    m = m+1;
                end
                for j = 1 : size(n4e,1)
                    curN4e = n4e{j};
                    for l = 1 : numel(curN4e)
                        if curN4e(l) > rN(k)
                            curN4e(l) = curN4e(l) - m;
                        end
                    end
                    n4e{j} = curN4e;
                end
                n4Db(n4Db(:)>rN(k))=n4Db(n4Db(:)>rN(k))-m;
                n4Nb(n4Nb(:)>rN(k))=n4Nb(n4Nb(:)>rN(k))-m;
                if ~isempty(n4Sb)
                    sbNodes = n4Sb(:,1:2);
                    sbNodes(sbNodes>rN(k)) = sbNodes(sbNodes>rN(k))-m;
                    n4Sb(:,1:2) = sbNodes;
                end
                rN = rN - m;
                k = k + m;
            end
        end
    end

end