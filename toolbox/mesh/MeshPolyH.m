classdef MeshPolyH < MeshPoly
    %MeshPolyH - Polygon mesh with hanging-node-free H-refinement
    %   OBJ = MeshPolyH(GEOM) creates a polygon mesh of H type from the
    %   struct or Geometry object GEOM (see MeshPoly), built with the
    %   type 'H' (the coarse H-type meshes of Geometry for the square,
    %   which 'Circle' also uses, the L-shape and the slit domain).
    %
    %   OBJ = MeshPolyH(NAME) reads the files NAMEH_*.dat if
    %   NAMEH_c4n.dat exists and otherwise NAME_*.dat, e.g. written by
    %   save.
    %
    %   OBJ = MeshPolyH(c4n,n4e) creates the mesh with the node
    %   coordinates c4n and the node ids n4e of the elements (a cell
    %   array or a matrix); all boundary edges are Dirichlet edges.
    %
    %   OBJ = MeshPolyH(...,n4Db) also specifies the node pairs of the
    %   Dirichlet boundary.
    %
    %   OBJ = MeshPolyH(...,n4Nb) also specifies the node pairs of the
    %   Neumann boundary.
    %
    %   OBJ = MeshPolyH(...,n4Sb) also specifies the node pairs of the
    %   symmetry boundary, whose third column is the fixed axis (1 or
    %   2, 0 for a general normal).
    %
    %   refine splits marked polygons such that no hanging nodes arise;
    %   it keeps no refinement data (no parent4e or prolongation).
    %
    %   MeshPolyH functions:
    %       refine       - Refine marked elements of the polygon mesh
    %       closure      - Enforce the closure condition on marked elements
    %       joinElements - Join pairs of quadrilaterals into polygons
    %
    %   See also MeshPoly, MeshPolyT, MeshPolyQ, Geometry

    % Copyright (c) 2026 by the fem4m authors

    % refine assembles the new mesh block by block of unknown size
    %#ok<*AGROW>

    methods
        function obj = MeshPolyH(c4n,n4e,n4Db,n4Nb,n4Sb)
            %MeshPolyH - Create an H-type polygon mesh object
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
            if ischar(c4n) && exist([c4n,'H_c4n.dat'],'file')
                % files of the H type, NAMEH_c4n.dat
                c4n = [c4n,'H'];
            elseif isstruct(c4n)
                c4n.type = 'H';
            elseif strcmpi(class(c4n),'Geometry')
                if c4n.opts.type ~= 'H'
                    opts = c4n.opts;
                    opts.type = 'H';
                    c4n = Geometry(opts);
                end
            end

            args = {c4n,n4e,n4Db,n4Nb,n4Sb};
            obj = obj@MeshPoly(args{1:given});
        end

        function obj = refine(obj,r4e)
            %REFINE - Refine marked elements of the polygon mesh
            %   OBJ = REFINE(OBJ,R4E) refines the elements marked by the logical or 0/1
            %   vector R4E (one entry per element) and returns the new mesh, which has
            %   no hanging nodes. closure first moves the marks of quadrilaterals whose
            %   longest edge is an inner edge to their neighbors across that edge.
            %   Every marked polygon is then split into an inner polygon with the same
            %   number of nodes and quadrilaterals between it and the old edges: the
            %   new nodes lie halfway between a center and the vertices. The center is
            %   the centroid of a hexagon, the midpoint of the longest edge of a
            %   quadrilateral or pentagon (that edge is divided into three parts), or
            %   the corner node of a polygon with 7 or 8 nodes at the reentrant corner
            %   of the L-shape or the slit (its two edges at that node are bisected). A
            %   quadrilateral with two adjacent boundary edges is split into three
            %   quadrilaterals by the midpoints of the segments from their common node
            %   to the other three vertices; two opposite boundary edges raise the
            %   error fem4m:MeshPolyH:cornerNotImplemented. Marked elements with
            %   another number of nodes stay as they are. The boundary lists follow the
            %   new edges, and joinElements finally merges pairs of quadrilaterals.

            r4e = obj.closure(r4e); %closure

            newN4e = obj.n4e(~r4e);
            newC4n = obj.c4n;
            newDb = obj.n4Db(~r4e(obj.e4ed(obj.ed4Db,1)),:);
            newNb = obj.n4Nb(~r4e(obj.e4ed(obj.ed4Nb,1)),:);
            % symmetry boundary with the fixed axis in the third column
            n4Sb = obj.n4Sb;
            if isempty(n4Sb)
                n4Sb = zeros(0,3);
            end
            newSb = n4Sb(~r4e(obj.e4ed(obj.ed4Sb,1)),:);

            e4ed2 = obj.e4ed(:,2);
            length4ed = obj.length4ed;

            for j = unique(obj.nrN4e)'

                parts = find(obj.nrN4e==j & r4e);
                if ~isempty(parts)
                    curN4e = cell2mat(obj.n4e(parts));
                    curEd4e = cell2mat(obj.ed4e(parts));
                    % rotate n4e and ed4e: the longest edge comes first
                    if size(curEd4e,1)==1
                        length = length4ed(curEd4e)';
                    else
                        length = length4ed(curEd4e);
                    end
                    [~,indx] = max(length,[],2);
                    for k = 2 : j
                        curN4e(indx==k,:) = curN4e(indx==k,mod((k:k+j-1)-1,j)+1);
                        curEd4e(indx==k,:) = curEd4e(indx==k,mod((k:k+j-1)-1,j)+1);
                    end
                    % refine
                    if j == 6
                        mid4e = obj.mid4e;
                        mid4e = mid4e(parts,:);
                        newNode = zeros(numel(parts),6);
                        for k = 1 : 6
                            newNode(:,k) = size(newC4n,1)+(1:numel(parts));
                            newC4n = [newC4n;(mid4e+obj.c4n(curN4e(:,k),:))/2];
                        end
                        newN4e = [newN4e;
                            mat2cell(newNode,ones(size(newNode,1),1))];
                        curNewN4e = [curN4e(:,[1 2]),newNode(:,[2 1]);
                            curN4e(:,[2 3]),newNode(:,[3 2]);
                            curN4e(:,[3 4]),newNode(:,[4 3]);
                            curN4e(:,[4 5]),newNode(:,[5 4]);
                            curN4e(:,[5 6]),newNode(:,[6 5]);
                            curN4e(:,[6 1]),newNode(:,[1 6])];
                        newN4e = [newN4e;
                            mat2cell(curNewN4e,ones(size(curNewN4e,1),1))];
                        % boundary
                        if size(curEd4e,1)==1
                            bdryE = find(any(e4ed2(curEd4e)'==0,2));
                        else
                            bdryE = find(any(e4ed2(curEd4e)==0,2));
                        end
                        for k = 1 : 6
                            [~,indx] = intersect(curEd4e(bdryE,k),obj.ed4Db);
                            indx = bdryE(indx);
                            newDb = [newDb;curN4e(indx,k),curN4e(indx,mod(k,6)+1)];
                            [~,indx] = intersect(curEd4e(bdryE,k),obj.ed4Nb);
                            indx = bdryE(indx);
                            newNb = [newNb;curN4e(indx,k),curN4e(indx,mod(k,6)+1)];
                            [~,indx,indSb] = intersect(curEd4e(bdryE,k),obj.ed4Sb);
                            indx = bdryE(indx);
                            newSb = [newSb;curN4e(indx,k),curN4e(indx,mod(k,6)+1),n4Sb(indSb,3)];
                        end

                    elseif j == 4
                        % corner element?
                        if size(curEd4e,1)==1
                            corner = 2==sum(e4ed2(curEd4e)'==0,2);
                        else
                            corner = 2==sum(e4ed2(curEd4e)==0,2);
                        end
                        %corner
                        for k = find(corner)'
                            if e4ed2(curEd4e(k,1))==0 && e4ed2(curEd4e(k,2))==0 % right corner
                                newNode = size(newC4n,1)+(1:3);
                                newC4n = [newC4n;
                                    (obj.c4n(curN4e(k,2),:)+obj.c4n(curN4e(k,1),:))/2;
                                    (obj.c4n(curN4e(k,2),:)+obj.c4n(curN4e(k,3),:))/2;
                                    (obj.c4n(curN4e(k,2),:)+obj.c4n(curN4e(k,4),:))/2;];
                                curNewN4e = [curN4e(k,[4 1]),newNode([1 3]);
                                    newNode(1),curN4e(k,2),newNode([2 3]);
                                    curN4e(k,[3 4]),newNode([3 2])];
                                newN4e = [newN4e;
                                    mat2cell(curNewN4e,ones(size(curNewN4e,1),1))];
                                %
                                if ~isempty(intersect(curEd4e(k,1),obj.ed4Db))
                                    newDb = [newDb;
                                        curN4e(k,1),newNode(1);
                                        newNode(1),curN4e(k,2)];
                                elseif ~isempty(intersect(curEd4e(k,1),obj.ed4Sb))
                                    sbAxis = n4Sb(obj.ed4Sb==curEd4e(k,1),3);
                                    newSb = [newSb;
                                        curN4e(k,1),newNode(1),sbAxis;
                                        newNode(1),curN4e(k,2),sbAxis];
                                else
                                    newNb = [newNb;
                                        curN4e(k,1),newNode(1);
                                        newNode(1),curN4e(k,2)];
                                end
                                if ~isempty(intersect(curEd4e(k,2),obj.ed4Db))
                                    newDb = [newDb;
                                        curN4e(k,2),newNode(2);
                                        newNode(2),curN4e(k,3)];
                                elseif ~isempty(intersect(curEd4e(k,2),obj.ed4Sb))
                                    sbAxis = n4Sb(obj.ed4Sb==curEd4e(k,2),3);
                                    newSb = [newSb;
                                        curN4e(k,2),newNode(2),sbAxis;
                                        newNode(2),curN4e(k,3),sbAxis];
                                else
                                    newNb = [newNb;
                                        curN4e(k,2),newNode(2);
                                        newNode(2),curN4e(k,3)];
                                end
                            elseif  e4ed2(curEd4e(k,1))==0 && e4ed2(curEd4e(k,4))==0% left corner
                                newNode = size(newC4n,1)+(1:3);
                                newC4n = [newC4n;
                                    (obj.c4n(curN4e(k,1),:)+obj.c4n(curN4e(k,2),:))/2;
                                    (obj.c4n(curN4e(k,1),:)+obj.c4n(curN4e(k,3),:))/2;
                                    (obj.c4n(curN4e(k,1),:)+obj.c4n(curN4e(k,4),:))/2;];
                                curNewN4e = [curN4e(k,[2 3]),newNode([2 1]);
                                    curN4e(k,1),newNode([1 2 3]);
                                    curN4e(k,[3 4]),newNode([3 2])];
                                newN4e = [newN4e;
                                    mat2cell(curNewN4e,ones(size(curNewN4e,1),1))];
                                if ~isempty(intersect(curEd4e(k,1),obj.ed4Db))
                                    newDb = [newDb;
                                        curN4e(k,1),newNode(1);
                                        newNode(1),curN4e(k,2)];
                                elseif ~isempty(intersect(curEd4e(k,1),obj.ed4Sb))
                                    sbAxis = n4Sb(obj.ed4Sb==curEd4e(k,1),3);
                                    newSb = [newSb;
                                        curN4e(k,1),newNode(1),sbAxis;
                                        newNode(1),curN4e(k,2),sbAxis];
                                else
                                    newNb = [newNb;
                                        curN4e(k,1),newNode(1);
                                        newNode(1),curN4e(k,2)];
                                end
                                if ~isempty(intersect(curEd4e(k,4),obj.ed4Db))
                                    newDb = [newDb;
                                        curN4e(k,4),newNode(3);
                                        newNode(3),curN4e(k,1)];
                                elseif ~isempty(intersect(curEd4e(k,4),obj.ed4Sb))
                                    sbAxis = n4Sb(obj.ed4Sb==curEd4e(k,4),3);
                                    newSb = [newSb;
                                        curN4e(k,4),newNode(3),sbAxis;
                                        newNode(3),curN4e(k,1),sbAxis];
                                else
                                    newNb = [newNb;
                                        curN4e(k,4),newNode(3);
                                        newNode(3),curN4e(k,1)];
                                end
                            elseif  e4ed2(curEd4e(k,3))==0 && e4ed2(curEd4e(k,4))==0% left corner
                                newNode = size(newC4n,1)+(1:3);
                                newC4n = [newC4n;
                                    (obj.c4n(curN4e(k,4),:)+obj.c4n(curN4e(k,1),:))/2;
                                    (obj.c4n(curN4e(k,4),:)+obj.c4n(curN4e(k,2),:))/2;
                                    (obj.c4n(curN4e(k,4),:)+obj.c4n(curN4e(k,3),:))/2;];
                                curNewN4e = [curN4e(k,[2 3]),newNode([3 2]);
                                    curN4e(k,[1 2]),newNode([2 1]);
                                    curN4e(k,4),newNode([1 2 3])];
                                newN4e = [newN4e;
                                    mat2cell(curNewN4e,ones(size(curNewN4e,1),1))];
                                if ~isempty(intersect(curEd4e(k,3),obj.ed4Db))
                                    newDb = [newDb;
                                        curN4e(k,3),newNode(3);
                                        newNode(3),curN4e(k,4)];
                                elseif ~isempty(intersect(curEd4e(k,3),obj.ed4Sb))
                                    sbAxis = n4Sb(obj.ed4Sb==curEd4e(k,3),3);
                                    newSb = [newSb;
                                        curN4e(k,3),newNode(3),sbAxis;
                                        newNode(3),curN4e(k,4),sbAxis];
                                else
                                    newNb = [newNb;
                                        curN4e(k,3),newNode(3);
                                        newNode(3),curN4e(k,4)];
                                end
                                if ~isempty(intersect(curEd4e(k,4),obj.ed4Db))
                                    newDb = [newDb;
                                        curN4e(k,4),newNode(1);
                                        newNode(1),curN4e(k,1)];
                                elseif ~isempty(intersect(curEd4e(k,4),obj.ed4Sb))
                                    sbAxis = n4Sb(obj.ed4Sb==curEd4e(k,4),3);
                                    newSb = [newSb;
                                        curN4e(k,4),newNode(1),sbAxis;
                                        newNode(1),curN4e(k,1),sbAxis];
                                else
                                    newNb = [newNb;
                                        curN4e(k,4),newNode(1);
                                        newNode(1),curN4e(k,1)];
                                end
                            elseif  e4ed2(curEd4e(k,2))==0 && e4ed2(curEd4e(k,3))==0% right corner
                                newNode = size(newC4n,1)+(1:3);
                                newC4n = [newC4n;
                                    (obj.c4n(curN4e(k,3),:)+obj.c4n(curN4e(k,1),:))/2;
                                    (obj.c4n(curN4e(k,3),:)+obj.c4n(curN4e(k,2),:))/2;
                                    (obj.c4n(curN4e(k,3),:)+obj.c4n(curN4e(k,4),:))/2;];
                                curNewN4e = [curN4e(k,[1 2]),newNode([2 1]);
                                    newNode(2),curN4e(k,3),newNode([3 1]);
                                    curN4e(k,[4 1]),newNode([1 3])];
                                newN4e = [newN4e;
                                    mat2cell(curNewN4e,ones(size(curNewN4e,1),1))];
                                if ~isempty(intersect(curEd4e(k,2),obj.ed4Db))
                                    newDb = [newDb;
                                        curN4e(k,2),newNode(2);
                                        newNode(2),curN4e(k,3)];
                                elseif ~isempty(intersect(curEd4e(k,2),obj.ed4Sb))
                                    sbAxis = n4Sb(obj.ed4Sb==curEd4e(k,2),3);
                                    newSb = [newSb;
                                        curN4e(k,2),newNode(2),sbAxis;
                                        newNode(2),curN4e(k,3),sbAxis];
                                else
                                    newNb = [newNb;
                                        curN4e(k,2),newNode(2);
                                        newNode(2),curN4e(k,3)];
                                end
                                if ~isempty(intersect(curEd4e(k,3),obj.ed4Db))
                                    newDb = [newDb;
                                        curN4e(k,3),newNode(3);
                                        newNode(3),curN4e(k,4)];
                                elseif ~isempty(intersect(curEd4e(k,3),obj.ed4Sb))
                                    sbAxis = n4Sb(obj.ed4Sb==curEd4e(k,3),3);
                                    newSb = [newSb;
                                        curN4e(k,3),newNode(3),sbAxis;
                                        newNode(3),curN4e(k,4),sbAxis];
                                else
                                    newNb = [newNb;
                                        curN4e(k,3),newNode(3);
                                        newNode(3),curN4e(k,4)];
                                end
                            else
                                error('fem4m:MeshPolyH:cornerNotImplemented','corner not implemented');
                            end
                        end
                        % no corner
                        mid4e = 0.5*(obj.c4n(curN4e(:,1),:)+obj.c4n(curN4e(:,2),:));
                        curN4e = curN4e(~corner,:);
                        curEd4e = curEd4e(~corner,:);
                        parts = parts(~corner);
                        mid4e = mid4e(~corner,:);
                        newNode = zeros(nnz(parts),4);
                        for k = 1 : 4
                            newNode(:,k) = size(newC4n,1)+(1:numel(parts));
                            newC4n = [newC4n;(mid4e+obj.c4n(curN4e(:,k),:))/2];
                        end
                        newN4e = [newN4e;
                            mat2cell(newNode,ones(size(newNode,1),1))];
                        curNewN4e = [curN4e(:,[2 3]),newNode(:,[3 2]);
                            curN4e(:,[3 4]),newNode(:,[4 3]);
                            curN4e(:,[4 1]),newNode(:,[1 4])];
                        newN4e = [newN4e;
                            mat2cell(curNewN4e,ones(size(curNewN4e,1),1))];
                        [~,indx] = intersect(curEd4e(:,1),obj.ed4Db);
                        newDb = [newDb;
                            curN4e(indx,1),newNode(indx,1);
                            newNode(indx,1),newNode(indx,2);
                            newNode(indx,2),curN4e(indx,2);];
                        [~,indx] = intersect(curEd4e(:,1),obj.ed4Nb);
                        newNb = [newNb;
                            curN4e(indx,1),newNode(indx,1);
                            newNode(indx,1),newNode(indx,2);
                            newNode(indx,2),curN4e(indx,2);];
                        [~,indx,indSb] = intersect(curEd4e(:,1),obj.ed4Sb);
                        newSb = [newSb;
                            curN4e(indx,1),newNode(indx,1),n4Sb(indSb,3);
                            newNode(indx,1),newNode(indx,2),n4Sb(indSb,3);
                            newNode(indx,2),curN4e(indx,2),n4Sb(indSb,3);];
                    elseif j == 5
                        mid4e = 0.5*(obj.c4n(curN4e(:,1),:)+obj.c4n(curN4e(:,2),:));
                        newNode = zeros(numel(parts),5);
                        for k = 1 : 5
                            newNode(:,k) = size(newC4n,1)+(1:numel(parts));
                            newC4n = [newC4n;(mid4e+obj.c4n(curN4e(:,k),:))/2];
                        end
                        newN4e = [newN4e;
                            mat2cell(newNode,ones(size(newNode,1),1))];
                        curNewN4e = [curN4e(:,[2 3]),newNode(:,[3 2]);
                            curN4e(:,[3 4]),newNode(:,[4 3]);
                            curN4e(:,[4 5]),newNode(:,[5 4]);
                            curN4e(:,[5 1]),newNode(:,[1 5])];
                        newN4e = [newN4e;
                            mat2cell(curNewN4e,ones(size(curNewN4e,1),1))];
                        [~,indx] = intersect(curEd4e(:,1),obj.ed4Db);
                        newDb = [newDb;
                            curN4e(indx,1),newNode(indx,1);
                            newNode(indx,1),newNode(indx,2);
                            newNode(indx,2),curN4e(indx,2);];
                        for k = [2,5]
                            [~,indx] = intersect(curEd4e(:,k),obj.ed4Db);
                            newDb = [newDb; curN4e(indx,[k,mod(k,5)+1])];
                        end
                        [~,indx] = intersect(curEd4e(:,1),obj.ed4Nb);
                        newNb = [newNb;
                            curN4e(indx,1),newNode(indx,1);
                            newNode(indx,1),newNode(indx,2);
                            newNode(indx,2),curN4e(indx,2);];
                        for k = [2,5]
                            [~,indx] = intersect(curEd4e(:,k),obj.ed4Nb);
                            newNb = [newNb; curN4e(indx,[k,mod(k,5)+1])];
                        end
                        [~,indx,indSb] = intersect(curEd4e(:,1),obj.ed4Sb);
                        newSb = [newSb;
                            curN4e(indx,1),newNode(indx,1),n4Sb(indSb,3);
                            newNode(indx,1),newNode(indx,2),n4Sb(indSb,3);
                            newNode(indx,2),curN4e(indx,2),n4Sb(indSb,3);];
                        for k = [2,5]
                            [~,indx,indSb] = intersect(curEd4e(:,k),obj.ed4Sb);
                            newSb = [newSb; curN4e(indx,[k,mod(k,5)+1]),n4Sb(indSb,3)];
                        end

                    elseif j == 7

                        %Lshape corners
                        bdEdges = e4ed2(curEd4e)==0;
                        if size(curEd4e,1)==1
                            bdEdges = bdEdges';
                        end
                        for k = 1 : size(curN4e,1)
                            if bdEdges(k,1)==1 && bdEdges(k,7)==1
                                indx = 1;
                            else
                                indx = find(bdEdges(k,:),1,'last');
                            end
                            curN4e(k,:) = curN4e(k,mod(indx+(0:6)-1,7)+1);
                            curEd4e(k,:) = curEd4e(k,mod(indx+(0:6)-1,7)+1);
                        end
                        % first node is corner node
                        mid4e = obj.c4n(curN4e(:,1),:);
                        newNode = zeros(numel(parts),6);
                        for k = 1:6
                            newNode(:,k) = size(newC4n,1)+(1:numel(parts));
                            newC4n = [newC4n;(mid4e+obj.c4n(curN4e(:,k+1),:))/2];
                        end
                        newN4e = [newN4e;
                            mat2cell([curN4e(:,1),newNode],ones(size(newNode,1),1))];
                        curNewN4e = [
                            curN4e(:,[2 3]),newNode(:,[2 1]);
                            curN4e(:,[3 4]),newNode(:,[3 2]);
                            curN4e(:,[4 5]),newNode(:,[4 3]);
                            curN4e(:,[5 6]),newNode(:,[5 4]);
                            curN4e(:,[6 7]),newNode(:,[6 5])];
                        newN4e = [newN4e;
                            mat2cell(curNewN4e,ones(size(curNewN4e,1),1))];
                        [~,indx] = intersect(curEd4e(:,1),obj.ed4Db);
                        newDb = [newDb;
                            curN4e(indx,1),newNode(indx,1);
                            newNode(indx,1),curN4e(indx,2);];
                        [~,indx] = intersect(curEd4e(:,end),obj.ed4Db);
                        newDb = [newDb;
                            curN4e(indx,7),newNode(indx,6);
                            newNode(indx,6),curN4e(indx,1);];
                        [~,indx] = intersect(curEd4e(:,1),obj.ed4Nb);
                        newNb = [newNb;
                            curN4e(indx,1),newNode(indx,1);
                            newNode(indx,1),curN4e(indx,2);];
                        [~,indx] = intersect(curEd4e(:,end),obj.ed4Nb);
                        newNb = [newNb;
                            curN4e(indx,7),newNode(indx,6);
                            newNode(indx,6),curN4e(indx,1);];
                        [~,indx,indSb] = intersect(curEd4e(:,1),obj.ed4Sb);
                        newSb = [newSb;
                            curN4e(indx,1),newNode(indx,1),n4Sb(indSb,3);
                            newNode(indx,1),curN4e(indx,2),n4Sb(indSb,3);];
                        [~,indx,indSb] = intersect(curEd4e(:,end),obj.ed4Sb);
                        newSb = [newSb;
                            curN4e(indx,7),newNode(indx,6),n4Sb(indSb,3);
                            newNode(indx,6),curN4e(indx,1),n4Sb(indSb,3);];
                    elseif j == 8

                        %slit corners
                        bdEdges = e4ed2(curEd4e)==0;
                        if size(curEd4e,1)==1
                            bdEdges = bdEdges';
                        end
                        for k = 1 : size(curN4e,1)
                            if bdEdges(k,1)==1 && bdEdges(k,8)==1
                                indx = 1;
                            else
                                indx = find(bdEdges(k,:),1,'last');
                            end
                            curN4e(k,:) = curN4e(k,mod(indx+(0:7)-1,8)+1);
                            curEd4e(k,:) = curEd4e(k,mod(indx+(0:7)-1,8)+1);
                        end
                        % first node is corner node
                        mid4e = obj.c4n(curN4e(:,1),:);
                        newNode = zeros(numel(parts),7);
                        for k = 1:7
                            newNode(:,k) = size(newC4n,1)+(1:numel(parts));
                            newC4n = [newC4n;(mid4e+obj.c4n(curN4e(:,k+1),:))/2];
                        end
                        newN4e = [newN4e;
                            mat2cell([curN4e(:,1),newNode],ones(size(newNode,1),1))];
                        curNewN4e = [
                            curN4e(:,[2 3]),newNode(:,[2 1]);
                            curN4e(:,[3 4]),newNode(:,[3 2]);
                            curN4e(:,[4 5]),newNode(:,[4 3]);
                            curN4e(:,[5 6]),newNode(:,[5 4]);
                            curN4e(:,[6 7]),newNode(:,[6 5]);
                            curN4e(:,[7 8]),newNode(:,[7 6])];
                        newN4e = [newN4e;
                            mat2cell(curNewN4e,ones(size(curNewN4e,1),1))];
                        [~,indx] = intersect(curEd4e(:,1),obj.ed4Db);
                        newDb = [newDb;
                            curN4e(indx,1),newNode(indx,1);
                            newNode(indx,1),curN4e(indx,2);];
                        [~,indx] = intersect(curEd4e(:,end),obj.ed4Db);
                        newDb = [newDb;
                            curN4e(indx,8),newNode(indx,7);
                            newNode(indx,7),curN4e(indx,1);];
                        [~,indx] = intersect(curEd4e(:,1),obj.ed4Nb);
                        newNb = [newNb;
                            curN4e(indx,1),newNode(indx,1);
                            newNode(indx,1),curN4e(indx,2);];
                        [~,indx] = intersect(curEd4e(:,end),obj.ed4Nb);
                        newNb = [newNb;
                            curN4e(indx,8),newNode(indx,7);
                            newNode(indx,7),curN4e(indx,1);];
                        [~,indx,indSb] = intersect(curEd4e(:,1),obj.ed4Sb);
                        newSb = [newSb;
                            curN4e(indx,1),newNode(indx,1),n4Sb(indSb,3);
                            newNode(indx,1),curN4e(indx,2),n4Sb(indSb,3);];
                        [~,indx,indSb] = intersect(curEd4e(:,end),obj.ed4Sb);
                        newSb = [newSb;
                            curN4e(indx,8),newNode(indx,7),n4Sb(indSb,3);
                            newNode(indx,7),curN4e(indx,1),n4Sb(indSb,3);];
                    else
                        newN4e = [newN4e;
                            mat2cell(curN4e,ones(size(curN4e,1),1))];
                        newDb = [newDb;
                            obj.n4Db(obj.nrN4e(obj.e4ed(obj.ed4Db,1))==j & r4e(obj.e4ed(obj.ed4Db,1)),:)];
                        newNb = [newNb;
                            obj.n4Nb(obj.nrN4e(obj.e4ed(obj.ed4Nb,1))==j & r4e(obj.e4ed(obj.ed4Nb,1)),:)];
                        newSb = [newSb;
                            n4Sb(obj.nrN4e(obj.e4ed(obj.ed4Sb,1))==j & r4e(obj.e4ed(obj.ed4Sb,1)),:)];
                    end
                end
            end

            obj = MeshPolyH(newC4n,newN4e,newDb,newNb,newSb);
            obj = obj.joinElements;

        end

        function r4e = closure(obj,r4e)
            %CLOSURE - Enforce the closure condition on marked elements
            %   R4E = CLOSURE(OBJ,R4E) moves the mark of every marked
            %   quadrilateral whose longest edge is an inner edge and
            %   which has at most one boundary edge to its neighbor
            %   across that edge (the quadrilateral itself is unmarked)
            %   and repeats this until no such quadrilateral is marked.
            %   R4E is a logical or 0/1 vector with one entry per
            %   element; refine then creates no hanging nodes.

            e4ed2 = obj.e4ed(:,2);
            length4ed = obj.length4ed;
            indx = 1;
            while nnz(indx)>0
                parts = find(obj.nrN4e==4 & r4e);
                if ~isempty(parts)
                    curEd4e = cell2mat(obj.ed4e(parts));
                    % rotate ed4e: the longest edge comes first
                    if size(curEd4e,1)==1
                        length = length4ed(curEd4e)';
                    else
                        length = length4ed(curEd4e);
                    end
                    [~,indy] = max(length,[],2);
                    for k = 2 : 4
                        curEd4e(indy==k,:) = curEd4e(indy==k,mod((k:k+4-1)-1,4)+1);
                    end
                    if size(curEd4e,1)==1
                        indx = e4ed2(curEd4e(:,1))~=0 & (1>=sum(e4ed2(curEd4e)'==0,2));
                    else
                        indx = e4ed2(curEd4e(:,1))~=0 & (1>=sum(e4ed2(curEd4e)==0,2));
                    end
                    if nnz(indx)>0
                        r4e(obj.e4ed(curEd4e(indx,1),:))=1;
                        r4e(parts(indx)) = 0;
                    end
                else
                    indx = [];
                end
            end
        end

        function obj = joinElements(obj)
            %joinElements - Join pairs of quadrilaterals into polygons
            %   OBJ = joinElements(OBJ) merges every two quadrilaterals of
            %   equal diameter that share an inner edge, the first edge
            %   (in n4e) of one of them, into one polygon: a hexagon when
            %   one of the two has no boundary edge, otherwise a pentagon
            %   without the boundary node of the shared edge, whose two
            %   boundary edges become one. Nodes no longer used are
            %   removed; without such pairs OBJ is returned unchanged.

            newC4n = obj.c4n;
            newN4e = obj.n4e;
            newDb = obj.n4Db;
            newNb = obj.n4Nb;
            newSb = obj.n4Sb;
            if isempty(newSb)
                newSb = zeros(0,3);
            end
            diam4e = obj.diam4e;

            curEd4e = cell2mat(obj.ed4e(obj.nrN4e==4));
            innerParts = find(obj.e4ed(curEd4e(:,1),2)~=0);
            joinParts = obj.nrN4e(obj.e4ed(curEd4e(innerParts,1),1))==4 ...
                & obj.nrN4e(obj.e4ed(curEd4e(innerParts,1),2))==4 ...
                & abs(diam4e(obj.e4ed(curEd4e(innerParts,1),1))...
                -diam4e(obj.e4ed(curEd4e(innerParts,1),2)))...
                ./abs(diam4e(obj.e4ed(curEd4e(innerParts,1),1)))<100*eps;
            joinEdges = unique(curEd4e(innerParts(joinParts),1));

            if ~isempty(joinEdges)
                curN4e1 = cell2mat(obj.n4e(obj.e4ed(joinEdges,1)));
                curN4e2 = cell2mat(obj.n4e(obj.e4ed(joinEdges,2)));
                curEd4e1 = cell2mat(obj.ed4e(obj.e4ed(joinEdges,1)));
                curEd4e2 = cell2mat(obj.ed4e(obj.e4ed(joinEdges,2)));
                e4ed2 = obj.e4ed(:,2);
                if size(curEd4e1,1)==1
                    indx = all(e4ed2(curEd4e1)'~=0,2) | all(e4ed2(curEd4e2)'~=0,2);
                else
                    indx = all(e4ed2(curEd4e1)~=0,2) | all(e4ed2(curEd4e2)~=0,2);
                end
                curN4e = [curN4e1(indx,[2 3 4 1]),curN4e2(indx,[3 4])];
                newN4e(obj.e4ed(joinEdges(indx),1)) = mat2cell(curN4e,ones(size(curN4e,1),1));
                indy = e4ed2(curEd4e1(:,2))==0;
                curN4e = [curN4e2(~indx & indy,4),curN4e1(~indx & indy,[3 4 1]),curN4e2(~indx & indy,3)];
                newN4e(obj.e4ed(joinEdges(~indx & indy),1)) = mat2cell(curN4e,ones(size(curN4e,1),1));

                boundaryNode = curN4e1(~indx & indy,2);
                boundary = [curN4e2(~indx & indy,4),curN4e1(~indx & indy,3)];
                for k = 1 : numel(boundaryNode)
                    curNode = boundaryNode(k);
                    indz = find(newDb(:,1)==curNode | newDb(:,2)==curNode);
                    if ~isempty(indz)
                        newDb(indz(1),:) = boundary(k,:);
                        newDb(indz(2),:)= [];
                    end
                    indz = find(newNb(:,1)==curNode | newNb(:,2)==curNode);
                    if ~isempty(indz)
                        newNb(indz(1),:) = boundary(k,:);
                        newNb(indz(2),:)= [];
                    end
                    indz = find(newSb(:,1)==curNode | newSb(:,2)==curNode);
                    if ~isempty(indz)
                        newSb(indz(1),1:2) = boundary(k,:);
                        newSb(indz(2),:)= [];
                    end
                end
                curN4e = [curN4e1(~indx & ~indy,4),curN4e2(~indx & ~indy,[3 4 1]),curN4e1(~indx & ~indy,3)];
                newN4e(obj.e4ed(joinEdges(~indx & ~indy),1)) = mat2cell(curN4e,ones(size(curN4e,1),1));

                boundaryNode = curN4e2(~indx & ~indy,2);
                boundary = [curN4e1(~indx & ~indy,4),curN4e2(~indx & ~indy,3)];
                for k = 1 : numel(boundaryNode)
                    curNode = boundaryNode(k);
                    indz = find(newDb(:,1)==curNode | newDb(:,2)==curNode);
                    if ~isempty(indz)
                        newDb(indz(1),:) = boundary(k,:);
                        newDb(indz(2),:)= [];
                    end
                    indz = find(newNb(:,1)==curNode | newNb(:,2)==curNode);
                    if ~isempty(indz)

                        newNb(indz(1),:) = boundary(k,:);
                        newNb(indz(2),:)= [];
                    end
                    indz = find(newSb(:,1)==curNode | newSb(:,2)==curNode);
                    if ~isempty(indz)
                        newSb(indz(1),1:2) = boundary(k,:);
                        newSb(indz(2),:)= [];
                    end
                end

                newN4e(obj.e4ed(joinEdges,2)) = [];
                [newC4n,newN4e,newDb,newNb,newSb] = MeshPoly.removeGhostNodes(newC4n,newN4e,newDb,newNb,newSb);
                obj = MeshPolyH(newC4n,newN4e,newDb,newNb,newSb);
            end
        end
    end

end