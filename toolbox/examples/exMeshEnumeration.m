%[text] # Enumeration of nodes, edges and elements of a mesh
%[text] Introduces the arrays in which fem4m stores a triangulation: the coordinates of the nodes, the nodes of the elements, the edges with their neighbouring elements, the incidences of the nodes and the boundary edges, together with their orientation conventions. A local refinement towards a corner shows that the refined meshes stay regular and satisfy the Euler relation.
%[text] **Builds on:** none, this is the first example. **Next:** `exElementBasisFunctions`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### Regular triangulations
%[text] Let $\\Omega\\subset\\mathbb{R}^2$ be a bounded polygonal domain. A regular triangulation $\\mathcal{T}$ of $\\Omega$ is a finite set of closed triangles $T$ with $\\overline\\Omega = \\bigcup\_{T\\in\\mathcal{T}} T$ such that two different triangles intersect in the empty set, in a common vertex or in a common edge. In particular there are no hanging nodes: no vertex of a triangle lies in the interior of an edge of another triangle. The vertices of the triangles are the nodes $\\mathcal{N}$ and their edges the edges $\\mathcal{E}$ of $\\mathcal{T}$. The interior edges form $\\mathcal{E}(\\Omega)$, the edges on the boundary $\\mathcal{E}(\\partial\\Omega)$, which splits into the Dirichlet edges $\\mathcal{E}(\\Gamma\_D)$ and the Neumann edges $\\mathcal{E}(\\Gamma\_N)$.
%[text] ### Data arrays
%[text] A `Mesh` numbers the nodes, edges and elements from 1 to `nrNodes`, `nrEdges` and `nrElems` and stores the mesh in arrays whose names read "x for y": row $y$ of `x4y` holds the x of the object $y$, e.g. `c4n(z,:)` the coordinates of the node $z$ and `n4e(T,:)` the nodes of the element $T$. The incidences `ed4n` and `e4n` are sparse matrices.
%[text:table]
%[text] | Array | Size | Content |
%[text] | --- | --- | --- |
%[text] | `c4n` | nrNodes-by-2 | coordinates of the nodes |
%[text] | `n4e` | nrElems-by-3 | nodes of the elements, counterclockwise |
%[text] | `n4ed` | nrEdges-by-2 | nodes of the edges |
%[text] | `ed4e` | nrElems-by-3 | edges of the elements: `ed4e(T,:)` are the edges `n4e(T,[1 2])`, `n4e(T,[2 3])` and `n4e(T,[3 1])` |
%[text] | `e4ed` | nrEdges-by-2 | the elements $T\_+$ and $T\_-$ of an edge, $T\_- = 0$ on the boundary |
%[text] | `ed4n` | nrNodes-by-nrNodes, sparse | `ed4n(i,j) = ed4n(j,i)` is the edge between the nodes $i$ and $j$, zero if there is none |
%[text] | `e4n` | nrNodes-by-nrElems, sparse | `e4n(z,T) = T` if $z$ is a node of $T$: row $z$ holds the elements of the patch of $z$ |
%[text] | `n4Db`, `n4Nb` | nrDbEdges-by-2, nrNbEdges-by-2 | nodes of the Dirichlet and of the Neumann edges |
%[text] | `ed4Db`, `ed4Nb` | column vectors | ids of the Dirichlet and of the Neumann edges, in the order of `n4Db` and `n4Nb` |
%[text:table]
%[text] `c4n`, `n4e` and the boundary arrays `n4Db` and `n4Nb` define the mesh, `Mesh` computes the other arrays from them, following the usual data structures.
%[text] ### Orientation conventions
%[text] - The nodes of every element are counterclockwise, so `area4e` is positive. The first edge `n4e(T,[1 2])` is the refinement edge of $T$, and the third node, opposite to it, is the newest vertex.
%[text] - An interior edge $E$ is stored with increasing node ids, `n4ed(E,1) < n4ed(E,2)`. The element $T\_+$ = `e4ed(E,1)` traverses $E$ in this direction, from `n4ed(E,1)` to `n4ed(E,2)`, and $T\_-$ = `e4ed(E,2)` in the opposite direction.
%[text] - The unit normal $n\_E$ of an edge (`normal4ed`) is its unit tangent from `n4ed(E,1)` to `n4ed(E,2)` rotated clockwise by a quarter turn. It is the outer normal of $T\_+$ and points from $T\_+$ into $T\_-$.
%[text] - A boundary edge keeps the orientation of its element, in `n4ed` as well as in `n4Db` and `n4Nb`: the domain lies to the left of it, $T\_+$ is its only element, and $n\_E$ is the outer unit normal of $\\Omega$. \
%[text] ### The Euler relation
%[text] Every triangle has three edges, an interior edge belongs to two triangles and a boundary edge to one, so $3|\\mathcal{T}| = 2|\\mathcal{E}| - |\\mathcal{E}(\\partial\\Omega)|$. For a simply connected domain, Euler's polyhedron formula for the planar graph of the nodes and edges gives
%[text]{"align":"center"} $ |\\mathcal{N}| - |\\mathcal{E}| + |\\mathcal{T}| = 1. $
%[text] For a domain with $m$ holes the right-hand side is $1-m$. On fine meshes, where the boundary edges are few, the two relations give $|\\mathcal{E}|\\approx 3|\\mathcal{N}|$ and $|\\mathcal{T}|\\approx 2|\\mathcal{N}|$. This determines the size of finite element spaces; the space of continuous piecewise quadratic functions, for instance, has one degree of freedom per node and per edge, about $4|\\mathcal{N}|$ in total.
%[text] ### Refinement
%[text] `mesh.refine(marked4ed)` refines the edges marked by a 0/1 vector with one entry per edge; a vector with one entry per element marks all edges of the selected elements, and `refine('uniform')` marks all edges. To keep the mesh regular, the closure marks the refinement edge of every element with a marked edge, repeatedly. The new nodes are the midpoints of the marked edges; they are appended to `c4n`, so the old nodes keep their ids, and `parent4e` records the parent of every new element.
%[text] - Red-green-blue refinement (`'RGB'`, the default) divides an element with three marked edges into four similar triangles (red), an element with only its refinement edge marked into two (green), and an element with two marked edges into three (blue).
%[text] - Newest vertex bisection (`'NVB'`) bisects an element at the midpoint of its refinement edge. The midpoint becomes the newest vertex of both children, so their refinement edges are the edges opposite to it; an element with two or three marked edges is bisected once more in one or both children. \
%[text] Both refinements keep the angles of the triangles bounded away from zero, so the refined meshes stay shape regular.
%%
%[text] ## Parameters
reftype = 'crisscross';  % coarse mesh of the unit square: both diagonals, four triangles
% reftype = 'criss';     % one diagonal, two triangles
refinement = 'RGB';      % red-green-blue refinement, the default of refine
% refinement = 'NVB';    % newest vertex bisection
nrLevels = 4;            % number of local refinements towards the corner (0,0)
%%
%[text] ## Initial mesh
%[text] `Mesh` creates the coarse mesh of the unit square from the options of a `Geometry`; all boundary edges are Dirichlet edges. Without a semicolon, a mesh with at most 100 nodes displays the arrays that define it: `c4n`, `n4e` and `n4Db`.
mesh = Mesh(struct('name','Square','reftype',reftype))
%[text] `plotEnumeration` writes the ids of the nodes (green), of the elements at their centroids (yellow) and of the edges at their midpoints (purple) into the plot of the mesh.
plot(mesh);
hold on
mesh.plotEnumeration;
hold off
axis equal
%%
%[text] ## Local refinement
%[text] Every step marks the elements of the patch of node 1, the corner (0,0), which are the nonzero entries of row 1 of `e4n`, and refines them; the closure refines the neighbouring elements as far as needed to keep the mesh regular. The loop records the numbers of nodes, edges, elements and boundary edges, the Dirichlet and Neumann edges, of every mesh.
refinedMesh = mesh;
nrNodes4lvl = zeros(1,nrLevels+1);
nrEdges4lvl = zeros(1,nrLevels+1);
nrElems4lvl = zeros(1,nrLevels+1);
nrBdEdges4lvl = zeros(1,nrLevels+1);
for lvl = 1:nrLevels+1
    if lvl > 1
        % mark
        marked4e = full(refinedMesh.e4n(1,:)) > 0;
        % refine
        refinedMesh = refinedMesh.refine(refinement,marked4e);
    end
    % count
    nrNodes4lvl(lvl) = refinedMesh.nrNodes;
    nrEdges4lvl(lvl) = refinedMesh.nrEdges;
    nrElems4lvl(lvl) = refinedMesh.nrElems;
    nrBdEdges4lvl(lvl) = numel(refinedMesh.ed4Db) + numel(refinedMesh.ed4Nb);
end
%%
%[text] ## Results
%[text] ### Elements
%[text] Row $T$ lists the nodes of the element counterclockwise and its edges from the first to the second, from the second to the third and from the third to the first node.
table((1:mesh.nrElems)',mesh.n4e,mesh.ed4e,VariableNames=["T" "n4e" "ed4e"])
%%
%[text] ### Edges
%[text] The boundary edges have only the element $T\_+$; the interior edges have increasing node ids and two elements.
table((1:mesh.nrEdges)',mesh.n4ed,mesh.e4ed,VariableNames=["E" "n4ed" "e4ed"])
%%
%[text] ### Edges and elements at a node
%[text] The rows of the sparse arrays `ed4n` and `e4n` give the edges at node 1 and the elements of its patch.
edges4node1 = nonzeros(mesh.ed4n(1,:))'
elems4node1 = nonzeros(mesh.e4n(1,:))'
%%
%[text] ### Boundary edges
%[text] `ed4Db` gives the ids of the Dirichlet edges; the rows of `n4ed` for these ids are the rows of `n4Db`.
table(mesh.ed4Db,mesh.n4ed(mesh.ed4Db,:),VariableNames=["E" "n4ed"])
%%
%[text] ### Edge normals
clf
plot(mesh);
hold on
mid4ed = mesh.point4ed(1/2);
normal4ed = mesh.normal4ed;
quiver(mid4ed(:,1),mid4ed(:,2),normal4ed(:,1),normal4ed(:,2),0.3,'LineWidth',1.5);
hold off
axis equal
%[text] The normals of the interior edges point from $T\_+$ into $T\_-$, those of the boundary edges out of the domain.
%%
%[text] ### Refined mesh
clf
plot(refinedMesh);
axis equal
%[text] The refinement concentrates at the corner (0,0). The closure has refined further elements, so that the mesh has no hanging nodes.
%%
%[text] ### Euler relation
table((0:nrLevels)',nrNodes4lvl',nrEdges4lvl',nrElems4lvl',nrBdEdges4lvl', ...
    (nrNodes4lvl-nrEdges4lvl+nrElems4lvl)',(2*nrEdges4lvl-nrBdEdges4lvl-3*nrElems4lvl)', ...
    VariableNames=["level" "nodes" "edges" "elements" "boundary edges" "N-E+T" "2E-EB-3T"])
%[text] On every mesh $|\\mathcal{N}| - |\\mathcal{E}| + |\\mathcal{T}| = 1$ and $2|\\mathcal{E}| - |\\mathcal{E}(\\partial\\Omega)| - 3|\\mathcal{T}| = 0$, as for every regular triangulation of the square. A hanging node would violate the second relation: the edges at it would belong to one element only without being boundary edges. The closure guarantees that no hanging nodes arise.
%%
%[text] ### Dirichlet and Neumann boundary
%[text] `Mesh(NAME)` loads a mesh from the files `NAME_c4n.dat`, `NAME_n4e.dat`, `NAME_n4Db.dat` and `NAME_n4Nb.dat` on the path, here Cook's membrane from the folder `geometries/cooks`, refined once uniformly by newest vertex bisection. `plotBd` draws the Dirichlet edges as solid red and the Neumann edges as dashed blue lines into the current axes and leaves them on hold.
clf
meshCooks = Mesh('Cooks').refine('uniform','NVB');
plot(meshCooks);
meshCooks.plotBd;
hold off
axis equal
%%
%[text] ## Possibilities for extensions
%[text] - Set `refinement = 'NVB'` and compare the plot of the refined mesh with that of red-green-blue refinement; for this marking the numbers in the Euler table are the same.
%[text] - Mark edges instead of elements, e.g. the boundary edges with `refinedMesh.refine(refinedMesh.e4ed(:,2) == 0)`, and watch the closure refine the neighbouring elements.
%[text] - `mesh.patch4n` collects the elements of all node patches from `e4n`, `mesh.area4e` and `mesh.length4ed` compute the areas of the elements and the lengths of the edges.
%[text] - With `refinement = 'NVB'`, `refinedMesh.coarse(elems)` undoes the newest vertex bisections at the newest vertices of the elements `elems`; it does not undo red-green-blue refinements.
%[text] - `Mesh3d` stores tetrahedral meshes in the same way, with the faces in addition (`n4f`, `f4e`, `e4f`). \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
