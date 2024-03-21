function P = prolongationMatrix(fe,feOld)
%prolongationMatrix - Prolongation between finite element spaces on refined meshes
%   P = prolongationMatrix(FE,FEOLD) returns the sparse matrix that maps
%   the coefficients of a function of FEOLD, defined on the mesh before
%   the last refine call, to the coefficients of the same function in FE.
%   P has one row per degree of freedom of FE and one column per degree
%   of freedom of FEOLD. FE is of the class of FEOLD or of a subclass of
%   it with the same kind of values (an enrichment such as HHJB of
%   HHJk). P is exact for nested spaces and otherwise a local L2
%   projection, see below. Spaces whose degrees of freedom are not
%   numbered per element, such as the face space DGPkFace, and spaces
%   whose basis is not defined on the reference element (isRefElemBased
%   false, e.g. PDGk) have no prolongation.
%
%   The routine uses only the finite element independent refinement data
%   of FE.mesh: parent4e, the ancestor of every element in the old mesh,
%   and the composite P1 prolongation, whose rows hold the barycentric
%   coordinates of the new nodes. Their entries give for every element
%   the barycentric map B of its vertices into the ancestor. Elements are
%   grouped by B, so the local matrices are computed only once per group.
%
%   Elements that are affine equivalent on the reference element (all 3d
%   elements and most 2d elements) use reference matrices,
%       L(e) = A(o_F(e)) * R(B(e)) / A(o_C(e)),
%   where R is the local L2 projection from the reference element onto
%   the reference child with map B and A(o) maps the reference basis to
%   the basis of an element of orientation class o (rank pattern of the
%   vertex ids and, in 2d, the boundary edges), both computed from the
%   element class itself on reference simplices, where
%   feval(class(FE),FE.degree,MESH) builds the same space. Other
%   elements (isAffineEquivalent false, or element-wise degrees) use the
%   local L2 projection on the physical elements (with zero rows for
%   the unused degrees of freedom nonDoFNodes). For spaces that are
%   not nested and whose basis is mapped by a Piola-type transform
%   (bubble spaces such as CurlBk, and HHJB), the reference matrices give
%   the L2 projection in reference coordinates.
%   The local matrices are averaged over the elements that contain a
%   fine degree of freedom.
%
%   prolongationMatrix raises an error if FE.mesh is not the result of
%   a refine call (no parent4e or prolongation), if FEOLD does not live
%   on the mesh before that call, if the classes do not match as stated
%   above, if feval(class(FE),FE.degree,MESH) does not rebuild the space
%   (reference matrices only), or if the basis is not a polynomial of
%   degree at most the degree of the space plus 4.
%
%   See also FE, FEVector, Mesh, Mesh3d

% Copyright (c) 2026 by the fem4m authors

meshF = fe.mesh;
meshC = feOld.mesh;
if ~isprop(meshF,'parent4e') || isempty(meshF.parent4e) ...
        || ~isprop(meshF,'prolongation') || isempty(meshF.prolongation)
    error('fem4m:FE:noParentMesh', ...
        'The mesh of %s stores no parent elements: it is not the result of a refine call.',class(fe));
end
par = double(meshF.parent4e(:));
P1 = meshF.prolongation;
if numel(par) ~= meshF.nrElems || max(par) > meshC.nrElems ...
        || ~isequal(size(P1),[meshF.nrNodes,meshC.nrNodes])
    error('fem4m:FE:notParentMesh', ...
        'The old space must live on the mesh before the last refine call of the new mesh.');
end
if ~isa(fe,class(feOld)) || (~strcmp(class(fe),class(feOld)) ...
        && ~isequal(valueSize(fe),valueSize(feOld)))
    error('fem4m:FE:prolongationClassMismatch', ...
        ['The new space must be of the class of the old one or of a subclass with ', ...
        'the same kind of values, not %s and %s.'],class(fe),class(feOld));
end
if ~fe.isRefElemBased
    error('fem4m:FE:prolongationUnsupported', ...
        'No prolongation for %s: its basis is not defined on the reference element.',class(fe));
end
if size(fe.dof4e,1) ~= meshF.nrElems || size(feOld.dof4e,1) ~= meshC.nrElems
    error('fem4m:FE:prolongationUnsupported', ...
        'No prolongation for %s: its degrees of freedom are not numbered per element (dof4e).',class(fe));
end

% barycentric map of every element into its ancestor,
% B(e,i,j) = j-th barycentric coordinate of vertex i of element e
n = size(meshF.n4e,2);
nE = meshF.nrElems;
n4eF = double(meshF.n4e);
n4eC = double(meshC.n4e(par,:));
B = zeros(nE,n,n);
for i = 1 : n
    for j = 1 : n
        B(:,i,j) = full(P1(sub2ind(size(P1),n4eF(:,i),n4eC(:,j))));
    end
end
if any(abs(sum(B,3)-1) > 1e-10,'all') || any(B(:) < -1e-10)
    error('fem4m:FE:notParentMesh', ...
        'The old space must live on the mesh before the last refine call of the new mesh.');
end
[mapsB,~,group] = unique(reshape(B,nE,n*n),'rows');

if fe.isAffineEquivalent && feOld.isAffineEquivalent && uniformDegree(fe) && uniformDegree(feOld)
    L = factoredLocal(fe,feOld,mapsB,group,par);
else
    L = physicalLocal(fe,feOld,mapsB,group,par);
end

% assembly with averaging over the elements of every fine dof
nRows = max([fe.nrDoF;double(fe.dof4e(:))]);
nCols = max([feOld.nrDoF;double(feOld.dof4e(:))]);
S = permute(L,[3 2 1]);
[I,J] = localDoFtoGlobalDoF(feOld.dof4e(par,:),fe.dof4e);
k = find(S);
count = accumarray(double(fe.dof4e(:)),1,[nRows,1]);
P = spdiags(1./max(count,1),0,nRows,nRows)*sparse(I(k),J(k),S(k),nRows,nCols);
end

function val = uniformDegree(fe)
%uniformDegree - True unless the element has element-wise degrees
%   VAL = uniformDegree(FE) is false if FE has a property degree4e with
%   an entry other than FE.degree, and true otherwise.
val = ~isprop(fe,'degree4e') || all(fe.degree4e(:) == fe.degree);
end

function L = factoredLocal(fe,feOld,mapsB,group,par)
%factoredLocal - Local matrices from reference matrices
%   L = factoredLocal(FE,FEOLD,MAPSB,GROUP,PAR) returns the local
%   matrices L(:,:,e) = A(o_F(e))*R(B(e))/A(o_C(e)) of all elements e of
%   FE.mesh (FE.nrBasisFunc-by-FEOLD.nrBasisFunc-by-nrElems), see above.
%   MAPSB holds the distinct barycentric maps B (one flattened n-by-n
%   map per row), GROUP the row of MAPSB of every element and PAR its
%   ancestor in FEOLD.mesh.
n = size(fe.mesh.n4e,2);
d = n - 1;
X = [zeros(1,d);eye(d)];
[oF,bF] = orientationClass(fe.mesh);
[oC,bC] = orientationClass(feOld.mesh);
oC = oC(par);
bC = bC(par);
AF = orientationMatrices(fe);
AC = orientationMatrices(feOld);
% children of all maps, placed disjointly in the reference simplex
K = size(mapsB,1);
c4n = zeros(n*K,d);
for k = 1 : K
    c4n((k-1)*n+(1:n),:) = reshape(mapsB(k,:),n,n)*X;
end
R = refProlongation(sameSpace(fe,simplexMesh(c4n,reshape(1:n*K,n,K)')), ...
    sameSpace(feOld,simplexMesh(X,1:n)));
[keys,~,ik] = unique([group,oF,bF,oC,bC],'rows');
Lk = zeros(fe.nrBasisFunc,feOld.nrBasisFunc,size(keys,1));
for k = 1 : size(keys,1)
    Lk(:,:,k) = AF(:,:,keys(k,2),keys(k,3))*R(:,:,keys(k,1))/AC(:,:,keys(k,4),keys(k,5));
end
Lk(abs(Lk) < 1e-12) = 0;
L = Lk(:,:,ik);
end

function [o,b] = orientationClass(mesh)
%orientationClass - Rank pattern of the vertex ids and boundary edges
%   [O,B] = orientationClass(MESH) returns for every element of MESH the
%   index O of the rank pattern of its vertex ids in
%   sortrows(perms(1:n)) and the boundary class B: B-1 encodes in 2d
%   the boundary edges of the element (bit j for edge j), B = 1 in 3d.
[nE,n] = size(mesh.n4e);
[~,order] = sort(double(mesh.n4e),2);
ranks = zeros(nE,n);
ranks(sub2ind([nE,n],repmat((1:nE)',1,n),order)) = repmat(1:n,nE,1);
[~,o] = ismember(ranks,sortrows(perms(1:n)),'rows');
if n == 3
    boundary = reshape(mesh.e4ed(mesh.ed4e,2) == 0,nE,3);
    b = boundary*[1;2;4] + 1;
else
    b = ones(nE,1);
end
end

function A = orientationMatrices(fe)
%orientationMatrices - Reference basis to the basis of every orientation
%   A = orientationMatrices(FE) returns the matrices A(:,:,o,b) for all
%   orientation classes (o,b) of the element class of FE. A(:,:,o,b)
%   maps the reference basis (ids in increasing order, all edges on the
%   boundary) to the basis of an element of the same geometry with
%   orientation class (o,b), see orientationClass. In 2d the interior
%   edges get a neighbour triangle.
n = size(fe.mesh.n4e,2);
d = n - 1;
X = [zeros(1,d);eye(d)];
Q = sortrows(perms(1:n));
nb = 1 + 7*(d == 2);
c4n = zeros(0,d);
n4e = zeros(0,n);
centre = zeros(size(Q,1),nb);
for o = 1 : size(Q,1)
    for b = 1 : nb
        offset = size(c4n,1);
        ids = offset + Q(o,:);
        c4n(ids,:) = X;
        centre(o,b) = size(n4e,1) + 1;
        n4e = [n4e;ids]; %#ok<AGROW>
        for j = 1 : 3*(d == 2)
            if bitget(b-1,j)
                continue % edge j on the boundary
            end
            j2 = mod(j,3) + 1;
            j3 = mod(j+1,3) + 1;
            c4n = [c4n;X(j,:)+X(j2,:)-X(j3,:)]; %#ok<AGROW>
            n4e = [n4e;ids([j2 j]),size(c4n,1)]; %#ok<AGROW>
        end
    end
end
Ap = refProlongation(sameSpace(fe,simplexMesh(c4n,n4e)),sameSpace(fe,simplexMesh(X,1:n)));
A = reshape(Ap(:,:,centre(:)),size(Ap,1),size(Ap,2),size(Q,1),nb);
end

function sz = valueSize(fe)
%valueSize - Sizes of the components of the basis values of FE
%   SZ = valueSize(FE) returns the sizes of the dimensions 3 and 4 of
%   the basis values of FE at the barycentre of the first element, e.g.
%   [1 1] for a scalar and [2 1] for a 2d vector-valued space.
x = ones(1,size(fe.mesh.n4e,2)-1)/size(fe.mesh.n4e,2);
pts = fe.mesh.point4e(x);
sz = size(fe.getBasis(pts(1,:),x,1),[3 4]);
end

function feM = sameSpace(fe,mesh)
%sameSpace - The space of FE on another mesh
%   FEM = sameSpace(FE,MESH) returns feval(class(FE),FE.degree,MESH) and
%   raises an error if FEM has another number of basis functions per
%   element than FE.
feM = feval(class(fe),fe.degree,mesh);
if feM.nrBasisFunc ~= fe.nrBasisFunc
    error('fem4m:FE:prolongationUnsupported', ...
        'No prolongation for %s: %s(DEGREE,MESH) does not rebuild the same space.', ...
        class(fe),class(fe));
end
end

function mesh = simplexMesh(c4n,n4e)
%simplexMesh - Mesh of given simplices with complete boundary lists
%   MESH = simplexMesh(C4N,N4E) returns the Mesh (2d) or Mesh3d (3d) with
%   the node coordinates C4N and the simplices N4E, whose boundary edges
%   or faces are all Dirichlet boundary.
if size(c4n,2) == 2
    edges = [n4e(:,[1 2]);n4e(:,[2 3]);n4e(:,[3 1])];
    [~,~,j] = unique(sort(edges,2),'rows');
    once = accumarray(j,1) == 1;
    mesh = Mesh(c4n,n4e,edges(once(j),:),[]);
else
    mesh = Mesh3d(c4n,n4e,Mesh3d.getBoundary(n4e),[],[],'t4e',zeros(size(n4e,1),1));
end
end

function R = refProlongation(newFE,oldFE)
%refProlongation - Local L2 projections onto the elements of a mesh
%   R = refProlongation(NEWFE,OLDFE) projects the basis of OLDFE, which
%   lives on the reference simplex alone, onto the basis of NEWFE on
%   every element of its mesh (in reference coordinates of OLDFE).
%   R(:,:,e) is the NEWFE.nrBasisFunc-by-OLDFE.nrBasisFunc matrix of
%   element e. The basis of OLDFE is evaluated at arbitrary points by
%   exact polynomial interpolation.
mesh = newFE.mesh;
nE = mesh.nrElems;
evalOld = oldEvaluator(oldFE,size(mesh.n4e,2)-1);
flat = @(v)reshape(v,size(v,1),size(v,2),[]);
mix = @(a,b)permute(pagemtimes(permute(flat(a),[2 3 1]),permute(flat(b),[3 2 1])),[3 1 2]);
degree = 2*max(newFE.degree,oldFE.degree) + 2;
M = integrate(@(x,xr,i)mix(newFE.getBasis(x,xr,i),newFE.getBasis(x,xr,i)), ...
    mesh.n4e,1:nE,mesh,degree);
Bm = integrate(@(x,xr,i)mix(newFE.getBasis(x,xr,i),evalOld(x)),mesh.n4e,1:nE,mesh,degree);
R = pagemldivide(permute(M,[2 3 1]),permute(Bm,[2 3 1]));
R(abs(R) < 1e-12) = 0;
end

function f = oldEvaluator(oldFE,d)
%oldEvaluator - Evaluate the basis of a one-element space anywhere
%   F = oldEvaluator(OLDFE,D) returns a function handle F(Y) that
%   evaluates the basis of OLDFE, a space on the reference simplex of
%   dimension D alone, at the points Y (one row per point): one row per
%   point, one column per basis function, the components in the third
%   dimension. The basis is interpolated in the Dubiner basis of the
%   lowest degree P = OLDFE.degree+1,...,OLDFE.degree+4 that reproduces
%   it at a test point; otherwise an error is raised.
nb = oldFE.nrBasisFunc;
for p = oldFE.degree + (1:4)
    E = exponents(d,p);
    Z = E/p;
    V = zeros(size(Z,1),0);
    for l = 1 : size(Z,1)
        b = oldFE.getBasis(Z(l,:),Z(l,:),1);
        V(l,1:numel(b)) = b(:)';
    end
    C = dubiner(Z,p)\V;
    f = @(Y)reshape(dubiner(Y,p)*C,size(Y,1),nb,[]);
    y = [0.21 0.17 0.13];
    y = y(1:d);
    b = oldFE.getBasis(y,y,1);
    if norm(reshape(f(y),1,[]) - b(:)',inf) < 1e-9*max(1,norm(b(:),inf))
        return
    end
end
error('fem4m:FE:prolongationUnsupported', ...
    'No prolongation for %s: its basis is not a polynomial of degree <= %d.',class(oldFE),p);
end

function E = exponents(d,p)
%EXPONENTS - Exponents of all monomials of total degree <= P
%   E = EXPONENTS(D,P) returns the exponents of all monomials in D <= 3
%   variables of total degree <= P, one row per monomial.
[a,b,c] = ndgrid(0:p);
E = unique([a(:),b(:),c(:)],'rows');
E = unique(E(:,1:d),'rows');
E = E(sum(E,2) <= p,:);
end

function V = dubiner(Y,p)
%DUBINER - Evaluate the Dubiner basis of degree P at the points Y
%   V = DUBINER(Y,P) evaluates the Dubiner basis of total degree P
%   (DGDk in 2d, DGDk3d in 3d) at the points Y (N-by-2 or N-by-3) and
%   returns one row per point and one column per basis function.
if size(Y,2) == 2
    V = DGDk.getElementBasis(p,Y);
else
    V = permute(DGDk3d.getElemBasis(p,Y,(p+1)*(p+2)*(p+3)/6),[3 1 2]);
end
V = reshape(V,size(Y,1),[]);
end

function L = physicalLocal(fe,feOld,mapsB,group,par)
%physicalLocal - Local L2 projections on the physical elements
%   L = physicalLocal(FE,FEOLD,MAPSB,GROUP,PAR) returns the local
%   matrices L(:,:,e) (FE.nrBasisFunc-by-FEOLD.nrBasisFunc-by-nrElems)
%   of the L2 projection of the basis of the ancestor PAR(e) onto the
%   basis of every element e of FE.mesh, computed on the physical
%   element with a Jacobi-scaled mass matrix. MAPSB and GROUP are as in
%   factoredLocal; the rows of the unused degrees of freedom
%   (nonDoFNodes) are zero.
meshF = fe.mesh;
n = size(meshF.n4e,2);
nE = meshF.nrElems;
nbF = fe.nrBasisFunc;
nbC = feOld.nrBasisFunc;
degree = 2*max(fe.degree,feOld.degree);
Bflat = mapsB(group,:);
[pC,C] = coarsePolynomials(feOld,n-1);
C = permute(C(:,par,:),[1 3 2]);
MB = integrate(@localProducts,meshF.n4e,1:nE,meshF,degree);
used = true(max([fe.nrDoF;double(fe.dof4e(:))]),1);
if isprop(fe,'nonDoFNodes')
    used(fe.nonDoFNodes) = false;
end
% unrefined elements (identity map), whose local matrices are the identity
% up to round-off
isId = all(abs(mapsB - reshape(eye(n),1,[])) < 1e-14,2);
L = zeros(nbF,nbC,nE);
for e = 1 : nE
    u = used(fe.dof4e(e,:));
    M = reshape(MB(e,:,1:nbF),nbF,nbF);
    Bm = reshape(MB(e,:,nbF+1:end),nbF,nbC);
    % Jacobi scaling that keeps M symmetric
    D = 1./sqrt(diag(M(u,u)));
    Lt = (M(u,u).*(D*D'))\(D.*Bm(u,:));
    if isId(group(e))
        % round-off removed column-wise
        Lt(abs(Lt) < 1e-10*max(abs(Lt),[],1)) = 0;
    end
    L(u,:,e) = D.*Lt;
end

    function val = localProducts(pts,pts_ref,idx)
        %localProducts - Products of fine with fine and coarse basis
        %   VAL = localProducts(PTS,PTS_REF,IDX) returns at the reference
        %   point PTS_REF (1-by-d) the products VAL(e,i,j) of the fine
        %   basis function i with the fine basis function j (j <= nbF)
        %   or with the coarse basis function j-nbF, summed over the
        %   components, for all elements e of FE.mesh.
        phiF = FormBasis.elementRows(fe.getBasis(pts,pts_ref,idx),idx,fe);
        phiF = reshape(phiF,nE,nbF,[]);
        lambda = [1-sum(pts_ref),pts_ref];
        lambdaParent = zeros(nE,n);
        for j = 1 : n
            lambdaParent(:,j) = Bflat(:,(j-1)*n+(1:n))*lambda(:);
        end
        Mon = dubiner(lambdaParent(:,2:end),pC);
        phiC = pagemtimes(reshape(Mon.',1,size(C,1),nE),C);
        phiC = reshape(permute(phiC,[3 2 1]),nE,nbC,[]);
        A = permute(phiF,[2 3 1]);
        val = permute(pagemtimes(A,permute(cat(2,phiF,phiC),[3 2 1])),[3 1 2]);
    end
end

function [p,C] = coarsePolynomials(feOld,d)
%coarsePolynomials - Dubiner coefficients of the coarse basis
%   [P,C] = coarsePolynomials(FEOLD,D) returns the lowest degree
%   P = FEOLD.degree,...,FEOLD.degree+4 whose Dubiner basis reproduces
%   the basis of FEOLD (checked at a test point; otherwise an error is
%   raised), and the coefficients C: C(m,c,k) is the coefficient of the
%   Dubiner function m in the basis value k (basis function and
%   component) of the coarse element c, in its reference coordinates.
nC = feOld.mesh.nrElems;
for p = feOld.degree + (0:4)
    E = exponents(d,p);
    Z = E/max(p,1);
    V = [];
    for l = 1 : size(Z,1)
        v = feOld.getBasis([],Z(l,:),[]);
        v = reshape(v,size(v,1),[]);
        if size(v,1) == 1
            v = repmat(v,nC,1);
        end
        V(l,:,:) = v; %#ok<AGROW>
    end
    C = reshape(dubiner(Z,p)\reshape(V,size(Z,1),[]),size(E,1),nC,[]);
    y = [0.21 0.17 0.13];
    y = y(1:d);
    b = feOld.getBasis([],y,[]);
    b = reshape(b,size(b,1),[]);
    if size(b,1) == 1
        b = repmat(b,nC,1);
    end
    f = reshape(sum(dubiner(y,p).'.*C,1),nC,[]);
    if norm(f(:) - b(:),inf) < 1e-9*max(1,norm(b(:),inf))
        return
    end
end
error('fem4m:FE:prolongationUnsupported', ...
    'No prolongation for %s: its basis is not a polynomial of degree <= %d.',class(feOld),p);
end