function val = integrate(integrand,parts,indices,mesh,degree)
%INTEGRATE - Vectorized Gauss integration over mesh parts
%   VAL = INTEGRATE(INTEGRAND,PARTS,INDICES,MESH,DEGREE) integrates the
%   function handle INTEGRAND over the parts INDICES of MESH (edges,
%   faces or elements) with a Gauss rule that is exact for polynomials
%   of degree DEGREE on affine parts, in a vectorized loop over the
%   Gauss points. PARTS holds the node ids of these parts, one row per
%   entry of INDICES (e.g. MESH.n4e(INDICES,:)). Its number of columns
%   selects the reference part and the weights:
%       2 - edges, interval [0,1], length4ed
%       3 - triangles of a 2d mesh (area4e) or faces of a 3d mesh
%           (area4f), reference triangle (0,0),(1,0),(0,1)
%       4 - tetrahedra of a 3d mesh (area4e), reference tetrahedron,
%           or quadrilaterals of a 2d mesh, unit square
%   The columns of PARTS are mapped in this order to the vertices
%   (0,0),(1,0),(0,1) of the triangle, (0,0,0),(1,0,0),(0,1,0),(0,0,1)
%   of the tetrahedron and (0,0),(1,0),(1,1),(0,1) of the square
%   (bilinear map; exact only for parallelograms, since the weight is
%   the constant area4e). The rule has ceil((DEGREE+1)/2) Gauss points
%   in every direction (a conical product rule on triangles and
%   tetrahedra). integrate raises an error for other numbers of
%   columns.
%
%   INTEGRAND(PTS,PTSREF,IDX) is called once per Gauss point with the
%   point PTSREF (1-by-d) of the reference part, its images PTS in the
%   parts (one row per part) and the ids IDX = INDICES of the parts. It
%   returns an array with one row per part. VAL is the weighted sum of
%   these arrays and has the same size: VAL(i,...) is the integral over
%   the part INDICES(i).
%
%   On a MeshSurface, PTS and the weights are point4ed and length4ed
%   (PARTS with two columns) or point4e and area4e of MESH at PTSREF,
%   i.e. on the curved surface for a surface degree > 1; otherwise only
%   the number of columns of PARTS is used.
%
%   On a MeshPoly with polygonal elements, PARTS is the cell array of
%   the node lists of the polygons, MESH.n4e(INDICES). Every polygon is
%   integrated on a fan of triangles from its point mid4e (MESH.convex
%   true) or from its first node (MESH.convex false). INTEGRAND is then
%   called with PTSREF = [] and the ids IDX of the polygons in INDICES
%   with the same number of nodes, and returns one row per entry of
%   IDX.
%
%   See also BilinearForm, LinearForm, ErrorForm, EstimatorForm

% Copyright (c) 2026 by the fem4m authors

if size(parts,2) ~= 2 && isa(mesh,'MeshPoly')
    if mesh.convex
        val = integratePolyMid(integrand,parts,indices,mesh,degree);
    else
        val = integratePoly(integrand,parts,indices,mesh,degree);
    end
elseif isa(mesh,'MeshSurface')
    val = integrateSurface(integrand,parts,indices,mesh,degree);
else
    val = integrateMesh(integrand,parts,indices,mesh,degree);
end

end

%% integrate mesh
function val = integrateMesh(integrand,parts,indices,mesh,degree)
%integrateMesh - Quadrature over mesh elements
%   VAL = integrateMesh(INTEGRAND,PARTS,INDICES,MESH,DEGREE) integrates
%   INTEGRAND over the edges, triangles, faces, quadrilaterals or
%   tetrahedra given by the node table PARTS, see integrate.

nrGaussPoints = ceil((degree+1)/2);

% 1d or 2d or 3d
if size(parts,2) == 2
    % 1d line
    length4ed = mesh.length4ed;
    [refGaussPoints,refGaussWeights] = getGaussPoints(nrGaussPoints);
    weightFactor = length4ed(indices);
elseif size(parts,2) == 3
    % 2d triangle (element in 2d, face in 3d)
    [refGaussPoints,refGaussWeights] = getConProdGaussPoints(nrGaussPoints);
    if size(mesh.c4n,2) == 2
        area4e = mesh.area4e();
        weightFactor = 2*area4e(indices);
    else
        area4f = mesh.area4f();
        weightFactor = 2*area4f(indices);
    end
elseif (size(parts,2) == 4) && (size(mesh.c4n,2) == 2)
    % 2d quadrilateral (reference unit square)
    area4e = mesh.area4e;
    [refGaussPoints,refGaussWeights] = getCompositeGaussPoints(nrGaussPoints);
    weightFactor = area4e(indices);
elseif (size(parts,2) == 4) && (size(mesh.c4n,2) == 3)
    % 3d tetrahedron
    area4e = mesh.area4e;
    [refGaussPoints,refGaussWeights] = getConProdGaussPoints3D(nrGaussPoints);
    weightFactor = 6*area4e(indices);
else
    error('fem4m:integrate:unknownDomain','Could not integrate because no Domain is specified.');
end

% Integrate
for curPoint = 1 : numel(refGaussWeights)
    %transformation
    curCoords4e = ref2arbitrary(parts,refGaussPoints(curPoint,:),mesh.c4n);
    % evaluate function handle
    integrandVal = integrand(curCoords4e,refGaussPoints(curPoint,:),indices);
    if curPoint == 1; val = zeros(size(integrandVal)); end
    % weighted sum
    val = val + (refGaussWeights(curPoint)*weightFactor).*integrandVal;
end

end

%% integrate surface
function val = integrateSurface(integrand,parts,indices,mesh,degree)
%integrateSurface - Quadrature over surface elements
%   VAL = integrateSurface(INTEGRAND,PARTS,INDICES,MESH,DEGREE)
%   integrates INTEGRAND over the edges (PARTS with two columns) or the
%   elements INDICES of the MeshSurface MESH with the points and the
%   metric of its parametrization, see integrate.

nrGaussPoints = ceil((degree+1)/2);

if size(parts,2) == 2 % 1d line
    [refGaussPoints,refGaussWeights] = getGaussPoints(nrGaussPoints);
else
    % 2d surface triangle
    [refGaussPoints,refGaussWeights] = getConProdGaussPoints(nrGaussPoints);
end

for curPoint = 1 : numel(refGaussWeights)
    % surface transformation
    if size(parts,2) == 2
        curCoords4e = mesh.point4ed(refGaussPoints(curPoint,:));
    else
        curCoords4e = mesh.point4e(refGaussPoints(curPoint,:));
    end
    curCoords4e = curCoords4e(indices,:);
    % evaluate function handle
    integrandVal = integrand(curCoords4e,refGaussPoints(curPoint,:),indices);
    if curPoint == 1; val = zeros(size(integrandVal)); end
    % surface weights
    if size(parts,2) == 2
        weightFactor = mesh.length4ed(refGaussPoints(curPoint,:));
    elseif size(parts,2) == 3
        weightFactor = 2*mesh.area4e(refGaussPoints(curPoint,:));
    else
        weightFactor = mesh.area4e(refGaussPoints(curPoint,:));
    end
    weightFactor = weightFactor(indices);
    % weighted sum
    val = val + (refGaussWeights(curPoint)*weightFactor).*integrandVal;
end

end

%% integrate poly midpoint
function val = integratePolyMid(integrand,parts,indices,mesh,degree)
%integratePolyMid - Quadrature over convex polygons
%   VAL = integratePolyMid(INTEGRAND,PARTS,INDICES,MESH,DEGREE)
%   integrates INTEGRAND over the polygons INDICES of the MeshPoly MESH,
%   with the node lists PARTS (cell array), on the triangles that
%   connect every edge with the point mid4e of the polygon.

nrN4e = mesh.nrN4e(indices);
[refGaussPoints,refGaussWeights] = getConProdGaussPoints(ceil((degree+1)/2));
mid4e = mesh.mid4e;
c4n = [mesh.c4n;mid4e(indices,:)];
valIsEmpty = true;

for j = unique(nrN4e)'
    curParts = cell2mat(parts(nrN4e==j));
    for k = 1:j
        curN4e = [curParts(:,[k mod(k,j)+1]), mesh.nrNodes+find(nrN4e==j)];
        area4e = computeArea4e(c4n,curN4e);
        for curPoint = 1 : numel(refGaussWeights)
            % points and weights of subtriangulation
            points = ref2arbitrary(curN4e,refGaussPoints(curPoint,:),c4n);
            weights = refGaussWeights(curPoint)*2*area4e;
            % evaluate function handle
            integrandVal = integrand(points,[],indices(nrN4e==j));
            if valIsEmpty
                dim = size(integrandVal);
                dim(1) = numel(indices);
                val = zeros(dim);
                valIsEmpty = false;
            end
            % weighted sum
            val(nrN4e==j,:,:,:) = val(nrN4e==j,:,:,:) + weights.*integrandVal;
        end
    end
end

end

%% integrate nonconvex poly
function val = integratePoly(integrand,parts,indices,mesh,degree)
%integratePoly - Quadrature over nonconvex polygons
%   VAL = integratePoly(INTEGRAND,PARTS,INDICES,MESH,DEGREE) integrates
%   INTEGRAND over the polygons INDICES of the MeshPoly MESH, with the
%   node lists PARTS (cell array), on the fan of triangles from the
%   first node of every polygon, weighted with their signed areas.
% connectingNode = 1 !
nrN4e = mesh.nrN4e(indices);
[refGaussPoints,refGaussWeights] = getConProdGaussPoints(ceil((degree+1)/2));
valIsEmpty = true;

for j = unique(nrN4e)'
    curParts = cell2mat(parts(nrN4e==j));
    for k = 1:j-2
        curN4e = curParts(:,[1,mod(k,j)+1,mod(k+1,j)+1]);
        area4e = computeArea4e(mesh.c4n,curN4e);
        for curPoint = 1 : numel(refGaussWeights)
            % points and weights of subtriangulation
            points = ref2arbitrary(curN4e,refGaussPoints(curPoint,:),mesh.c4n);
            weights = refGaussWeights(curPoint)*2*area4e;
            % evaluate function handle
            integrandVal = integrand(points,[],indices(nrN4e==j));
            if valIsEmpty
                dim = size(integrandVal);
                dim(1) = numel(indices);
                val = zeros(dim);
                valIsEmpty = false;
            end
            % weighted sum
            val(nrN4e==j,:,:,:) = val(nrN4e==j,:,:,:) + weights.*integrandVal;
        end
    end
end

end

%% ref2arbitrary
function val = ref2arbitrary(parts,curPoint,c4n)
%ref2arbitrary - Map reference points to elements
%   VAL = ref2arbitrary(PARTS,CURPOINT,C4N) maps the point CURPOINT
%   (1-by-d) of the reference part to every part of the node table PARTS
%   with the node coordinates C4N and returns one row per row of PARTS.

%1d, 2d or 3d
if size(parts,2) == 2
    %1d ref = conv(0,1)
    c1 = c4n(parts(:,1),:);
    c2 = c4n(parts(:,2),:);

    val = c1 + curPoint(1)*(c2-c1);

elseif size(parts,2) == 3
    % 2d ref = conv{ (0,0),(1,0),(0,1)}
    c1 = c4n(parts(:,1),:);
    c2 = c4n(parts(:,2),:);
    c3 = c4n(parts(:,3),:);

    val = c1 + curPoint(1)*(c2-c1) + ...
        curPoint(2)*(c3-c1);

elseif size(parts,2) == 4 && size(c4n,2) == 2
    % 2d ref = unit square, PARTS columns at (0,0),(1,0),(1,1),(0,1)
    c1 = c4n(parts(:,1),:);
    c2 = c4n(parts(:,2),:);
    c3 = c4n(parts(:,3),:);
    c4 = c4n(parts(:,4),:);

    val = c1 + curPoint(1)*(c2-c1) + ...
        curPoint(2)*(c4-c1 + curPoint(1)*(c3-c4-c2+c1));

elseif size(parts,2) == 4 && size(c4n,2) == 3
    % 3d ref = conv{(0,0,0),(1,0,0),(0,1,0),(0,0,1)}
    c1 = c4n(parts(:,1),:);
    c2 = c4n(parts(:,2),:);
    c3 = c4n(parts(:,3),:);
    c4 = c4n(parts(:,4),:);
    val = c1 + curPoint(1)*(c2-c1) + curPoint(2)*(c3-c1) + curPoint(3)*(c4-c1);
end
end

function area4e = computeArea4e(c4n,n4e)
%computeArea4e - Areas of triangle rows
%   AREA4E = computeArea4e(C4N,N4E) returns the signed areas of the
%   triangles N4E with the node coordinates C4N, positive for
%   counterclockwise vertices.
x1 = c4n(n4e(:,1),1);
x2 = c4n(n4e(:,2),1);
x3 = c4n(n4e(:,3),1);
y1 = c4n(n4e(:,1),2);
y2 = c4n(n4e(:,2),2);
y3 = c4n(n4e(:,3),2);
area4e = ( x1.*(y2 - y3) + x2.*(y3 - y1) + x3.*(y1 - y2) )/2;
end

%% Gauss points

function [x,w] = getConProdGaussPoints(n)
%getConProdGaussPoints - Conical product rule on triangle
%   [X,W] = getConProdGaussPoints(N) returns the N^2 points X (N^2-by-2)
%   and weights W of the conical product rule on the reference triangle
%   (0,0),(1,0),(0,1), exact for polynomials of degree 2N-1; the
%   weights sum to 1/2.

%% Find n Gauss_Legendre Points for Interval [-1,1]
gamma = (1 : n-1) ./ sqrt(4*(1 : n-1).^2 - ones(1,n-1) );
[V,D] = eig( diag(gamma,1) + diag(gamma,-1) );
r = diag(D);
a = 2*V(1,:).^2;% norm-factor -> int(1,-1,1)=2

%% Find n Gauss_Jacobi Points for Interval [-1,1]
delta = -1./(4*(1 : n).^2-ones(1,n));
gamma = sqrt((2 : n).*(1 : n-1)) ./ (2*(2 : n)-ones(1,n-1));
[V,D] = eig( diag(delta)+diag(gamma,1)+diag(gamma,-1) );
s = diag(D);
b = 2*V(1,:).^2; % norm-factor -> int((1-x),-1,1)=2

%% linear map to Interval [0,1]
% w(x) = 1 changes norm-factor from 2 to 1
r = .5 * r + .5;
a = .5 * a';
% w(x) = (1-x) changes norm-factor from 2 to 1/2
s = .5 * s + .5;
b = .25 * b';

%% conical product [ s_j , r_i(1-s_j) ]  a_i*b_j
s = repmat(s',n,1); s = s(:);
r = repmat(r,n,1);
x = [ s , r.*(ones(n^2,1)-s) ];
w = a*b';
w = w(:);
end

function [x,w] = getGaussPoints(n)
%getGaussPoints - Gauss-Legendre rule on unit interval
%   [X,W] = getGaussPoints(N) returns the N Gauss-Legendre points X
%   (N-by-1) and weights W on [0,1], exact for polynomials of degree
%   2N-1; the weights sum to 1.

%% Find n Gauss_Legendre Points for Interval [-1,1]
gamma = (1 : n-1) ./ sqrt(4*(1 : n-1).^2 - ones(1,n-1) );
[V,D] = eig( diag(gamma,1) + diag(gamma,-1) );
x = diag(D);
w = 2*V(1,:).^2;

%% linear map to Interval [0,1]
x = .5 * x + .5;
w = .5 * w';
end

function [x,w] = getCompositeGaussPoints(n)
%getCompositeGaussPoints - Tensor product Gauss rule on unit square
%   [X,W] = getCompositeGaussPoints(N) returns the N^2 points X
%   (N^2-by-2) and weights W of the tensor product Gauss-Legendre rule
%   on [0,1]^2, exact for polynomials of degree 2N-1 in each variable;
%   the weights sum to 1.

%% Find n Gauss_Legendre Points for Interval [-1,1]
gamma = (1 : n-1) ./ sqrt(4*(1 : n-1).^2 - ones(1,n-1) );
[V,D] = eig( diag(gamma,1) + diag(gamma,-1) );
x = diag(D);
w = 2*V(1,:).^2;

%% linear map to square [0,1] times [0,1]
x1 = .5*x + .5;
x2 = .5*x + .5;
xi1 = repmat(x1',n,1);
x = [xi1(:),repmat(x2,n,1)];
w = w'*w;
w = .25 *w(:);
end

function [x,w] = getConProdGaussPoints3D(n)
%getConProdGaussPoints3D - Conical product rule on tetrahedron
%   [X,W] = getConProdGaussPoints3D(N) returns the N^3 points X
%   (N^3-by-3) and weights W of the conical product rule on the
%   reference tetrahedron (0,0,0),(1,0,0),(0,1,0),(0,0,1), exact for
%   polynomials of degree 2N-1; the weights sum to 1/6.

% Find n Gauss_Legendre Points for Interval [-1,1]
gamma = (1 : n-1) ./ sqrt(4*(1 : n-1).^2 - ones(1,n-1) );
[V,D] = eig( diag(gamma,1) + diag(gamma,-1) );
r = diag(D);
a = 2*V(1,:).^2;% norm-factor -> int(1,-1,1)=2

% Find n Gauss_Jacobi Points for Interval [-1,1]
delta = -1./(4*(1 : n).^2-ones(1,n));
gamma = sqrt((2 : n).*(1 : n-1)) ./ (2*(2 : n)-ones(1,n-1));
[V,D] = eig( diag(delta)+diag(gamma,1)+diag(gamma,-1) );
s = diag(D);
b = 2*V(1,:).^2; % norm-factor -> int((1-x),-1,1)=2

% Find n Gauss_Jacobi Points for Interval [-1,1]
delta = -1./((1:n).*(2:n+1));
gamma = sqrt(((2:n).^2-1).^2./((2:n).^2.*(4*(2:n).^2-1)));
[V,D] = eig( diag(delta)+diag(gamma,1)+diag(gamma,-1) );
t = diag(D);
c = 8/3*V(1,:).^2; % norm-factor -> int((1-x)^2,-1,1)=8/3

% linear map to Interval [0,1]
% w(x) = 1 changes norm-factor from 2 to 1
r = .5 * r + .5;
a = .5 * a';
% w(x) = (1-x) changes norm-factor from 2 to 1/2
s = .5 * s + .5;
b = .25 * b';
% w(x) = (1-x)^2 changes norm-factor from 8/3 to 1/3
t = .5 * t + .5;
c = .125 * c';

% conical product [r_i(1-s_j)(1-t_k),t_k, s_j(1-t_k)]  a_i*b_j*c_k
r = r*ones(1,n^2);r=r';r=r(:);
s = s*ones(1,n);s=s';s=s(:)*ones(1,n);s=s(:);
t = t*ones(1,n^2);t=t(:);

a = a*ones(1,n^2);a=a';a=a(:);
b = b*ones(1,n);b=b';b=b(:)*ones(1,n);b=b(:);
c = c*ones(1,n^2);c=c(:);

x = [r.*(1-s).*(1-t),t,s.*(1-t)];
w = a.*b.*c;
end