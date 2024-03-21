function val = getEdgeValues4ElemFunktion(pts,parts,f,mesh)
%getEdgeValues4ElemFunktion - Element function values on both elements of edges
%   VAL = getEdgeValues4ElemFunktion(PTS,PARTS,F,MESH) evaluates the
%   function F on the two neighboring elements of the edges PARTS (a
%   vector of N edge ids) of the 2d, polygonal or surface mesh MESH.
%   PTS holds one point per edge (one row per edge). F(X,ELEMS) returns
%   the values at the points X (one row per point) on the elements
%   ELEMS, as an array with one row per point. VAL has the size
%   [size(F(PTS,E1)) 2] with E1 = MESH.e4ed(PARTS,1): VAL(:,...,1)
%   holds the values on the first element and VAL(:,...,2) those on
%   the second element MESH.e4ed(PARTS,2), zero on boundary edges.
%
%   See also getEdgeValues4RefElemFunktion, getFaceValues4ElemFunktion,
%   FormValue, FormBasis, Mesh

% Copyright (c) 2026 by the fem4m authors

e4ed1 = mesh.e4ed(parts,1);
e4ed2 = mesh.e4ed(parts,2);
val1 = f(pts,e4ed1);
originDim = size(val1);
val1 = val1(:,:);
val2 = zeros(size(val1));
if any(e4ed2~=0)
    tmp = f(pts(e4ed2~=0,:),e4ed2(e4ed2~=0));
    val2(e4ed2~=0,:) = tmp(:,:);
end
val = [val1(:);val2(:)];
val = reshape(val,[originDim,2]);

end