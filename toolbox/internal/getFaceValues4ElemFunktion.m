function val = getFaceValues4ElemFunktion(pts,parts,f,mesh)
%getFaceValues4ElemFunktion - Element function values on both elements of faces
%   VAL = getFaceValues4ElemFunktion(PTS,PARTS,F,MESH) evaluates the
%   function F on the two neighboring elements of the faces PARTS (a
%   vector of N face ids) of the 3d mesh MESH. PTS (N-by-3) holds one
%   point per face. F(X,ELEMS) returns the values at the points X (one
%   row per point) on the elements ELEMS, as an array with one row per
%   point. VAL has the size [size(F(PTS,E1)) 2] with
%   E1 = MESH.e4f(PARTS,1): VAL(:,...,1) holds the values on the first
%   element and VAL(:,...,2) those on the second element
%   MESH.e4f(PARTS,2), zero on boundary faces.
%
%   See also getFaceValues4RefElemFunktion, getEdgeValues4ElemFunktion,
%   FormValue, FormBasis, Mesh3d

% Copyright (c) 2026 by the fem4m authors

e4f1 = mesh.e4f(parts,1);
e4f2 = mesh.e4f(parts,2);
val1 = f(pts,e4f1);
originDim = size(val1);
val1 = val1(:,:);
val2 = zeros(size(val1));
if any(e4f2~=0)
    tmp = f(pts(e4f2~=0,:),e4f2(e4f2~=0));
    val2(e4f2~=0,:) = tmp(:,:);
end
val = [val1(:);val2(:)];
val = reshape(val,[originDim,2]);

end