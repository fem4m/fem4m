function val = getEdgeValues4RefElemFunktion(refp,f,mesh)
%getEdgeValues4RefElemFunktion - Reference element function values on all edges
%   VAL = getEdgeValues4RefElemFunktion(REFP,F,MESH) evaluates the
%   function F at one point of every edge of the 2d or surface mesh
%   MESH, on both neighboring elements. The scalar REFP in [0,1] is the
%   parameter of the point along the edge as a local edge of the first
%   element MESH.e4ed(ED,1), 0 at its first and 1 at its second node;
%   only REFP(1) is used. F(XREF) returns the values of all elements at
%   the point XREF (1-by-2) of the reference triangle, as an array with
%   one row per element. VAL is nrEdges-by-...-by-2, the middle
%   dimensions those of F(XREF): VAL(ED,...,1) is the value of the
%   element MESH.e4ed(ED,1) and VAL(ED,...,2) the value of the element
%   MESH.e4ed(ED,2) at the same point (its edge parameter 1-REFP), zero
%   on boundary edges.
%
%   See also getEdgeValues4ElemFunktion, getFaceValues4RefElemFunktion,
%   FormValue, FormBasis, Mesh

% Copyright (c) 2026 by the fem4m authors

refp = [refp(1),0,1-refp(1),refp(1),0,1-refp(1)];

f4e1Ed1 = f([refp(1),refp(2)]);
f4e1Ed2 = f([refp(3),refp(4)]);
f4e1Ed3 = f([refp(5),refp(6)]);
f4e2Ed1 = f([1-refp(1),refp(2)]);
f4e2Ed2 = f([1-refp(3),1-refp(4)]);
f4e2Ed3 = f([refp(5),1-refp(6)]);

dims = size(f4e2Ed3);
f4e1Ed1 = reshape(f4e1Ed1,dims(1),[]);
f4e1Ed2 = reshape(f4e1Ed2,dims(1),[]);
f4e1Ed3 = reshape(f4e1Ed3,dims(1),[]);
f4e2Ed1 = reshape(f4e2Ed1,dims(1),[]);
f4e2Ed2 = reshape(f4e2Ed2,dims(1),[]);
f4e2Ed3 = reshape(f4e2Ed3,dims(1),[]);

edNumber = mesh.edNr4ed;
ind11 = find(edNumber(:,1)==1);
ind21 = find(edNumber(:,1)==2);
ind31 = find(edNumber(:,1)==3);
ind12 = find(edNumber(:,2)==1);
ind22 = find(edNumber(:,2)==2);
ind32 = find(edNumber(:,2)==3);

e4ed = mesh.e4ed;
val = zeros(size(e4ed,1),size(f4e1Ed1,2),2);
if ~isempty(ind11)
    val(ind11,:,1) = f4e1Ed1(e4ed(ind11,1),:);
end
if ~isempty(ind21)
    val(ind21,:,1) = f4e1Ed2(e4ed(ind21,1),:);
end
if ~isempty(ind31)
    val(ind31,:,1) = f4e1Ed3(e4ed(ind31,1),:);
end
if ~isempty(ind12)
    val(ind12,:,2) = f4e2Ed1(e4ed(ind12,2),:);
end
if ~isempty(ind22)
    val(ind22,:,2) = f4e2Ed2(e4ed(ind22,2),:);
end
if ~isempty(ind32)
    val(ind32,:,2) = f4e2Ed3(e4ed(ind32,2),:);
end

val = reshape(val,[mesh.nrEdges,dims(2:end),2]);

end