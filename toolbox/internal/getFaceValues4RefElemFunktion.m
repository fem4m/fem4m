function val = getFaceValues4RefElemFunktion(refp,f,mesh)
%getFaceValues4RefElemFunktion - Reference element function values on all faces
%   VAL = getFaceValues4RefElemFunktion(REFP,F,MESH) evaluates the
%   function F at one point of every face of the 3d mesh MESH, on both
%   neighboring elements. REFP = [S T] (1-by-2) gives the point
%   P1 + S*(P2-P1) + T*(P3-P1) of the face with the nodes
%   MESH.n4f(K,:) = [P1 P2 P3] in increasing order of their node ids, so
%   that both elements evaluate the same point. F(XREF) returns the
%   values of all elements at the point XREF (1-by-3) of the reference
%   tetrahedron, as an array with one row per element. VAL is
%   nrFaces-by-...-by-2, the middle dimensions those of F(XREF):
%   VAL(K,...,L) is the value of the element MESH.e4f(K,L), zero for
%   L = 2 on boundary faces.
%
%   See also getFaceValues4ElemFunktion, getEdgeValues4RefElemFunktion,
%   FormValue, FormBasis, Mesh3d, FE

% Copyright (c) 2026 by the fem4m authors

faceNodes = [1 2 3; 1 2 4; 4 2 3; 1 3 4];
fperm = [ 1 2 3; 3 1 2; 2 3 1; 2 1 3; 1 3 2; 3 2 1];
refc = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
dims = [];

for j = 1 : 4
    [~,sortIdx] = sort(mesh.n4e(:,faceNodes(j,:)),2);
    [~,perm] = ismember([0 0 0;sortIdx],[0 0 0;fperm],'rows');
    for k = 1 : 6
        z = (refc(faceNodes(j,fperm(k,1)),:) + ...
            refp(1)*(refc(faceNodes(j,fperm(k,2)),:)-refc(faceNodes(j,fperm(k,1)),:)) + ...
            refp(2)*(refc(faceNodes(j,fperm(k,3)),:)-refc(faceNodes(j,fperm(k,1)),:)) );
        for l = 1 : 2
            indx = mesh.fNr4f(:,l) == j & perm(mesh.e4f(:,l)+1)-1 == k;
            if isempty(dims)
                fz = f(z);
                dims = size(fz);
                fz = reshape(fz,dims(1),[]);
                val = zeros(mesh.nrFaces,prod(dims(2:end)),2);
            else
                fz = reshape(f(z),dims(1),[]);
            end
            val(indx,:,l) = fz(mesh.e4f(indx,l),:);
        end
    end

end

val = reshape(val,[mesh.nrFaces,dims(2:end),2]);

end