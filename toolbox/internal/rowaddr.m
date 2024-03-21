function val = rowaddr(A,indx,indy)
%ROWADDR - Entries of a matrix at row and column index pairs
%   VAL = ROWADDR(A,INDX,INDY) returns the entries A(INDX(i),INDY(i))
%   of the dense or sparse matrix A as a full array, computed by linear
%   indexing, not the submatrix A(INDX,INDY); e.g. rowaddr(N4E,E,K)
%   with column vectors E and K returns the global id of the local node
%   K(i) of the element E(i). INDX and INDY are expanded implicitly: a
%   column INDX and a row INDY give the matrix of the entries
%   A(INDX(i),INDY(j)). VAL has the size of the expanded indices, except
%   for a vector A and vector indices, where MATLAB linear indexing gives
%   VAL the orientation of A; a single-row A gives a column if INDX is a
%   column with more than one row.
%
%   See also sub2ind, sparse

% Copyright (c) 2026 by the fem4m authors

indz = uint64(indx) + size(A,1).*(uint64(indy)-1);

val = full(A(indz));

% A with a single row (e.g. one tetrahedron): A(INDZ) is a row
if size(val,1)~=size(indx,1) && size(A,1)==1
    val = val';
end

end