function val = matMul(A,B)
%matMul - Elementwise matrix product along first dimension
%   VAL = matMul(A,B) returns the N-by-P-by-R array of the matrix
%   products VAL(k,:,:) = A(k,:,:)*B(k,:,:) of the N-by-P-by-Q array A
%   and the N-by-Q-by-R array B, where A(k,:,:) is read as a P-by-Q and
%   B(k,:,:) as a Q-by-R matrix. A or B may have a single row (size 1
%   in the first dimension), which is expanded implicitly. An N-by-P
%   matrix is an N-by-P-by-1 array. matMul raises an error if
%   size(A,3) differs from size(B,2).
%
%   See also pagemtimes, mtimes

% Copyright (c) 2026 by the fem4m authors

if size(A,3) ~= size(B,2)
    dimA = size(A);
    if numel(dimA)==2; dimA(3)=1; end
    dimB = size(B);
    if numel(dimB)==2; dimB(3)=1; end
    error('fem4m:matMul:dimensionMismatch', ...
        'matMul: dimension mismatch for matrix product A(k,:,:)*B(k,:,:) ([%s]*[%s])', ...
        num2str(dimA(2:end)),num2str(dimB(2:end)));
end

val = sum(permute(A,[1 2 4 3]).*permute(B,[1 4 3 2]),4);