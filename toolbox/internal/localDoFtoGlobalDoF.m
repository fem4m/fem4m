function [I,J] = localDoFtoGlobalDoF(dof4e1,dof4e2)
%localDoFtoGlobalDoF - Global index pairs of local element dofs
%   I = localDoFtoGlobalDoF(DOF4E1) returns the row indices of the
%   local matrices of all elements, ordered like the entries of an
%   nrElems-by-n1-by-n1 array S of local matrices: the column vector I
%   holds DOF4E1(e,b) at the position of S(e,a,b). DOF4E1
%   (nrElems-by-n1) holds the global dof ids of the elements and may be
%   a cell array (joined with cell2mat).
%
%   I = localDoFtoGlobalDoF(DOF4E1,DOF4E2) takes the row indices from
%   DOF4E2 (nrElems-by-n2, numeric or cell array) and orders them like
%   an nrElems-by-n1-by-n2 array S: I holds DOF4E2(e,b) at the position
%   of S(e,a,b).
%
%   [I,J] = localDoFtoGlobalDoF(...) also returns the column indices J,
%   which hold DOF4E1(e,a) at the position of S(e,a,b), so that
%   sparse(I,J,S(:)) assembles the local matrices. I and J are double.
%
%   See also prolongationMatrix, BilinearForm, sparse

% Copyright (c) 2026 by the fem4m authors

if isa(dof4e1,'cell')
    dof4e1 = cell2mat(dof4e1);
end

dof4e1 = double(dof4e1);
if nargin < 2
    dof4e2 = dof4e1;
else
    if isa(dof4e2,'cell')
        dof4e2 = cell2mat(dof4e2);
    end
    dof4e2 = double(dof4e2);
end

I = permute(dof4e2,[1 3 2]).*ones(1,size(dof4e1,2));
I = I(:);

J = dof4e1(:).*ones(1,size(dof4e2,2));
J = J(:);

end