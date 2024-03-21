function spyEps(A,tol)
%spyEps - Display sparsity pattern ignoring tiny entries
%   spyEps(A) displays the sparsity pattern of A with spy, treating
%   entries with absolute value at most 2*eps as zero.
%
%   spyEps(A,TOL) uses the threshold TOL instead ([] for 2*eps).
%
%   See also spy, spfun

% Copyright (c) 2026 by the fem4m authors

if nargin < 2 || isempty(tol)
    tol = 2*eps;
end

A = spfun(@(z)spEps(z,tol),A);
spy(A);

end

function z = spEps(z,tol)
%spEps - Zero out tiny entries
%   Z = spEps(Z,TOL) sets the entries of Z with absolute value at most
%   TOL to zero.

z(abs(z) <= tol) = 0;

end