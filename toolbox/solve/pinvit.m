function [x,lambda,iter] = pinvit(x,A,B,tol,mg,mgIter,maxIter)
%PINVIT - Inexact preconditioned inverse iteration
%   X = PINVIT(X0,A,B,TOL,MG,mgIter,maxIter) approximates an
%   eigenvector X of the generalized eigenvalue problem
%   A*X = LAMBDA*B*X for symmetric positive definite A and B by
%   preconditioned inverse iteration from the start vector X0,
%   usually for the smallest eigenvalue. Each step sets
%   X = X - Z, normalized to X'*B*X = 1, where Z is the result of
%   at most mgIter cycles of the Multigrid solver MG for the matrix
%   A (its finest level) applied to the residual
%   R = A*X - LAMBDA*B*X of the Rayleigh quotient LAMBDA. The
%   iteration stops when sqrt(2*R'*Z) <= TOL or after maxIter-1
%   steps. A and B are the matrices of the finest level of MG
%   restricted to its free entries, e.g. A(fe.freeNodes,fe.freeNodes),
%   and X0 is a column vector of the matching size.
%
%   [X,LAMBDA] = PINVIT(...) also returns the Rayleigh quotient
%   LAMBDA of X.
%
%   [X,LAMBDA,ITER] = PINVIT(...) also returns ITER, the number of
%   steps plus 1 (ITER = 1 when X0 already meets TOL).
%
%   See also eigs, Multigrid

% Copyright (c) 2026 by the fem4m authors

iter = 1;
x = x./sqrt(x'*B*x);
Ax = A*x;
Bx = B*x;
lambda = (x'*Ax)/(x'*Bx);
r = Ax - lambda*Bx;
z = mg.solve(r,0,mgIter,[],0);
eNew = sqrt(2*r'*z);
while (eNew > tol) && (iter < maxIter)
    x = x - z;
    x = x./sqrt(x'*B*x);
    Ax = A*x;
    Bx = B*x;
    lambda = (x'*Ax)/(x'*Bx);
    r = Ax - lambda*Bx;
    z = mg.solve(r,0,mgIter,[],0);
    eNew = sqrt(2*r'*z);
    iter = iter + 1;
end
end