function [x,lambda,k,iter] = pdaSet(A,b,low,up,iterMax,x,isVertex,display,P,tol,finaltol,fe,iterGMRES,restartGMRES)
%pdaSet - Primal-dual active set algorithm
%   X = pdaSet(A,B) solves the quadratic problem
%   min 1/2*X'*A*X - B'*X subject to LOW <= X <= UP for a symmetric
%   positive definite matrix A with the primal-dual active set
%   method (cf. Hintermueller et al. 2003) without bounds (default
%   LOW = -Inf, UP = Inf). The solution satisfies A*X = B on the
%   inactive entries, and LAMBDA = A*X - B >= 0 where X = LOW and
%   LAMBDA <= 0 where X = UP. The active sets are predicted with the
%   constant c = 1e3, and the iteration stops when they no longer
%   change. The systems on the inactive entries are solved with
%   backslash.
%
%   X = pdaSet(A,B,LOW) also specifies the lower bounds LOW, a
%   column vector with one entry per row of A (default -Inf).
%
%   X = pdaSet(A,B,LOW,UP) also specifies the upper bounds UP, a
%   column vector with one entry per row of A (default Inf).
%
%   X = pdaSet(A,B,LOW,UP,iterMax) also specifies the maximum
%   number of active set iterations (default 1000). When it is
%   reached without convergence, a message is printed at every
%   DISPLAY level and the last iterate is returned.
%
%   X = pdaSet(A,B,...,X0) also specifies the initial iterate X0,
%   which defines the first active sets (default zeros).
%
%   X = pdaSet(A,B,...,isVertex) constrains only the entries where
%   the mask isVertex is true (default all entries).
%
%   X = pdaSet(A,B,...,DISPLAY) also selects the output: 0 (default)
%   none, 1 an iteration counter, 2 the counter without P and the
%   GMRES messages with P, 3 also the residuals and the contact
%   error of each step and the number of iterations, 4 also a plot
%   of the active set in figure 1 in each step.
%
%   X = pdaSet(A,B,...,P) solves the systems on the inactive entries
%   with GMRES preconditioned by the function handle P, which
%   applies an approximation of the inverse of A to a column vector
%   with one entry per row of A (zero on the active entries); only
%   the inactive entries of its result are used. A final GMRES
%   solve with tolerance FINALTOL follows the active set iteration.
%
%   X = pdaSet(A,B,...,TOL) also specifies the tolerance of the
%   GMRES solves in the active set iteration (default 1e-6).
%
%   X = pdaSet(A,B,...,FINALTOL) also specifies the tolerance of
%   the final GMRES solve (default eps).
%
%   X = pdaSet(A,B,...,FE) also specifies the finite element space
%   of the plot, with X on the entries FE.freeNodes; it is required
%   for DISPLAY > 3 (error fem4m:pdaSet:missingFE).
%
%   X = pdaSet(A,B,...,iterGMRES) also specifies the maximum number
%   of outer GMRES iterations (default 100).
%
%   X = pdaSet(A,B,...,restartGMRES) also specifies the number of
%   inner GMRES iterations between restarts (default 100). Both GMRES
%   limits are capped by the number of inactive entries. TOL,
%   FINALTOL, iterGMRES and restartGMRES are used only with P. Empty
%   arguments select the defaults.
%
%   [X,LAMBDA] = pdaSet(...) also returns the multipliers
%   LAMBDA = A*X - B, zero on the inactive and the unconstrained
%   entries. With P, LAMBDA is computed before the final GMRES solve
%   and holds only up to TOL.
%
%   [X,LAMBDA,K] = pdaSet(...) also returns the number K of active
%   set iterations.
%
%   [X,LAMBDA,K,ITER] = pdaSet(...) also returns ITER, the sum of the
%   GMRES iteration counts [outer inner] (1-by-2) of the active set
%   iteration. It is 0 without P, with DISPLAY > 1 or when no GMRES
%   solve was needed, and does not count the final GMRES solve.
%
%   See also gmres, exObstacleP1Adaptive

% Copyright (c) 2026 by the fem4m authors

k = 0;
iter = 0;
c = 1e3;

if nargin < 3 || isempty(low)
    low = -inf(size(A,1),1);
end

if nargin < 4 || isempty(up)
    up = inf(size(A,1),1);
end

if nargin < 5 || isempty(iterMax)
    iterMax = 1000;
end

if nargin < 6 || isempty(x)
    x = zeros(size(A,1),1);
end

if nargin < 7 || isempty(isVertex)
    isVertex = ones(size(A,1),1);
end

if nargin < 8 || isempty(display)
    display = 0;
end

if nargin < 9 || isempty(P)
    P = [];
end

if nargin < 10 || isempty(tol)
    tol = 1e-6;
end

if nargin < 11 || isempty(finaltol)
    finaltol = eps;
end

if nargin < 12 || isempty(fe)
    if display > 3
        error('fem4m:pdaSet:missingFE','For display>3 a FE needs to be specified for plotting.');
    end
end

if nargin < 13 || isempty(iterGMRES)
    iterGMRES = 100;
end

if nargin < 14 || isempty(restartGMRES)
    restartGMRES = 100;
end

lambda = A*x - b;
activeLb = (lambda + c*(low - x)) > 0;
activeUb = (lambda + c*(up - x)) < 0;
% only vertices are constrained
activeLb(~isVertex) = false;
activeUb(~isVertex) = false;
active = activeLb | activeUb;
activeOld = ~active;

while(k < iterMax && any(active~=activeOld))
    % A_ii*u_i = b_i - A_il*low_l - A_iu*up_u
    % u_l = low_l
    % u_u = up_u
    % lambda = A*u - b
    % lambda_i = 0
    if nnz(~active) > 0
        if ~isempty(P)
            rhs = b(~active)-A(~active,activeLb)*low(activeLb)-A(~active,activeUb)*up(activeUb);
            if display > 1
                x(~active) = gmres(A(~active,~active),rhs,min(nnz(~active),restartGMRES),...
                    tol,min(nnz(~active),iterGMRES),@(z)(applyP(z,P,active)),[],x(~active));
            else
                [x(~active),~,~,it] = gmres(A(~active,~active),rhs,min(nnz(~active),restartGMRES),...
                    tol,min(nnz(~active),iterGMRES),@(z)(applyP(z,P,active)),[],x(~active));
                iter = iter + it;
            end
        else
            x(~active) = A(~active,~active)\(b(~active)-A(~active,activeLb)*low(activeLb)-A(~active,activeUb)*up(activeUb));
        end
    end
    x(activeLb) = low(activeLb);
    x(activeUb) = up(activeUb);
    lambda = A*x - b;
    lambda(~active) = 0;

    activeOld = active;
    activeLbOld = activeLb;
    activeUbOld = activeUb;
    activeLb = (lambda + c*(low - x)) > 0;
    activeUb = (lambda + c*(up - x)) < 0;
    % only vertices are constrained
    activeLb(~isVertex) = false;
    activeUb(~isVertex) = false;
    active = activeLb | activeUb;

    k = k + 1;

    % diagnostic output
    if display==1 || (display > 1 && isempty(P))
        if k > 1
            fprintf(repmat('\b',1,numel(int2str(k-1))));
        end
        fprintf('%d \b', k);
    end
    if display > 2
        psi = zeros(numel(x),1);
        psi(activeLbOld) = low(activeLbOld);
        psi(activeUbOld) = up(activeUbOld);
        disp([' local res:', num2str(norm(A(~activeOld,~activeOld)*x(~activeOld)-(b(~activeOld)-A(~activeOld,activeLbOld)*low(activeLbOld)-A(~activeOld,activeUbOld)*up(activeUbOld)))),...
            ' global res:',num2str(norm(A*x-lambda-b)),...
            ' contact error: ', num2str(norm(x(activeOld) - psi(activeOld))),...
            ' lambda(I): ', num2str(norm(lambda(~activeOld)))]);
    end
    if display > 3
        figure(1);clf;
        plotActive = zeros(numel(fe.fixedNodes)+numel(fe.freeNodes),1);
        plotActive(fe.freeNodes(active)) = -1;
        fe.plot(plotActive);
        shading interp;
        colormap(gray);
        pause(0.5);
    end
end

% final GMRES solve with tolerance finaltol (only with P)
if ~isempty(P)
    rhs = b(~active)-A(~active,activeLb)*low(activeLb)-A(~active,activeUb)*up(activeUb);
    if display > 1
        x(~active) = gmres(A(~active,~active),rhs,min(nnz(~active),restartGMRES),...
            finaltol,min(nnz(~active),iterGMRES),@(z)(applyP(z,P,active)),[],x(~active));
    else
        [x(~active),~,~,~] = gmres(A(~active,~active),rhs,min(nnz(~active),restartGMRES),...
            finaltol,min(nnz(~active),iterGMRES),@(z)(applyP(z,P,active)),[],x(~active));
    end
end

if display == 1 || (display > 1 && isempty(P))
    fprintf(repmat('\b',1,numel(int2str(k))));
end
if display > 2
    disp(['Iterations = ',num2str(k)]);
end

if nnz(active~=activeOld)>0
    disp('Stopped without converging (max iterations)');
end

end

function val = applyP(z,P,active)
%applyP - Apply the preconditioner on the inactive entries
%   VAL = applyP(Z,P,ACTIVE) extends Z by zeros on the entries of
%   the logical mask ACTIVE, applies the function handle P, and
%   returns the inactive entries of the result.
val = zeros(numel(active),1);
val(~active) = z;
val = P(val);
val = val(~active);
end