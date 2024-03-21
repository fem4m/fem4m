%[text] # Smoothing property of the Jacobi and Richardson iterations
%[text] Applies the damped Jacobi iteration to the linear system of the $P\_1$ discretization of the Poisson problem and expands its error in the eigenvectors of the iteration. The oscillating components of the error vanish after a few steps, while the smooth components hardly change: the iteration is a poor solver but a good smoother, the first ingredient of the two-grid and multigrid methods.
%[text] **Builds on:** `exPoissonP1`. **Next:** `exSolverTwoGrid`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] On the unit square $\\Omega = (0,1)^2$ the $P\_1$ discretization of $-\\Delta u = f$ in $\\Omega$, $u = 0$ on $\\partial\\Omega$ (`exPoissonP1`) leads to the linear system $Ax = b$ for the coefficients of the $M$ interior nodes. The problem object `problemScalarEllipticSquareOscExact(10)` provides the source of the oscillating exact solution $u = \\sin(10\\pi x\_1)\\sin(10\\pi x\_2)$. Uniform refinement of the coarse mesh with two triangles gives right isosceles triangles whose nodes form a uniform grid with the spacing $h$. The iterations below approximate the discrete solution $x = A^{-1}b$.
%[text] ### Jacobi and Richardson iteration
%[text] Starting from $x\_0$, both iterations compute
%[text]{"align":"center"} $ x\_{n+1} = x\_n + \\omega D^{-1}(b - Ax\_n), \\qquad n = 0, 1, 2, \\dots $
%[text] with a damping factor $\\omega \> 0$ and $D = \\mathrm{diag}(A)$ for the Jacobi iteration (`Multigrid.Jacobi`) or $D = I$ for the Richardson iteration (`Multigrid.Richardson`). On these meshes $A$ is the five-point stencil: every diagonal entry is 4 and every node couples with $-1$ to its horizontal and vertical neighbours. Hence $D = 4I$, and the Jacobi iteration with the damping factor $\\omega$ performs exactly the same steps as the Richardson iteration with $\\omega/4$; the pairs used in the solver examples, Jacobi with $1/2$ and Richardson with $1/8$, coincide.
%[text] ### Error expansion in eigenvectors
%[text] The error $e\_n = x\_n - x$ satisfies $e\_{n+1} = S e\_n$ with the iteration matrix $S = I - \\omega D^{-1}A$, so $e\_n = S^n e\_0$ does not depend on $b$. Let $v\_1,\\dots,v\_M$ be the $D$-orthonormal eigenvectors of $Av\_j = \\lambda\_j Dv\_j$ with the eigenvalues $0 \< \\lambda\_1\\le\\dots\\le\\lambda\_M$. Then $Sv\_j = (1-\\omega\\lambda\_j)\\,v\_j$, and
%[text]{"align":"center"} $ e\_0 = \\sum\_{j=1}^M c\_j v\_j \\quad\\Longrightarrow\\quad e\_n = \\sum\_{j=1}^M (1-\\omega\\lambda\_j)^n\\,c\_j v\_j. $
%[text] Every component is multiplied by its own factor $1-\\omega\\lambda\_j$ in every step; the iteration converges if and only if all factors lie in $(-1,1)$, i.e. for $0 \< \\omega \< 2/\\lambda\_M$.
%[text] ### Smooth and oscillating eigenvectors
%[text] For the five-point stencil the eigenvectors are the grid functions $\\sin(k\\pi z\_1)\\sin(l\\pi z\_2)$ at the interior nodes $z$, $1\\le k,l\\le 1/h-1$, with the eigenvalues (for $D = 4I$)
%[text]{"align":"center"} $ \\lambda\_{kl} = 1 - \\frac{\\cos(k\\pi h) + \\cos(l\\pi h)}{2} \\in (0,2). $
%[text] Small frequencies $k$ and $l$ give smooth eigenvectors with small eigenvalues, down to $\\lambda\_{11}\\approx\\pi^2h^2/2$. The upper half of the spectrum, $\\lambda\_{kl}\\ge 1$, requires $\\cos(k\\pi h)\\le 0$ or $\\cos(l\\pi h)\\le 0$, i.e. $\\max(k,l)\\ge 1/(2h)$: these eigenvectors oscillate on the scale of the mesh and cannot be represented on the mesh with the double spacing $2h$.
%[text] ### Smoothing property
%[text] For $\\omega = 1/2$ all factors $1 - \\lambda\_{kl}/2$ lie in $(0,1)$, and those of the upper half of the spectrum are at most $1/2$: every step at least halves the oscillating part of the error, independently of $h$; this is the smoothing property. The smooth part is reduced only by the factor $1 - \\lambda\_{11}/2\\approx 1 - \\pi^2h^2/4$ per step, so the iteration needs of the order of $h^{-2}$ steps to reduce the whole error, in accordance with the condition number $\\lambda\_M/\\lambda\_1\\sim h^{-2}$ of $D^{-1}A$. The undamped Jacobi iteration $\\omega = 1$ is no smoother: the most oscillating eigenvector, the checkerboard pattern $k = l = 1/h - 1$, has the factor $1-\\lambda\_{kl} = -\\cos(\\pi h)\\approx -1$.
%[text] ### Expected behaviour
%[text] A random initial iterate contains all eigenvectors in its error. After a few steps the oscillating components are gone and the error is smooth, but not small: its norm decreases fast at first and then only very slowly. The frequency $k = l = 10$ of the exact solution lies in between; its factor $1 - (1-\\cos(10\\pi h))/2\\approx 0.78$ for $h = 2^{-5}$ removes it within a few dozen steps.
%%
%[text] ## Parameters
nrRefinements = 5;     % uniform refinements of the coarse mesh, grid spacing h = 2^-5
smoother = 'Jacobi';   % Jacobi iteration, D = diag(A)
damp = 1/2;            % damping factor omega
% damp = 1;            % undamped Jacobi: the checkerboard component is not damped
% smoother = 'Richardson'; damp = 1/8;  % D = I: the same steps as Jacobi with 1/2, since diag(A) = 4
nrIterations = 100;    % number of iterations
snapshots = [0 1 3 10 30 100];  % iterations at which the error is plotted
%%
%[text] ## Problem data and mesh
p = problemScalarEllipticSquareOscExact(10);
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform',nrRefinements);
plot(mesh);
%%
%[text] ## Linear system
%[text] The system of the interior nodes and its discrete solution $x$, the limit of the iterations. The diagonal entries of $A$ are all equal to 4. $D$ is the diagonal of $A$ for the Jacobi and the identity for the Richardson iteration.
fe = Pk(1,mesh);
free = fe.freeNodes;
a = BilinearForm('grad(u) * grad(v) dx');
A = a.assemble(fe);
A = A(free,free);
l = LinearForm('f(x) v dx');
b = l.assemble(fe,'f',@p.f);
b = b(free);
x = A\b;
diagonalEntries = unique(round(full(diag(A)),12))'
diagonal = struct('Jacobi',diag(A),'Richardson',ones(size(A,1),1));
D = spdiags(diagonal.(smoother),0,size(A,1),size(A,1));
%%
%[text] ## Smoothing iterations
%[text] The iteration starts from a random vector with entries in $(0,1)$; `rng(0)` makes it reproducible. The loop stores the norm of the error after every step and the error at the snapshot iterations.
rng(0);
xn = rand(numel(free),1);
e0 = xn - x;
error4iter = zeros(1,nrIterations+1);
error4iter(1) = norm(e0);
errors4snapshot = zeros(numel(free),numel(snapshots));
errors4snapshot(:,snapshots == 0) = e0;
for n = 1:nrIterations
    % smooth
    xn = xn + damp*(D\(b - A*xn));
    % errors
    error4iter(n+1) = norm(xn - x);
    if ismember(n,snapshots)
        errors4snapshot(:,snapshots == n) = xn - x;
    end
end
%%
%[text] ## Results
%[text] ### Error after n iterations
clf
tiledlayout('flow');
for k = 1:numel(snapshots)
    nexttile
    e = zeros(fe.nrDoF,1);
    e(free) = errors4snapshot(:,k);
    fe.plot(e,[],'EdgeColor','none');
    shading interp;
    axis equal tight
    colorbar;
    title("n = " + snapshots(k));
end
%[text] The noise of the initial error disappears in the first steps, the oscillation with the frequency of the exact solution within a few dozen steps. What remains is a smooth error that the iteration reduces only slowly.
%%
%[text] ### Oscillating and smooth parts of the error
%[text] `eig` computes the eigenvalues $\\lambda\_j$ and the $D$-orthonormal eigenvectors $v\_j$ of $Av = \\lambda Dv$. The coefficients $c = V^{-1}e$ of an error $e$ in the basis $V = (v\_1,\\dots,v\_M)$ split it into its oscillating part, the components of the upper half $\\lambda\_j\\ge\\lambda\_M/2$ of the spectrum, and its smooth part.
[V,lambda] = eig(full(A),full(D),'vector');
high = lambda >= max(lambda)/2;
c4snapshot = V\errors4snapshot;
table(snapshots',vecnorm(errors4snapshot)', ...
    vecnorm(V(:,high)*c4snapshot(high,:))',vecnorm(V(:,~high)*c4snapshot(~high,:))', ...
    VariableNames=["n" "error" "oscillating part" "smooth part"])
%[text] The oscillating part decreases at least by the factor $(1/2)^n$ until it reaches the level of the rounding errors; the smooth part decreases only very slowly.
%%
%[text] ### Coefficients over the eigenvalues
clf
semilogy(lambda,abs(c4snapshot),'.');
xlabel('\lambda_j');
ylabel('|c_j|');
legend("n = " + snapshots,'Location','southwest');
%[text] Each step multiplies the coefficient $c\_j$ by $1-\\omega\\lambda\_j$: the coefficients at the right end of the spectrum drop by many orders of magnitude, those at the left end stay almost unchanged.
%%
%[text] ### Error norms per iteration
%[text] The expansion $e\_n = \\sum\_j (1-\\omega\\lambda\_j)^n c\_j v\_j$ with the coefficients $c\_j$ of $e\_0$ gives the two parts for all iterations in exact arithmetic; the total error is the one of the iteration above.
c0 = V\e0;
nAll = 0:nrIterations;
oscillating4iter = vecnorm(V(:,high)*(c0(high).*(1 - damp*lambda(high)).^nAll));
smooth4iter = vecnorm(V(:,~high)*(c0(~high).*(1 - damp*lambda(~high)).^nAll));
clf
semilogy(nAll,error4iter,'-',nAll,oscillating4iter,'--',nAll,smooth4iter,':','LineWidth',1.5);
xlabel('iteration n');
ylim([1e-16,1e2]);
legend('error','oscillating part','smooth part');
%[text] After the first steps the error equals its smooth part, whose slowest component the iteration reduces only by the factor $1-\\omega\\lambda\_1$ close to one per step. A correction on a coarser mesh, which represents the smooth error well, removes it (`exSolverTwoGrid`).
%%
%[text] ## Possibilities for extensions
%[text] - Set `damp = 1`, the undamped Jacobi iteration: the checkerboard component has the factor $-\\cos(\\pi h)\\approx -1$, so it changes its sign in every step and the oscillating part decreases only slowly.
%[text] - Set `damp = 2/3`: the factors of the upper half of the spectrum are at most $1/3$ instead of $1/2$.
%[text] - Set `smoother = 'Richardson'; damp = 1/8;`: the iterates are the same up to rounding errors, because $D = 4I$; only the eigenvalues in the plot are those of $A$, four times larger. The two iterations differ where the diagonal entries vary, e.g. for the diffusion coefficient $\\kappa = 1$ for $x\_1 \< 1/2$ and $\\kappa = 100$ otherwise, `A = BilinearForm('kappa(x)*grad(u)*grad(v) dx').assemble(fe,'kappa',@(x)1+99*(x(:,1)>1/2));`. Richardson then needs $\\omega \< 2/\\lambda\_M(A)$, e.g. `damp = 1/eigs(A,1);`, and hardly smooths the half with $\\kappa = 1$, while Jacobi with $\\omega = 1/2$ smooths both halves.
%[text] - Set `nrRefinements = 6`, the limit for the dense eigenvalue decomposition: the oscillating part decreases with the same factor, while the smooth part decreases even more slowly, since $\\lambda\_1\\sim h^2$. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
