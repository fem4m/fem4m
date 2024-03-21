%[text] # Multigrid solver for the Poisson problem
%[text] Solves the linear systems of the $P\_1$ discretization of the Poisson problem with a geometric multigrid V-cycle on a hierarchy of uniformly refined meshes. The number of cycles does not grow with the number of unknowns $N$, so the cost of a solve is proportional to $N$.
%[text] **Builds on:** `exPoissonP1`, `exSolverSmoothing`, `exSolverTwoGrid`. **Next:** `exSolverSchwarz`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] On the unit square $\\Omega = (0,1)^2$ find $u$ with $-\\Delta u = f$ in $\\Omega$ and $u = 0$ on $\\partial\\Omega$. The problem object `problemScalarEllipticSquareOscExact` provides the source of the exact solution $u = \\sin(\\pi x\_1)\\sin(\\pi x\_2)$. Uniform refinement of a coarse mesh gives the meshes $\\mathcal{T}\_1, \\mathcal{T}\_2, \\dots, \\mathcal{T}\_L$ with nested $P\_1$ spaces $V\_1\\subset V\_2\\subset\\dots\\subset V\_L$ and the linear systems $A\_\\ell x\_\\ell = b\_\\ell$ of `exPoissonP1` on every level $\\ell$.
%[text] ### Why multigrid
%[text] The condition number of $A\_\\ell$ grows like $h\_\\ell^{-2}$, so classical iterations need more and more steps on finer meshes: Jacobi and Richardson iterations about $h\_\\ell^{-2}$, the conjugate gradient method about $h\_\\ell^{-1}$. A few Jacobi or Richardson steps, however, reduce the oscillating components of the error fast (`exSolverSmoothing`). The remaining smooth error can be represented on the coarser mesh, where it is cheaper to compute (`exSolverTwoGrid`). Multigrid applies this idea recursively on all levels.
%[text] ### The V-cycle
%[text] The prolongation $P\_\\ell$ maps the coefficients of a function of $V\_{\\ell-1}$ to its coefficients in $V\_\\ell$ (`fe.prolongation(feOld)`), the restriction is its transpose $P\_\\ell^T$. As $V\_{\\ell-1}\\subset V\_\\ell$, the matrix assembled on the coarser mesh is the Galerkin matrix $A\_{\\ell-1} = P\_\\ell^T A\_\\ell P\_\\ell$. One V-cycle on level $\\ell$ for the iterate $x$ and the right-hand side $b$
%[text] 1. smooths with $m$ steps $x \\leftarrow x + \\omega D\_\\ell^{-1}(b - A\_\\ell x)$, the damped Jacobi method with the diagonal $D\_\\ell$ of $A\_\\ell$ ($\\omega = 1/2$) or the Richardson method with $D\_\\ell = I$ ($\\omega = 1/8$),
%[text] 2. restricts the residual, $d = P\_\\ell^T(b - A\_\\ell x)$, computes a correction $e$ with $A\_{\\ell-1}e \\approx d$ by one V-cycle on level $\\ell-1$ from $e = 0$ (on the coarsest level by an LU factorization), and updates $x \\leftarrow x + P\_\\ell e$,
%[text] 3. smooths again with $m$ steps. \
%[text] `Multigrid(smoother)` builds the solver with $m = 2$; `addLevel` adds the matrices $A\_\\ell$ and prolongations $P\_\\ell$ from the coarsest to the finest level, restricted to the free nodes, and `solve` runs V-cycles on the finest level. On the uniform meshes of this example every diagonal entry of $A\_\\ell$ equals 4, so the two smoothers perform the same steps.
%[text] ### Expected convergence
%[text] The V-cycle is a contraction with a contraction number $\\rho \< 1$ that does not depend on the level. The number of cycles that reduce the residual by a fixed factor is therefore bounded independently of $N$. One cycle costs $O(N)$ operations, since the numbers of unknowns on the coarser levels decrease geometrically, so the cost of a solve is $O(N)$.
%%
%[text] ## Parameters
nrLevels = 8;         % number of levels
smoother = 'Jacobi';  % damped Jacobi smoother with omega = 1/2
% smoother = 'Richardson';  % Richardson smoother with omega = 1/8, the same steps on these meshes
relTol = 1e-8;        % stop when the residual norm is reduced by this factor
%%
%[text] ## Problem data and initial mesh
p = problemScalarEllipticSquareOscExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Multigrid on a hierarchy of meshes
%[text] Every level adds its matrix and the prolongation from the previous level to the solver and solves its system with V-cycles from the zero initial iterate.
mg = Multigrid(smoother);
nrDoF4lvl = zeros(1,nrLevels);
cycles4lvl = zeros(1,nrLevels);
reduction4lvl = zeros(1,nrLevels);
H1Error4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element space
    fe = Pk(1,mesh);
    nrDoF4lvl(lvl) = fe.nrDoF;
    free = fe.freeNodes;
    % assemble
    a = BilinearForm('grad(u) * grad(v) dx');
    A = a.assemble(fe);
    l = LinearForm('f(x) v dx');
    b = l.assemble(fe,'f',@p.f);
    % solve: add the level, then run V-cycles
    if lvl == 1
        mg = mg.addLevel(A(free,free));
    else
        P = fe.prolongation(feOld);
        mg = mg.addLevel(A(free,free),P(free,feOld.freeNodes));
    end
    x = zeros(fe.nrDoF,1);
    [x(free),cycles4lvl(lvl)] = mg.solve(b(free),relTol*norm(b(free)),100,[],0);
    residual = norm(b(free) - A(free,free)*x(free))/norm(b(free));
    reduction4lvl(lvl) = residual^(1/cycles4lvl(lvl));
    % errors
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(lvl) = h1.assemble(x,fe,'exact',@p.exactGradU);
    % refine
    feOld = fe;
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Cycles per level
%[text] The reduction $(\\|r\_n\\|/\\|r\_0\\|)^{1/n}$ is the average factor by which one of the $n$ cycles reduces the norm of the residual $r = b - Ax$.
table(nrDoF4lvl',cycles4lvl',reduction4lvl',H1Error4lvl', ...
    VariableNames=["N" "V-cycles" "reduction" "H1 error"])
%[text] While $N$ grows by a factor of about four from level to level, the number of cycles and the reduction factor stay bounded: the solver is optimal.
%%
%[text] ### Discrete solution
clf
plot(fe,x);
view([40,35]);
colorbar;
shading interp;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(u-u_h)||_{L^2}');
%[text] The error of the multigrid solution decays with the rate $N^{-1/2}$ of the discretization error: the algebraic error after the V-cycles is negligible.
%%
%[text] ## Possibilities for extensions
%[text] - Change the damping factor, e.g. `Multigrid(smoother,0.8)`: stronger damping of the oscillating error components needs fewer cycles, and their number still does not grow with $N$.
%[text] - Use a W-cycle with two coarse grid corrections per level, `Multigrid(smoother,[],[],[],2)`, or one pre- and one post-smoothing step, `Multigrid(smoother,[],1)`.
%[text] - Solve with `pcg(A(free,free),b(free),relTol,1000)` instead and compare the numbers of iterations, which grow like $h^{-1}$.
%[text] - `mg.solve(...,2)` also estimates the contraction number of the V-cycle.
%[text] - `exSolverSchwarz` uses multilevel and domain decomposition preconditioners for the conjugate gradient method. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
