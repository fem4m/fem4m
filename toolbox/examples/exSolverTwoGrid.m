%[text] # Two-grid method for the Poisson problem
%[text] Solves the linear system of the $P\_1$ discretization of the Poisson problem with the two-grid method: a few Richardson steps smooth the error on the fine mesh, and an exact solve on the next coarser mesh removes the remaining smooth error. The error contracts in every step by a factor that does not depend on the mesh size, so the number of steps stays bounded under refinement.
%[text] **Builds on:** `exPoissonP1`, `exSolverSmoothing`. **Next:** `exSolverMultigrid`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] On the unit square $\\Omega = (0,1)^2$ find $u$ with $-\\Delta u = f$ in $\\Omega$ and $u = 0$ on $\\partial\\Omega$. The problem object `problemScalarEllipticSquareExact` provides the source $f = 2x\_1(1-x\_1) + 2x\_2(1-x\_2)$ of the exact solution $u = x\_1(1-x\_1)\\,x\_2(1-x\_2)$. The uniform refinement $\\mathcal{T}\_h$ of a mesh $\\mathcal{T}\_H$ gives nested $P\_1$ spaces $V\_H\\subset V\_h$ and the linear systems $A x = b$ on the fine and $A\_H x\_H = b\_H$ on the coarse mesh, both for the coefficients of the interior nodes, as in `exPoissonP1`.
%[text] ### Smoothing and coarse grid correction
%[text] The Richardson iteration $x \\leftarrow x + \\omega(b - Ax)$ quickly removes the oscillating components of the error, but hardly changes the smooth ones (`exSolverSmoothing`). A smooth error is well approximated on the coarse mesh. The prolongation $P$ (`fe.prolongation(feCoarse)`) maps the coefficients of a function of $V\_H$ to its coefficients in $V\_h$, and its transpose $P^T$ restricts residuals to the coarse mesh. As $V\_H\\subset V\_h$, the matrix assembled on the coarse mesh is the Galerkin matrix $A\_H = P^T A P$. One step of the two-grid method from the iterate $x$
%[text] 1. smooths with $m$ steps $x \\leftarrow x + \\omega(b - Ax)$ of the Richardson iteration with $\\omega = 1/8$ (`Multigrid.Richardson`),
%[text] 2. restricts the residual, $d = P^T(b - Ax)$, solves the coarse system $A\_H e\_H = d$ exactly with an LU factorization of $A\_H$, which is computed once, and corrects $x \\leftarrow x + P e\_H$. \
%[text] On the uniform meshes of this example every diagonal entry of $A$ is 4, so the Jacobi smoother with $\\omega = 1/2$ (`Multigrid.Jacobi`) performs exactly the same steps as the Richardson smoother with $\\omega = 1/8$.
%[text] ### Two-grid error operator
%[text] The error $e = x - A^{-1}b$ changes in one step by the smoothing operator $S = I - \\omega A$ and the coarse grid correction $K = I - PA\_H^{-1}P^TA$, so the new error is $Ee$ with the two-grid operator
%[text]{"align":"center"} $ E = (I - PA\_H^{-1}P^TA)\\,S^m = (A^{-1} - PA\_H^{-1}P^T)\\,AS^m. $
%[text] $K$ is the $A$-orthogonal projection onto the complement of the range of $P$: it removes the error components of $V\_H$ and never increases the energy norm $\\|e\\|\_A = (e^TAe)^{1/2}$. Hackbusch's convergence theory bounds the two factors of the second form separately:
%[text] - Smoothing property: for $\\omega\\lambda\_{\\max}(A)\\le 1$, $\\|AS^m\\|\_2 = \\max\_j \\lambda\_j(1-\\omega\\lambda\_j)^m \\le \\eta(m)/\\omega$ with $\\eta(m) = m^m/(m+1)^{m+1} \< 1/(em)$.
%[text] - Approximation property: $\\|A^{-1} - PA\_H^{-1}P^T\\|\_2\\le C\_A/\\lambda\_{\\max}(A)$, a consequence of the $H^2$ regularity of the Poisson problem on the convex square and of the $L^2$ error estimate of the coarse Galerkin solution, with a constant $C\_A$ that does not depend on $h$. \
%[text] Together, $\\|E\\|\_2\\le C\_A\\,\\eta(m)/(\\omega\\lambda\_{\\max}(A))$. Here $\\lambda\_{\\max}(A) \< 8$ approaches 8 under refinement, so $\\omega = 1/8$ satisfies the condition, and the bound does not depend on $h$ and is smaller than one for sufficiently many smoothing steps $m$.
%[text] ### Expected convergence
%[text] The residual $r = b - Ax$ contracts in every step by about the same factor $\\rho \< 1$ on all meshes, so the number of steps that reduce the residual by the factor `relTol` stays bounded as $N$ grows. One step costs $O(N)$ operations plus the coarse solve. The coarse system, however, has about $N/4$ unknowns, so its LU factorization dominates the cost on fine meshes; the multigrid method of `exSolverMultigrid` replaces it by the same method on the next coarser level, recursively. The discrete solution itself converges to $u$ with the rate $N^{-1/2}$ in the energy norm.
%%
%[text] ## Parameters
nrLevels = 9;             % number of two-grid solves, on the meshes refined 2,...,nrLevels+1 times
nrSmoothingSteps = 5;     % number m of smoothing steps
smoother = 'Richardson';  % Richardson smoother x - omega*(A*x - b)
damp = 1/8;               % damping factor omega
% smoother = 'Jacobi'; damp = 1/2;  % Jacobi smoother: the same steps on these meshes, since diag(A) = 4
relTol = 1e-6;            % stop when the residual norm is reduced by this factor
%%
%[text] ## Problem data and initial mesh
%[text] The coarse mesh of the square has two triangles; one uniform refinement gives the first coarse mesh with an interior node.
p = problemScalarEllipticSquareExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Two-grid method on a sequence of meshes
%[text] Every level refines the current mesh, the coarse mesh $\\mathcal{T}\_H$, once to the fine mesh $\\mathcal{T}\_h$ and solves the fine system with two-grid steps from the zero initial iterate until the residual norm is reduced by `relTol`. The loop stores the residual norms of all steps.
a = BilinearForm('grad(u) * grad(v) dx');
l = LinearForm('f(x) v dx');
smooth = str2func("Multigrid." + smoother);
nrDoF4lvl = zeros(1,nrLevels);
steps4lvl = zeros(1,nrLevels);
contraction4lvl = zeros(1,nrLevels);
H1Error4lvl = zeros(1,nrLevels);
residual4lvl = cell(1,nrLevels);
for lvl = 1:nrLevels
    % finite element spaces on the coarse and the fine mesh
    feCoarse = Pk(1,mesh);
    freeCoarse = feCoarse.freeNodes;
    mesh = mesh.refine('uniform');
    fe = Pk(1,mesh);
    free = fe.freeNodes;
    nrDoF4lvl(lvl) = fe.nrDoF;
    % assemble
    A = a.assemble(fe);
    A = A(free,free);
    b = l.assemble(fe,'f',@p.f);
    b = b(free);
    AH = a.assemble(feCoarse);
    [L,U,rowPerm,colPerm] = lu(AH(freeCoarse,freeCoarse));
    P = fe.prolongation(feCoarse);
    P = P(free,freeCoarse);
    % solve with two-grid steps
    xFree = zeros(numel(free),1);
    residual4lvl{lvl} = 1;
    while residual4lvl{lvl}(end) > relTol
        for step = 1:nrSmoothingSteps
            xFree = smooth(xFree,b,A,damp);
        end
        d = P'*(b - A*xFree);
        xFree = xFree + P*(colPerm*(U\(L\(rowPerm*d))));
        residual4lvl{lvl}(end+1) = norm(b - A*xFree)/norm(b);
    end
    steps4lvl(lvl) = numel(residual4lvl{lvl}) - 1;
    contraction4lvl(lvl) = residual4lvl{lvl}(end)^(1/steps4lvl(lvl));
    % errors
    x = zeros(fe.nrDoF,1);
    x(free) = xFree;
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(lvl) = h1.assemble(x,fe,'exact',@p.exactGradU);
end
%%
%[text] ## Results
%[text] ### Steps per level
%[text] The contraction $(\\|r\_n\\|/\\|r\_0\\|)^{1/n}$ is the average factor by which one of the $n$ two-grid steps reduces the norm of the residual.
table(nrDoF4lvl',steps4lvl',contraction4lvl',H1Error4lvl', ...
    VariableNames=["N" "two-grid steps" "contraction" "H1 error"])
%[text] While $N$ grows by a factor of about four from level to level, the number of steps and the contraction factor stay bounded.
%%
%[text] ### Residual histories
clf
for lvl = 1:nrLevels
    semilogy(0:steps4lvl(lvl),residual4lvl{lvl},'-o');
    hold on
end
hold off
xlabel('two-grid step n');
ylabel('||r_n|| / ||r_0||');
legend("N = " + nrDoF4lvl,'Location','northeast');
%[text] After the first step, in which the residual norm may even grow slightly, the residual decreases geometrically with almost the same factor on all meshes.
%%
%[text] ### Discrete solution
clf
plot(fe,x,[],'EdgeColor','none');
view([40,35]);
colorbar;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(u-u_h)||_{L^2}');
%[text] The error of the two-grid solution decays with the rate $N^{-1/2}$ of the discretization error: the algebraic error after the two-grid steps is negligible.
%%
%[text] ## Possibilities for extensions
%[text] - Vary `nrSmoothingSteps`: more smoothing steps reduce the contraction factor, as the bound $C\_A\\,\\eta(m)/(\\omega\\lambda\_{\\max})$ predicts, but every step becomes more expensive.
%[text] - Add post-smoothing: repeat the smoothing loop after the coarse grid correction. With the same number of pre- and post-smoothing steps the two-grid operator is symmetric in the energy inner product, and the method can serve as a preconditioner for the conjugate gradient method.
%[text] - Drop the coarse grid correction, the two lines with `d` and `P`, and set `nrLevels = 4`: the method becomes the plain Richardson iteration of `exSolverSmoothing`, and the number of steps grows like $h^{-2}$.
%[text] - `exSolverMultigrid` solves the coarse system approximately by the same method on the next coarser mesh, recursively, which reduces the cost of a step to $O(N)$. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
