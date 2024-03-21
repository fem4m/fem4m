%[text] # Additive and multiplicative Schwarz preconditioners
%[text] Preconditions the conjugate gradient method for $P\_k$ finite elements of the Poisson problem with multilevel Schwarz methods, which combine local solves on all levels of a hierarchy of meshes, either additively or multiplicatively. With the BPX and the domain decomposition preconditioners the number of iterations stays bounded under refinement, while it grows like $h^{-1}$ without a preconditioner.
%[text] **Builds on:** `exSolverMultigrid`. **Next:** none, the last solver example.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] On the unit square $\\Omega = (0,1)^2$ find $u$ with $-\\Delta u = f$ in $\\Omega$ and $u = 0$ on $\\partial\\Omega$. The problem object `problemScalarEllipticSquareExact` provides the source $f = 2x\_1(1-x\_1) + 2x\_2(1-x\_2)$ of the exact solution $u = x\_1(1-x\_1)\\,x\_2(1-x\_2)$. Uniform refinement gives the meshes $\\mathcal{T}\_1,\\dots,\\mathcal{T}\_L$ with nested spaces $V\_1\\subset\\dots\\subset V\_L$ of continuous piecewise polynomials of degree $k$ (`Pk(degree,mesh)`) and, on every level $\\ell$, the symmetric positive definite system $A\_\\ell x = b\_\\ell$ for the free degrees of freedom.
%[text] ### Preconditioned conjugate gradients
%[text] The conjugate gradient method reduces the energy norm of the error by the factor $2\\left((\\sqrt\\kappa-1)/(\\sqrt\\kappa+1)\\right)^n$ in $n$ iterations, where $\\kappa$ is the condition number of the system matrix. As $\\kappa(A\_L)\\sim h^{-2}$, the number of iterations grows like $h^{-1}\\sim N^{1/2}$. The preconditioned method (`pcg(A,b,tol,maxit,B)`) applies a symmetric positive definite $B\\approx A\_L^{-1}$ to the residual in every iteration and converges with $\\kappa(BA\_L)$ instead; a preconditioner is optimal if $\\kappa(BA\_L)$ does not depend on $h$ and $B$ costs $O(N)$ operations.
%[text] ### Subspace correction
%[text] Schwarz methods split $V\_L$ into subspaces of all levels and solve small problems in them. In matrix form, level $\\ell$ provides a local solve operator $B\_\\ell\\approx A\_\\ell^{-1}$, a damped inverse of the diagonal or a sum of inverses of small submatrices, and the prolongation $I\_\\ell$ from level $\\ell-1$ (`fe.prolongation(feOld)`); $P\_\\ell = I\_L I\_{L-1}\\cdots I\_{\\ell+1}$ maps level $\\ell$ to level $L$. `addLevel` stores $A\_\\ell$, $I\_\\ell$ and $B\_\\ell$ from the coarsest to the finest level.
%[text] - The additive Schwarz preconditioner (`applyAdd`) corrects on all levels independently, with the same residual $r$: $B\_{\\textrm{add}}\\,r = \\sum\_{\\ell=1}^{L} P\_\\ell B\_\\ell P\_\\ell^T r$.
%[text] - The multiplicative Schwarz preconditioner (`applyMult`) corrects level by level with the updated residual. On level $\\ell$, from $x = B\_\\ell r$, it adds the correction of the levels below for the residual $I\_\\ell^T(r - A\_\\ell x)$, transferred with $I\_\\ell$, and then $B\_\\ell(r - A\_\\ell x)$ once more. This is a V-cycle with one pre- and one post-smoothing step with $B\_\\ell$ (`exSolverMultigrid`), without an exact solve on the coarsest level. The symmetric post-smoothing step makes the preconditioner symmetric; it is positive definite if $\\lambda\_{\\max}(B\_\\ell A\_\\ell) \< 2$ on every level. \
%[text] The additive variant can treat the levels in parallel, the multiplicative variant usually needs fewer iterations.
%[text] ### The preconditioners
%[text] The classes differ in the local solve operators $B\_\\ell$, damped with the factor $\\omega$ of their constructor. The domain decomposition preconditioners use the elements as subdomains on every level.
%[text:table]
%[text] | Class | Local solve $B\_\\ell$ on level $\\ell$ | $\\kappa(B\_{\\textrm{add}}A\_L)$ |
%[text] | --- | --- | --- |
%[text] | `SchwarzHB` | hierarchical basis: $\\omega I$ on the free nodes of level $\\ell$ that are not nodes of level $\\ell-1$ (all free nodes on level 1), $\\omega = 1/8$ | $O(\\vert\\log h\\vert^2)$ |
%[text] | `SchwarzBPX` | Bramble-Pasciak-Xu: $\\omega I$ on all free degrees of freedom, $\\omega = 1/8$ | $O(1)$ |
%[text] | `SchwarzVO` | overlapping vertex patches: the inverses of $A\_\\ell$ on the degrees of freedom of every vertex patch without those on its boundary, scaled with the inverse number of patches that share a degree of freedom, $\\omega = 1/2$ | $O(1)$ |
%[text] | `SchwarzBPS` | Bramble-Pasciak-Schatz: $\\mathrm{diag}(A\_\\ell)^{-1}$ on the vertex degrees of freedom plus the inverses of $A\_\\ell$ on the interior degrees of freedom of every edge and every element, $\\omega = 1/2$ | $O(1)$ |
%[text] | `SchwarzNN` | Neumann-Neumann: $\\mathrm{diag}(A\_\\ell)^{-1}$ on the vertices, the inverses on the interior degrees of freedom of every element, and in addition the inverses of the element matrices (Neumann problems, regularized with the element mass matrix divided by the area) on all other degrees of freedom of every element, the vertices included, scaled with the inverse number of elements that share them, $\\omega = 1/2$ | $O(1)$ |
%[text] | `SchwarzBDDC` | balancing domain decomposition by constraints: $\\mathrm{diag}(A\_\\ell)^{-1}$ on the primal degrees of freedom (the vertices), the inverses on the interior degrees of freedom of every element, and the element problems on the remaining degrees of freedom with zero primal values, scaled with the inverse number of elements that share them, $\\omega = 1/2$ | $O(1)$ |
%[text:table]
%[text] The bound for HB holds in two space dimensions; the uniform bound for BPX rests on the stability of multilevel decompositions in $H^1$. The local solves of VO, BPS, NN and BDDC handle the edge and interior degrees of freedom of the polynomial degree $k$ within the elements and their patches, and the hierarchy of levels handles the smooth error components, so their condition numbers are bounded independently of $h$, with constants that depend on $k$. For $k = 1$ there are no edge or interior degrees of freedom: the local solves of VO, BPS and BDDC reduce to $\\omega\\,\\mathrm{diag}(A\_\\ell)^{-1}$, which coincides with BPX on the uniform meshes of the square, where every diagonal entry is 4. The multiplicative variants have the same asymptotic behaviour with smaller constants.
%[text] ### Expected convergence
%[text] Without a preconditioner the number of iterations doubles from level to level, since $h$ is halved. With HB it grows linearly in the number of levels, $|\\log h|$; with the other preconditioners it stays bounded. The discrete solution converges to $u$ with the rate $N^{-k/2}$ in the energy norm.
%%
%[text] ## Parameters
degree = 2;                     % polynomial degree k
nrLevels = 7;                   % number of levels
preconditioner = SchwarzBDDC;   % or SchwarzHB, SchwarzBPX, SchwarzVO, SchwarzBPS, SchwarzNN
schwarzType = 'Mult';           % 'Mult' (multiplicative, applyMult) or 'Add' (additive, applyAdd)
relTol = 1e-6;                  % stop when the residual norm is reduced by this factor
maxIter = 2000;                 % maximal number of iterations
%%
%[text] ## Problem data and initial mesh
%[text] One uniform refinement of the coarse mesh of the square gives the first level.
p = problemScalarEllipticSquareExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Preconditioned conjugate gradients on a hierarchy of meshes
%[text] Every level adds its space, matrix, element matrices `A4e` (needed by NN and BDDC) and prolongation to the preconditioner and solves its system twice from the zero initial iterate: with the plain and with the preconditioned conjugate gradient method. The outputs `flag` (0 if the tolerance was reached) and `iter` of `pcg` replace its messages.
a = BilinearForm('grad(u) * grad(v) dx');
l = LinearForm('f(x) v dx');
nrDoF4lvl = zeros(1,nrLevels);
cgIter4lvl = zeros(1,nrLevels);
pcgIter4lvl = zeros(1,nrLevels);
pcgFlag4lvl = zeros(1,nrLevels);
H1Error4lvl = zeros(1,nrLevels);
feOld = [];
for lvl = 1:nrLevels
    % finite element space
    fe = Pk(degree,mesh);
    nrDoF4lvl(lvl) = fe.nrDoF;
    free = fe.freeNodes;
    % assemble
    [A,A4e] = a.assemble(fe);
    b = l.assemble(fe,'f',@p.f);
    % solve: add the level, then run CG and PCG
    if lvl == 1
        P = [];
    else
        P = fe.prolongation(feOld);
    end
    preconditioner = preconditioner.addLevel(fe,feOld,P,A,A4e);
    B = @(r)preconditioner.("apply" + schwarzType)(r);
    [~,~,~,cgIter4lvl(lvl)] = pcg(A(free,free),b(free),relTol,maxIter);
    x = zeros(fe.nrDoF,1);
    [x(free),pcgFlag4lvl(lvl),~,pcgIter4lvl(lvl)] = pcg(A(free,free),b(free),relTol,maxIter,B);
    % errors
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(lvl) = h1.assemble(x,fe,'exact',@p.exactGradU);
    % refine
    feOld = fe;
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Iterations per level
table(nrDoF4lvl',cgIter4lvl',pcgIter4lvl',pcgFlag4lvl',H1Error4lvl', ...
    VariableNames=["N" "CG iterations" "PCG iterations" "PCG flag" "H1 error"])
%[text] The number of CG iterations roughly doubles from level to level, while the number of PCG iterations stays bounded: the multilevel preconditioner removes the dependence on the mesh size.
%%
%[text] ### Iterations over the number of unknowns
clf
loglog(nrDoF4lvl,cgIter4lvl,'-o',nrDoF4lvl,pcgIter4lvl,'-s','LineWidth',1.5);
xlabel('N');
ylabel('iterations');
legend('CG','PCG','Location','northwest');
%[text] In the log-log scale the CG iterations approach a line with the slope $1/2$, the PCG iterations a horizontal line.
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
%[text] The error decays with the rate $N^{-k/2}$ of the discretization error; the algebraic error of PCG is negligible.
%%
%[text] ## Possibilities for extensions
%[text] - Compare the preconditioners and both types: with `schwarzType = 'Add'` the iteration counts are larger and grow slowly, with decreasing increments; `preconditioner = SchwarzHB;` shows the growth like $|\\log h|$ of the hierarchical basis preconditioner.
%[text] - Set `degree = 1`: VO, BPS and BDDC give exactly the iteration counts of BPX, NN differs. For `degree = 3` the domain decomposition preconditioners need fewer iterations than BPX, whose damped identity treats all degrees of freedom alike.
%[text] - The multiplicative variant needs $\\lambda\_{\\max}(B\_\\ell A\_\\ell) \< 2$ on every level. For `SchwarzNN` with `degree = 2` the default damping factor $1/2$ violates it: the iteration counts grow rapidly and PCG does not reach the tolerance on the finest level (flag 1). `SchwarzNN(1/4)` restores bounded iteration counts.
%[text] - `SchwarzBPS` also works for `Pk3d` on tetrahedral meshes (`Mesh3d`), where it solves on the interior degrees of freedom of the faces as well, and `SchwarzVO` for discontinuous elements such as `DGDk` with the interior penalty matrix of `exPoissonDG`. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
