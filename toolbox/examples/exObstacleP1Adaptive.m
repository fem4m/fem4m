%[text] # Adaptive $P\_1$ finite elements for the obstacle problem
%[text] Solves the obstacle problem, a variational inequality whose solution touches a given obstacle on an unknown contact set, with continuous piecewise linear finite elements on adaptively refined meshes. New are the constraints at the nodes, the primal-dual active set method `pdaSet` for the resulting quadratic program and the discrete Lagrange multiplier.
%[text] **Builds on:** `exPoissonPkAdaptive`. **Next:** `exEllipticCRAdaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] An elastic membrane over $\\Omega$, clamped at the boundary and loaded by the force $f$, must stay above the obstacle $\\chi$. Its deflection $u$ satisfies the complementarity conditions
%[text]{"align":"center"} $ -\\Delta u \\ge f, \\qquad u \\ge \\chi, \\qquad (-\\Delta u - f)(u - \\chi) = 0 \\textrm{ in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\partial\\Omega. $
%[text] On the contact set $\\{u = \\chi\\}$ the membrane lies on the obstacle, elsewhere it solves the Poisson equation $-\\Delta u = f$; the free boundary between both parts is unknown. The Lagrange multiplier $\\lambda = -\\Delta u - f\\ge 0$ is the contact force, it vanishes outside the contact set.
%[text] The problem object `problemObstacleSquareExact` takes the square $\\Omega = (-1,1)^2$, the obstacle $\\chi = 0$ (`p.lower`; `p.upper` returns $\\infty$, no upper obstacle) and, with $r = |x|$ and $r\_0 = 0.7$, the exact solution
%[text]{"align":"center"} $ u = (r^2-r\_0^2)^2 \\textrm{ for } r \> r\_0, \\qquad u = 0 \\textrm{ for } r\\le r\_0, $
%[text] with the boundary data $u\_D = u$ on $\\partial\\Omega$. The source is $f = -\\Delta u = -8(2r^2-r\_0^2)$ for $r \> r\_0$ and $f = -8r\_0^2(1+r\_0^2-r^2) \< 0$ on the contact set $\\{r\\le r\_0\\}$, where it gives the multiplier $\\lambda = -f$. The solution is continuously differentiable, its second derivatives jump at the free boundary $r = r\_0$, so $u\\in H^2(\\Omega)$, but $u\\notin H^3(\\Omega)$.
%[text] ### Weak formulation
%[text] With the convex set $K = \\{v\\in H^1(\\Omega) : v = u\_D \\textrm{ on } \\partial\\Omega,\\ v\\ge\\chi \\textrm{ in } \\Omega\\}$ the problem is the variational inequality: find $u\\in K$ with
%[text]{"align":"center"} $ \\int\_\\Omega \\nabla u\\cdot\\nabla(v-u)\\,dx \\ge \\int\_\\Omega f(v-u)\\,dx \\quad \\textrm{for all } v\\in K. $
%[text] Equivalently, $u$ minimizes the energy $E(v) = \\frac{1}{2}\\int\_\\Omega|\\nabla v|^2\\,dx - \\int\_\\Omega f v\\,dx$ over $K$; the minimizer exists and is unique.
%[text] ### Discretization
%[text] The discrete set $K\_h$ contains the continuous piecewise linear functions $v\_h$ (`Pk(1,mesh)`) with $v\_h = u\_D$ at the boundary nodes and $v\_h(z)\\ge\\chi(z)$ at all free nodes $z$. A piecewise linear function that lies above the obstacle $\\chi = 0$ at the vertices lies above it everywhere, so every $v\_h\\in K\_h$ satisfies the constraint $v\_h\\ge\\chi$ in $\\Omega$; only its boundary values interpolate $u\_D$. For Lagrange elements of degree $k\\ge 2$ the constraints at the nodes do not imply $v\_h\\ge\\chi$ between the nodes, which is why this example uses $P\_1$. With the stiffness matrix $A$ and the load vector $b$ restricted to the free nodes, after subtracting the boundary values as in `exPoissonPkAdaptive`, the discrete problem is the quadratic program
%[text]{"align":"center"} $ \\min\_x\\ \\frac{1}{2}x^T A x - b^T x \\quad \\textrm{subject to } x\_j\\ge\\chi\_j = \\chi(z\_j) \\textrm{ for all free nodes } z\_j. $
%[text] Its Karush-Kuhn-Tucker conditions with the discrete multiplier $\\lambda = A x - b$ are
%[text]{"align":"center"} $ \\lambda\\ge 0, \\qquad x\\ge\\chi, \\qquad \\lambda\_j(x\_j-\\chi\_j) = 0 \\textrm{ for all } j, $
%[text] the discrete counterpart of the complementarity conditions; $\\lambda\_j = \\int\_\\Omega \\nabla u\_h\\cdot\\nabla\\varphi\_j\\,dx - \\int\_\\Omega f\\varphi\_j\\,dx$ is the residual of the Poisson equation at the node $z\_j$.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx$ | `grad(u) * grad(v) dx` |
%[text] | $\\int\_\\Omega f v\\,dx$ | `f(x) v dx` |
%[text:table]
%[text] ### Primal-dual active set method
%[text] `[x,lambda,k] = pdaSet(A,b,low,up)` solves the quadratic program with the bounds `low` $\\le x\\le$ `up` by the primal-dual active set method. Given the iterates $x$ and $\\lambda$, it predicts the active set $\\mathcal{A} = \\{j : \\lambda\_j + c(\\chi\_j - x\_j) \> 0\\}$ of the nodes in contact, with the constant $c = 10^3$, and the inactive set $\\mathcal{I}$ of the others, then solves
%[text]{"align":"center"} $ x\_j = \\chi\_j \\textrm{ for } j\\in\\mathcal{A}, \\qquad (A x - b)\_j = 0 \\textrm{ for } j\\in\\mathcal{I}, $
%[text] a linear system for $x\_{\\mathcal{I}}$, and sets $\\lambda = A x - b$ on $\\mathcal{A}$ and $\\lambda = 0$ on $\\mathcal{I}$. It stops when the active set does not change any more; then $x$ and $\\lambda$ satisfy the Karush-Kuhn-Tucker conditions exactly. The method is a semismooth Newton method for the complementarity conditions; for an M-matrix $A$, which the $P\_1$ stiffness matrices on the meshes of this example are (no positive off-diagonal entries), it stops after finitely many steps. The third output $k$ is their number. The iteration starts from $x = 0$.
%[text] ### Error estimator
%[text] The example uses the residual estimator of the Poisson problem,
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{T\\in\\mathcal{T}} h\_T^2 \\|f + \\Delta u\_h\\|^2\_{L^2(T)} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} h\_E \\|\[\\nabla u\_h\\cdot n\_E\]\\|^2\_{L^2(E)}, $
%[text] where $\\Delta u\_h = 0$ on every triangle for $P\_1$. In the contact set its volume term $h\_T\\|f\\|\_{L^2(T)}$ does not vanish, although the error is small there; estimators designed for the obstacle problem subtract the discrete contact force and come with proofs of reliability and efficiency.
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2\\Vert f+\\Delta u\_h\\Vert^2\_{L^2(T)}$ | `h^2*(f(x)+div(grad(u)))^2 dx` |
%[text] | $h\_E\\Vert\[\\nabla u\_h\\cdot n\_E\]\\Vert^2\_{L^2(E)}$ | `h^1*([grad(u)*n])^2 dsi` |
%[text:table]
%[text] ### Adaptive algorithm
%[text] Starting from the coarse mesh, the adaptive finite element method repeats
%[text] 1. SOLVE: compute $u\_h$ and $\\lambda$ with `pdaSet` on the current mesh,
%[text] 2. ESTIMATE: compute $\\eta\_T$ for all $T\\in\\mathcal{T}$,
%[text] 3. MARK: select a set $\\mathcal{M}\\subseteq\\mathcal{T}$ of minimal size with $\\theta\\,\\eta^2 \\le \\sum\_{T\\in\\mathcal{M}}\\eta\_T^2$ (`MarkBulk` with $\\theta = 1/2$),
%[text] 4. REFINE: refine the edges of the marked triangles, and further edges to keep the mesh conforming, by red-green-blue refinement \
%[text] until the number of degrees of freedom exceeds `minNrDoF`.
%[text] ### Expected convergence
%[text] For $u\\in H^2(\\Omega)$ the $P\_1$ approximation of the obstacle problem satisfies $\\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)} = O(h)$, i.e. $O(N^{-1/2})$ already on uniform meshes. The solution of this problem has no singularity, so adaptive refinement keeps the rate $N^{-1/2}$ and improves only the constant.
%%
%[text] ## Parameters
minNrDoF = 1e4;   % stop as soon as the number of degrees of freedom exceeds minNrDoF
mark = MarkBulk;  % bulk (Doerfler) marking with theta = 1/2
% mark = MarkUniform;  % uniform refinement, for comparison
%%
%[text] ## Problem data and initial mesh
%[text] Three uniform refinements of the coarse mesh of the square, two triangles, give the first mesh.
p = problemObstacleSquareExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform',3);
plot(mesh);
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the loop stops when the number of degrees of freedom exceeds `minNrDoF`.
nrDoF4lvl = [];
iter4lvl = [];
H1Error4lvl = [];
eta4lvl = [];
while true
    % finite element space
    fe = Pk(1,mesh);
    nrDoF4lvl(end+1) = fe.nrDoF;
    % assemble
    a = BilinearForm('grad(u) * grad(v) dx');
    A = a.assemble(fe);
    l = LinearForm('f(x) v dx');
    b = l.assemble(fe,'f',@p.f);
    % solve
    x = zeros(fe.nrDoF,1);
    x(fe.fixedNodes) = p.u4Db(fe.c4dof(fe.fixedNodes,:));
    b = b - A*x;
    free = fe.freeNodes;
    low = p.lower(fe.c4dof(free,:));  % obstacle chi at the free nodes
    up = p.upper(fe.c4dof(free,:));   % no upper obstacle: inf
    lambda = zeros(fe.nrDoF,1);
    [x(free),lambda(free),nrIter] = pdaSet(A(free,free),b(free),low,up);
    iter4lvl(end+1) = nrIter;
    % estimate
    eta = EstimatorForm(['h^2*(f(x)+div(grad(u)))^2 dx', ...
        '+ h^1*([grad(u)*n])^2 dsi']);
    eta4e = eta.assemble(x,fe,'f',@p.f);
    eta4lvl(end+1) = norm(eta4e);
    % errors
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(end+1) = h1.assemble(x,fe,'exact',@p.exactGradU);
    % stop
    if fe.nrDoF >= minNrDoF
        break
    end
    % mark
    mark = mark.mark(fe,eta4e);
    % refine
    mesh = mesh.refine(mark.marked4ed);
end
%%
%[text] ## Results
%[text] ### Convergence table
%[text] The second column is the number of active set iterations, the last one the ratio of the error and the estimator.
table(nrDoF4lvl',iter4lvl',H1Error4lvl',eta4lvl',(H1Error4lvl./eta4lvl)', ...
    VariableNames=["N" "active set iterations" "H1 error" "eta" "error/eta"])
%[text] The ratio of error and estimator stays almost constant. The number of active set iterations grows with $N$, because every solve starts from $x = 0$.
%%
%[text] ### Adaptive mesh
clf
plot(mesh);
%%
%[text] ### Discrete solution
clf
plot(fe,x);
view([40,35]);
colorbar;
shading interp;
%[text] The discrete solution lies on the obstacle $\\chi = 0$ in the disk $r\\le r\_0$ and rises towards the boundary.
%%
%[text] ### Discrete Lagrange multiplier
%[text] The discrete multiplier $\\lambda\_j$ approximates $\\int\_\\Omega\\lambda\\varphi\_j\\,dx$. Divided by $\\int\_\\Omega\\varphi\_j\\,dx$, the row sums of the mass matrix (mass lumping), it gives the nodal values of a piecewise linear approximation of the contact force $\\lambda = -\\Delta u - f$.
mass = BilinearForm('u * v dx');
M = mass.assemble(fe);
clf
plot(fe,lambda./sum(M,2));
view([40,35]);
colorbar;
shading interp;
%[text] The multiplier is nonnegative, vanishes outside the contact set and approximates $\\lambda = -f = 8r\_0^2(1+r\_0^2-r^2)$ in the disk $r\\le r\_0$; it jumps to zero at the free boundary.
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(u-u_h)||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,eta4lvl,'\eta');
hold off
%[text] Error and estimator decrease with the rate $1/2$, the optimal rate of $P\_1$ elements.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkUniform;` to refine uniformly: the rate is $N^{-1/2}$ as well, since $u\\in H^2(\\Omega)$; the adaptive meshes reach a slightly smaller error with the same number of degrees of freedom.
%[text] - Start the active set method from the solution on the previous mesh: set `x0 = zeros(mesh.nrNodes,1);` before the loop and `x0 = mesh.prolongation*x;` after `mesh = mesh.refine(...)`, and call `pdaSet(A(free,free),b(free),low,up,[],x0(free))`. The active set of the prolongated solution is almost right, and two or three iterations suffice on every mesh.
%[text] - Use `problemObstacleLshapeExact`, the obstacle problem on an L-shaped domain with the corner singularity $r^{2/3}$ of `exPoissonPkAdaptive`: adaptive refinement keeps the rate $N^{-1/2}$, uniform refinement converges more slowly, and its rate approaches $1/3$ only on much finer meshes.
%[text] - Use `problemObstacleSquareDist`, the upper obstacle $\\mathrm{dist}(x,\\partial\\Omega)$ and $f = 7$ of elastic-plastic torsion; it has no exact solution, so remove the step errors, the error columns of the table and the error curve, and compare $\\eta$ only. For an upper obstacle the multiplier is nonpositive. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
