%[text] # Adaptive $P\_k$ finite elements for the Laplace eigenvalue problem
%[text] Approximates the smallest Dirichlet eigenvalue of the Laplacian on the L-shaped domain, whose eigenfunction has a singular gradient at the reentrant corner, with conforming $P\_k$ finite elements on adaptively refined meshes. New are the generalized matrix eigenvalue problem and the error estimator with the right-hand side $\\lambda\_h u\_h$; the eigenvalue converges twice as fast as the eigenfunction.
%[text] **Builds on:** `exPoissonPkAdaptive`. **Next:** `exEigenvalueDiffusionPkAdaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] On the L-shaped domain $\\Omega = (-1,1)^2\\setminus(\[0,1\]\\times\[-1,0\])$ find the eigenpairs $(\\lambda,u)$ with $u\\neq 0$ and
%[text]{"align":"center"} $ -\\Delta u = \\lambda u \\textrm{ in } \\Omega, \\qquad u = 0 \\textrm{ on } \\partial\\Omega. $
%[text] The eigenvalues form a sequence $0 \< \\lambda\_1 \< \\lambda\_2 \\le \\lambda\_3 \\le \\dots \\to \\infty$, repeated according to their multiplicity. The problem object `problemScalarEllipticEigenvalueLshape` provides the geometry with the Dirichlet condition on the whole boundary and reference values `p.eigenvalue(m)` of the first ten eigenvalues: $\\lambda\_1 = 9.6397238440219\\dots$ to double precision from a computation with 100 digits, the exact values $\\lambda\_3 = 2\\pi^2$ and $\\lambda\_8 = \\lambda\_9 = 5\\pi^2$, whose eigenfunctions such as $\\sin(\\pi x\_1)\\sin(\\pi x\_2)$ vanish on all lines $x\_1, x\_2\\in\\{-1,0,1\\}$, and numerical approximations of the others. The first eigenfunction behaves like the solution of `exPoissonPkAdaptive` at the reentrant corner, $u\\sim r^{2/3}\\sin(2\\varphi/3)$ in polar coordinates $(r,\\varphi)$, and its gradient is unbounded there.
%[text] ### Weak formulation
%[text] With the bilinear forms $a(u,v) = \\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx$ and $b(u,v) = \\int\_\\Omega u v\\,dx$ the eigenvalue problem reads: find $\\lambda\\in\\mathbb{R}$ and $u\\in H^1\_0(\\Omega)$ with $\\|u\\|\_{L^2(\\Omega)} = 1$ and
%[text]{"align":"center"} $ a(u,v) = \\lambda\\, b(u,v) \\quad \\textrm{for all } v\\in H^1\_0(\\Omega). $
%[text] ### Discretization
%[text] The discrete space $V\_h = P\_k(\\mathcal{T})\\cap H^1\_0(\\Omega)$ (`Pk(degree,mesh)`, the boundary degrees of freedom `fe.fixedNodes` set to zero) gives the discrete eigenpairs $(\\lambda\_h,u\_h)$ with $a(u\_h,v\_h) = \\lambda\_h\\, b(u\_h,v\_h)$ for all $v\_h\\in V\_h$. With the stiffness matrix $A$ and the mass matrix $B$ restricted to the free degrees of freedom this is the generalized matrix eigenvalue problem
%[text]{"align":"center"} $ A x = \\lambda\_h B x $
%[text] with symmetric positive definite matrices $A$ and $B$; `eigs(A,B,m,'smallestabs')` computes its $m$ smallest eigenvalues. The coefficient vector $x$ of the $m$-th eigenvalue is scaled to $x^T B x = \\|u\_h\\|^2\_{L^2(\\Omega)} = 1$. By the min-max principle the conforming discrete eigenvalues are upper bounds, $\\lambda\_m \\le \\lambda\_{h,m}$.
%[text:table]
%[text] | Matrix | Integral | fem4m form |
%[text] | --- | --- | --- |
%[text] | $A$ | $\\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx$ | `grad(u) * grad(v) dx` |
%[text] | $B$ | $\\int\_\\Omega u v\\,dx$ | `u * v dx` |
%[text:table]
%[text] ### Error estimator
%[text] The discrete eigenfunction is the Galerkin approximation of the solution $w\\in H^1\_0(\\Omega)$ of the source problem $-\\Delta w = \\lambda\_h u\_h$, and $w$ differs from $u$ only by higher order terms. The residual estimator of this source problem, with the right-hand side $f = \\lambda\_h u\_h$, reads
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{T\\in\\mathcal{T}} h\_T^2 \\|\\lambda\_h u\_h + \\Delta u\_h\\|^2\_{L^2(T)} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} h\_E \\|\[\\nabla u\_h\\cdot n\_E\]\\|^2\_{L^2(E)}. $
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2\\Vert \\lambda\_h u\_h+\\Delta u\_h\\Vert^2\_{L^2(T)}$ | `h^2*(f(x)+div(grad(u)))^2 dx` |
%[text] | $h\_E\\Vert\[\\nabla u\_h\\cdot n\_E\]\\Vert^2\_{L^2(E)}$ | `h^1*([grad(u)*n])^2 dsi` |
%[text:table]
%[text] The function $f = \\lambda\_h u\_h$ in the form depends on the discrete solution. A function handle with the three arguments `(pts,ptsRef,parts)`, the points, the reference point and the element indices, can evaluate it with `fe.getValue(pts,ptsRef,parts,x)`. For a simple eigenvalue and the sign of $u\_h$ chosen such that $b(u,u\_h) \> 0$, the estimator bounds the error of the eigenfunction, $\\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)} \\le C\\,\\eta$, up to higher order terms. The eigenvalue error satisfies the identity
%[text]{"align":"center"} $ \\lambda\_h - \\lambda = \\|\\nabla(u-u\_h)\\|^2\_{L^2(\\Omega)} - \\lambda\\|u-u\_h\\|^2\_{L^2(\\Omega)}, $
%[text] in which the $L^2$ error is of higher order, so $\\lambda\_h - \\lambda$ behaves like $\\eta^2$.
%[text] ### Adaptive algorithm
%[text] Starting from the coarse mesh, the adaptive finite element method repeats
%[text] 1. SOLVE: compute the $m$-th discrete eigenpair $(\\lambda\_h,u\_h)$ on the current mesh,
%[text] 2. ESTIMATE: compute $\\eta\_T$ for all $T\\in\\mathcal{T}$,
%[text] 3. MARK: select a set $\\mathcal{M}\\subseteq\\mathcal{T}$ of minimal size with $\\theta\\,\\eta^2 \\le \\sum\_{T\\in\\mathcal{M}}\\eta\_T^2$ (`MarkBulk` with $\\theta = 1/2$),
%[text] 4. REFINE: refine the edges of the marked triangles, and further edges to keep the mesh conforming, by red-green-blue refinement \
%[text] until the number of degrees of freedom exceeds `minNrDoF`.
%[text] ### Expected convergence
%[text] On uniformly refined meshes the corner singularity limits the rates to $\\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)} = O(N^{-1/3})$ and $\\lambda\_h - \\lambda = O(N^{-2/3})$ for every degree $k$. With bulk marking, a sufficiently small $\\theta$ and a sufficiently fine initial mesh, the adaptive algorithm provably recovers the optimal rates $N^{-k/2}$ of the eigenfunction and $N^{-k}$ of the eigenvalue; in practice $\\theta = 1/2$ achieves them as well, and $\\eta^2$ decreases like the eigenvalue error.
%%
%[text] ## Parameters
degree = 1;       % polynomial degree k
m = 1;            % compute the m-th eigenvalue lambda_m
minNrDoF = 1e4;   % stop as soon as the number of degrees of freedom exceeds minNrDoF
mark = MarkBulk;  % bulk (Doerfler) marking with theta = 1/2
% mark = MarkUniform;  % uniform refinement, for comparison
%%
%[text] ## Problem data and initial mesh
p = problemScalarEllipticEigenvalueLshape();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform',2);
plot(mesh);
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the loop stops when the number of degrees of freedom exceeds `minNrDoF`.
nrDoF4lvl = [];
lambda4lvl = [];
lambdaError4lvl = [];
eta4lvl = [];
while true
    % finite element space
    fe = Pk(degree,mesh);
    nrDoF4lvl(end+1) = fe.nrDoF;
    % assemble
    a = BilinearForm('grad(u) * grad(v) dx');
    A = a.assemble(fe);
    b = BilinearForm('u * v dx');
    B = b.assemble(fe);
    % solve
    free = fe.freeNodes;
    [V,D] = eigs(A(free,free),B(free,free),m,'smallestabs');
    [lambda4mode,ind] = sort(diag(D));
    lambdaH = lambda4mode(m);
    x = zeros(fe.nrDoF,1);
    x(free) = V(:,ind(m));
    x = x/sqrt(x'*B*x);  % ||u_h||_L2 = 1
    lambda4lvl(end+1) = lambdaH;
    % estimate
    eta = EstimatorForm(['h^2*(f(x)+div(grad(u)))^2 dx', ...
        '+ h^1*([grad(u)*n])^2 dsi']);
    eta4e = eta.assemble(x,fe,'f',@(pts,ptsRef,parts)(lambdaH*fe.getValue(pts,ptsRef,parts,x)));
    eta4lvl(end+1) = norm(eta4e);
    % errors
    lambdaError4lvl(end+1) = lambdaH - p.eigenvalue(m);
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
%[text] The rate is the slope of the eigenvalue error over $N$ in a log-log scale between two consecutive meshes; the last column compares the eigenvalue error with $\\eta^2$.
rate = [NaN, -diff(log(lambdaError4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',lambda4lvl',lambdaError4lvl',rate',(eta4lvl.^2)',(lambdaError4lvl./eta4lvl.^2)', ...
    VariableNames=["N" "lambda_h" "lambda_h - lambda" "rate" "eta^2" "error/eta^2"])
%[text] The discrete eigenvalues decrease monotonically towards $\\lambda\_1$ from above, and the ratio of the eigenvalue error and $\\eta^2$ stays bounded.
%%
%[text] ### Adaptive mesh
clf
plot(mesh);
%[text] The refinement concentrates at the reentrant corner, where the gradient of the eigenfunction is singular.
%%
%[text] ### Discrete eigenfunction
clf
plot(fe,x);
view([40,35]);
colorbar;
shading interp;
%[text] The sign of an eigenfunction is arbitrary: `eigs` may return $-u\_h$.
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,lambdaError4lvl,'\lambda_h-\lambda');
hold on
plotConvergence(nrDoF4lvl,eta4lvl.^2,'\eta^2');
hold off
%[text] The eigenvalue error and $\\eta^2$ decrease with the same rate, close to $k$: twice the rate $k/2$ of the energy error of the eigenfunction, and well above the rate $2/3$ of uniform refinement.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkUniform;` to refine uniformly: the eigenvalue error decreases only like $N^{-2/3}$.
%[text] - Set `degree = 2` or `3`: the adaptive eigenvalue error decreases like $N^{-k}$; the reference value $\\lambda\_1$ is known to double precision, so the rate stays visible far below $10^{-8}$.
%[text] - Set `m = 3`: the eigenfunction $\\sin(\\pi x\_1)\\sin(\\pi x\_2)$ of $\\lambda\_3 = 2\\pi^2$ is smooth, so uniform refinement already converges with the rate $N^{-k}$.
%[text] - Use `problemScalarEllipticEigenvalueSlit`, the slit domain, whose first eigenfunction behaves like $r^{1/2}$ at the tip of the slit: uniform refinement converges only like $N^{-1/2}$. `problemScalarEllipticEigenvalueIsospectralDrumsA` and `...DrumsB` are two domains of different shape with the same eigenvalues.
%[text] - `exEigenvalueDiffusionPkAdaptive` computes a multiple eigenvalue of an anisotropic diffusion operator together with the errors of the eigenfunctions. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
