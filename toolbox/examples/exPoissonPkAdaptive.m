%[text] # Adaptive $P\_k$ finite elements for the Poisson problem
%[text] Solves the Poisson problem on the L-shaped domain, whose solution has a singular gradient at the reentrant corner, with conforming $P\_k$ finite elements on adaptively refined meshes. A residual-based error estimator steers the refinement towards the singularity, which recovers the optimal convergence rate.
%[text] **Builds on:** `exPoissonP1`. **Next:** `exPoissonHPAdaptive`, `exPoisson3dPkAdaptive`, `exEigenvaluePkAdaptive`, `exElasticityPkAdaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] The L-shaped domain $\\Omega = (-1,1)^2\\setminus(\[0,1\]\\times\[-1,0\])$ has a reentrant corner at the origin. Its boundary consists of the two edges $\\Gamma\_D$ at the reentrant corner and the remaining edges $\\Gamma\_N$. Given $f$, $u\_D$ and $g$, find $u$ with
%[text]{"align":"center"} $ -\\Delta u = f \\textrm{ in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\Gamma\_D, \\qquad \\nabla u\\cdot n = g \\textrm{ on } \\Gamma\_N. $
%[text] The problem object `problemScalarEllipticLshapeExact` provides $f = 0$ and the boundary data of the exact solution, in polar coordinates $(r,\\varphi)$ with $\\varphi\\in\[0,2\\pi)$,
%[text]{"align":"center"} $ u(r,\\varphi) = r^{2/3}\\sin(2\\varphi/3). $
%[text] It vanishes on $\\Gamma\_D$, so $u\_D = 0$. Its gradient behaves like $r^{-1/3}$ and is unbounded at the origin: $u$ belongs to $H^{1+s}(\\Omega)$ only for $s \< 2/3$.
%[text] ### Weak formulation
%[text] Find $u\\in H^1(\\Omega)$ with $u = u\_D$ on $\\Gamma\_D$ and
%[text]{"align":"center"} $ \\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx = \\int\_\\Omega f v\\,dx + \\int\_{\\Gamma\_N} g v\\,ds \\quad \\textrm{for all } v\\in H^1\_D(\\Omega) = \\{v\\in H^1(\\Omega) : v = 0 \\textrm{ on } \\Gamma\_D\\}. $
%[text] ### Discretization
%[text] The discrete space $V\_h = P\_k(\\mathcal{T})\\cap C(\\overline\\Omega)$ of continuous piecewise polynomials of degree $k$ (`Pk(degree,mesh)`) has one degree of freedom per Lagrange node `fe.c4dof`. The degrees of freedom on $\\Gamma\_D$ (`fe.fixedNodes`) take the values of $u\_D$ at their nodes, the others (`fe.freeNodes`) are the unknowns. The Neumann data enter the right-hand side through the boundary integral over $\\Gamma\_N$, the domain `dsn` of the form string; $g$ depends on the outer normal, `@(x,n)p.g(x,n)`.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx$ | `grad(u) * grad(v) dx` |
%[text] | $\\int\_\\Omega f v\\,dx + \\int\_{\\Gamma\_N} g v\\,ds$ | `f(x) v dx + g(x) v dsn` |
%[text:table]
%[text] ### Error estimator
%[text] The residual-based error estimator $\\eta^2 = \\sum\_{T\\in\\mathcal{T}}\\eta\_T^2$ collects the residuals of the equation on the triangles, the jumps $\[\\nabla u\_h\\cdot n\_E\]$ of the normal derivative across the interior edges $\\mathcal{E}(\\Omega)$ and the residuals of the Neumann condition on the edges $\\mathcal{E}(\\Gamma\_N)$:
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{T\\in\\mathcal{T}} h\_T^2 \\|f + \\Delta u\_h\\|^2\_{L^2(T)} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} h\_E \\|\[\\nabla u\_h\\cdot n\_E\]\\|^2\_{L^2(E)} + \\sum\_{E\\in\\mathcal{E}(\\Gamma\_N)} h\_E \\|g - \\nabla u\_h\\cdot n\\|^2\_{L^2(E)}. $
%[text] `EstimatorForm` assembles the local contributions $\\eta\_T$ (`eta4e`); an edge term counts one half for each of the two neighbouring triangles.
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2\\Vert f+\\Delta u\_h\\Vert^2\_{L^2(T)}$ | `h^2*(f(x)+div(grad(u)))^2 dx` |
%[text] | $h\_E\\Vert\[\\nabla u\_h\\cdot n\_E\]\\Vert^2\_{L^2(E)}$ | `h^1*([grad(u)*n])^2 dsi` |
%[text] | $h\_E\\Vert g-\\nabla u\_h\\cdot n\\Vert^2\_{L^2(E)}$ | `h^1*(g(x)-grad(u)*n)^2 dsn` |
%[text:table]
%[text] The estimator is reliable and efficient: $\\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)} \\le C\_{\\textrm{rel}}\\,\\eta$ and, up to oscillations of the data, $\\eta \\le C\_{\\textrm{eff}}\\,\\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)}$.
%[text] ### Adaptive algorithm
%[text] Starting from the coarse mesh, the adaptive finite element method repeats
%[text] 1. SOLVE: compute $u\_h$ on the current mesh,
%[text] 2. ESTIMATE: compute $\\eta\_T$ for all $T\\in\\mathcal{T}$,
%[text] 3. MARK: select a set $\\mathcal{M}\\subseteq\\mathcal{T}$ of minimal size with $\\theta\\,\\eta^2 \\le \\sum\_{T\\in\\mathcal{M}}\\eta\_T^2$ (bulk or Dörfler marking, `MarkBulk` with $\\theta = 1/2$),
%[text] 4. REFINE: refine the edges of the marked triangles, and further edges to keep the mesh conforming, by red-green-blue refinement \
%[text] until the number of degrees of freedom exceeds `minNrDoF`.
%[text] ### Expected convergence
%[text] On uniformly refined meshes the singularity limits the rate to $\\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)} = O(h^{2/3}) = O(N^{-1/3})$ for every degree $k$. With bulk marking and a sufficiently small $\\theta$ the adaptive algorithm provably converges with the optimal rate $N^{-k/2}$ of a smooth solution; in practice $\\theta = 1/2$ achieves this rate as well, and the estimator decreases with the same rate.
%%
%[text] ## Parameters
degree = 1;       % polynomial degree k
minNrDoF = 1e4;   % stop as soon as the number of degrees of freedom exceeds minNrDoF
mark = MarkBulk;  % bulk (Doerfler) marking with theta = 1/2
% mark = MarkUniform;  % uniform refinement, for comparison
%%
%[text] ## Problem data and initial mesh
p = problemScalarEllipticLshapeExact();
mesh = Mesh(p.geom);
plot(mesh);
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the loop stops when the number of degrees of freedom exceeds `minNrDoF`.
nrDoF4lvl = [];
H1Error4lvl = [];
eta4lvl = [];
while true
    % finite element space
    fe = Pk(degree,mesh);
    nrDoF4lvl(end+1) = fe.nrDoF;
    % assemble
    a = BilinearForm('grad(u) * grad(v) dx');
    A = a.assemble(fe);
    l = LinearForm('f(x) v dx + g(x) v dsn');
    b = l.assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
    % solve
    x = zeros(fe.nrDoF,1);
    x(fe.fixedNodes) = p.u4Db(fe.c4dof(fe.fixedNodes,:));
    b = b - A*x;
    free = fe.freeNodes;
    x(free) = A(free,free)\b(free);
    % estimate
    eta = EstimatorForm(['h^2*(f(x)+div(grad(u)))^2 dx', ...
        '+ h^1*([grad(u)*n])^2 dsi', ...
        '+ h^1*(g(x)-grad(u)*n)^2 dsn']);
    eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g(x,n));
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
%[text] The last column is the ratio of the error and the estimator.
table(nrDoF4lvl',H1Error4lvl',eta4lvl',(H1Error4lvl./eta4lvl)', ...
    VariableNames=["N" "H1 error" "eta" "error/eta"])
%[text] The ratio stays bounded and roughly constant: the estimator reflects the error on every mesh.
%%
%[text] ### Adaptive mesh
clf
plot(mesh);
%[text] The refinement concentrates at the reentrant corner, where the gradient of $u$ is singular.
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
hold on
plotConvergence(nrDoF4lvl,eta4lvl,'\eta');
hold off
%[text] Error and estimator decrease with the same rate. The least-squares slopes in the legend are close to $k/2$, the rate of a smooth solution, and clearly above the rate $1/3$ of uniform refinement.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkUniform;` to refine uniformly: the rate drops to $N^{-1/3}$.
%[text] - Set `degree = 2` or `3`: the adaptive rate grows to $N^{-k/2}$, while uniform refinement stays at $N^{-1/3}$.
%[text] - Use `problemScalarEllipticSlitExact`, the slit domain with the stronger singularity $u\\sim r^{1/4}$: uniform refinement converges only like $N^{-1/8}$.
%[text] - Mark with `MarkMaximum` instead of `MarkBulk` and compare the meshes.
%[text] - `exPoissonHPAdaptive` adapts the polynomial degree as well and converges exponentially; `exPoisson3dPkAdaptive` is the same algorithm in 3D. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
