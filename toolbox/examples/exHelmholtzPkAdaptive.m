%[text] # Adaptive $P\_k$ finite elements for the Helmholtz problem
%[text] Solves the Helmholtz equation with an impedance boundary condition on the L-shaped domain, whose solution is a Bessel function with a singularity at the reentrant corner, with conforming $P\_k$ finite elements on adaptively refined meshes. New are the complex-valued, indefinite linear system and the resolution of the wave length: the error decreases only after the mesh resolves the waves.
%[text] **Builds on:** `exPoissonPkAdaptive`. **Next:** `exObstacleP1Adaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] The Helmholtz equation describes time-harmonic waves $\\mathrm{Re}(u(x)e^{i\\omega t})$ with the angular frequency $\\omega$, which is the wave number for a medium with refractive index one, and the wave length $2\\pi/\\omega$. On the L-shaped domain $\\Omega = (-1,1)^2\\setminus(\[0,1\]\\times\[-1,0\])$, with the two edges $\\Gamma\_D$ at the reentrant corner and the outer edges $\\Gamma\_N$, find $u$ with
%[text]{"align":"center"} $ -\\Delta u - \\omega^2 u = f \\textrm{ in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\Gamma\_D, \\qquad \\nabla u\\cdot n + i\\omega u = g \\textrm{ on } \\Gamma\_N. $
%[text] The impedance condition on $\\Gamma\_N$ is the simplest absorbing boundary condition: it lets waves that hit the boundary perpendicularly leave the domain without reflection.
%[text] The problem object `problemHelmholtzLshapeExact(omega,2/3)` provides $f = 0$ and the data $u\_D = u$ and $g = \\nabla u\\cdot n + i\\omega u$ (`p.g4Rb`) of the exact solution, in polar coordinates $(r,\\varphi)$ with $\\varphi\\in\[0,3\\pi/2\]$,
%[text]{"align":"center"} $ u(r,\\varphi) = J\_{2/3}(\\omega r)\\sin(2\\varphi/3), $
%[text] with the Bessel function $J\_{2/3}$ of the first kind. Bessel's equation makes $J\_\\nu(\\omega r)\\sin(\\nu\\varphi)$ a solution of $-\\Delta u - \\omega^2 u = 0$ for every order $\\nu$, and the order $\\nu = 2/3$ makes it vanish on $\\Gamma\_D$, $u\_D = 0$. Near the corner $J\_{2/3}(\\omega r)\\approx(\\omega r/2)^{2/3}/\\Gamma(5/3)$, so $u$ has the corner singularity $r^{2/3}\\sin(2\\varphi/3)$ of `exPoissonPkAdaptive`; away from the corner it oscillates with the wave length $2\\pi/\\omega$ and decays like $r^{-1/2}$.
%[text] ### Weak formulation
%[text] Multiplication with the complex conjugate $\\overline{v}$ of a test function, integration by parts and the impedance condition $\\nabla u\\cdot n = g - i\\omega u$ on $\\Gamma\_N$ give: find $u\\in H^1(\\Omega)$ with $u = u\_D$ on $\\Gamma\_D$ and
%[text]{"align":"center"} $ \\int\_\\Omega \\nabla u\\cdot\\nabla\\overline{v}\\,dx - \\omega^2\\int\_\\Omega u\\overline{v}\\,dx + i\\omega\\int\_{\\Gamma\_N} u\\overline{v}\\,ds = \\int\_\\Omega f\\overline{v}\\,dx + \\int\_{\\Gamma\_N} g\\overline{v}\\,ds \\quad \\textrm{for all } v\\in H^1\_D(\\Omega). $
%[text] The sesquilinear form on the left is indefinite: the term $-\\omega^2\\int\_\\Omega u\\overline{v}\\,dx$ dominates the Dirichlet energy for functions that vary slowly compared with the wave length, such as constants. The imaginary impedance term makes the problem uniquely solvable for every $\\omega \> 0$, while with $u = u\_D$ on the whole boundary $\\omega^2$ must not be a Dirichlet eigenvalue of $-\\Delta$ (resonance).
%[text] ### Discretization
%[text] The discrete space $V\_h = P\_k(\\mathcal{T})\\cap C(\\overline\\Omega)$ (`Pk(degree,mesh)`) has real basis functions $\\varphi\_j$, so the conjugation drops out of the matrix: with the stiffness matrix $S$, the mass matrix $M$ and the boundary mass matrix $R$ on $\\Gamma\_N$ the system matrix
%[text]{"align":"center"} $ A = S - \\omega^2 M + i\\omega R $
%[text] is complex symmetric, $A^T = A$, but neither Hermitian nor positive definite. The right-hand side is complex as well, since $g$ is; backslash solves the sparse complex system directly. The constant coefficients enter the form strings as `omega2(x)` and `iomega(x)` with the values `omega^2` and `1i*omega` passed to `assemble`.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx - \\omega^2\\int\_\\Omega u v\\,dx + i\\omega\\int\_{\\Gamma\_N} u v\\,ds$ | `grad(u)*grad(v) dx - omega2(x)*u*v dx + iomega(x)*u*v dsn` |
%[text] | $\\int\_\\Omega f v\\,dx + \\int\_{\\Gamma\_N} g v\\,ds$ | `f(x) v dx + g(x) v dsn` |
%[text:table]
%[text] ### Resolution condition
%[text] The mesh has to resolve the waves. A fixed number of elements per wave length, $\\omega h/k$ small, bounds the error of the best approximation by $(\\omega h)^k$ relative to the solution. The Galerkin error is larger, by the pollution effect: the discrete waves travel with a slightly wrong speed, and the phase error accumulates over the domain. Before the asymptotic range the relative energy error behaves like $(\\omega h)^k + \\omega(\\omega h)^{2k}$, and the Galerkin solution is quasi-optimal once $\\omega^{k+1}h^k$ is small. Higher degrees $k$ reduce the pollution. Here $\\omega = 20$ gives the wave length $2\\pi/20\\approx 0.31$, and the first mesh, whose triangles have legs of length $1/4$, has about one element per wave length: its discrete solution has no accuracy at all.
%[text] ### Error estimator
%[text] The residual-based estimator contains the residual of the Helmholtz equation, the jumps of the normal derivative and the residual of the impedance condition, all measured in absolute values of complex numbers:
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{T\\in\\mathcal{T}} h\_T^2 \\|f + \\Delta u\_h + \\omega^2 u\_h\\|^2\_{L^2(T)} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} h\_E \\|\[\\nabla u\_h\\cdot n\_E\]\\|^2\_{L^2(E)} + \\sum\_{E\\in\\mathcal{E}(\\Gamma\_N)} h\_E \\|g - \\nabla u\_h\\cdot n - i\\omega u\_h\\|^2\_{L^2(E)}. $
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2\\Vert f+\\Delta u\_h+\\omega^2 u\_h\\Vert^2\_{L^2(T)}$ | `h^2*(f(x)+div(grad(u))+omega2(x)*u)^2 dx` |
%[text] | $h\_E\\Vert\[\\nabla u\_h\\cdot n\_E\]\\Vert^2\_{L^2(E)}$ | `h^1*([grad(u)*n])^2 dsi` |
%[text] | $h\_E\\Vert g-\\nabla u\_h\\cdot n-i\\omega u\_h\\Vert^2\_{L^2(E)}$ | `h^1*(g(x)-grad(u)*n-iomega(x)*u)^2 dsn` |
%[text:table]
%[text] Once the mesh resolves the waves, the estimator is reliable and efficient with constants that do not grow with $\\omega$; on coarse meshes it can underestimate the error, because it does not see the pollution.
%[text] ### Adaptive algorithm
%[text] Starting from the coarse mesh, the adaptive finite element method repeats
%[text] 1. SOLVE: compute $u\_h$ on the current mesh,
%[text] 2. ESTIMATE: compute $\\eta\_T$ for all $T\\in\\mathcal{T}$,
%[text] 3. MARK: select a set $\\mathcal{M}\\subseteq\\mathcal{T}$ of minimal size with $\\theta\\,\\eta^2 \\le \\sum\_{T\\in\\mathcal{M}}\\eta\_T^2$ (`MarkBulk` with $\\theta = 1/2$),
%[text] 4. REFINE: refine the edges of the marked triangles, and further edges to keep the mesh conforming, by red-green-blue refinement \
%[text] until the number of degrees of freedom exceeds `minNrDoF`.
%[text] ### Expected convergence
%[text] On the first meshes the error stays at the size of the solution: the pre-asymptotic phase, in which the refinement resolves the waves in the whole domain. Afterwards uniform refinement converges like $N^{-1/3}$ in the energy norm because of the corner singularity, and adaptive refinement with the optimal rate $N^{-k/2}$; the $L^2$ error converges faster.
%%
%[text] ## Parameters
omega = 20;       % angular frequency omega, the wave number
degree = 2;       % polynomial degree k
minNrDoF = 2e4;   % stop as soon as the number of degrees of freedom exceeds minNrDoF
mark = MarkBulk;  % bulk (Doerfler) marking with theta = 1/2
% mark = MarkUniform;  % uniform refinement, for comparison
%%
%[text] ## Problem data and initial mesh
%[text] The geometry of the problem object has Dirichlet conditions on all edges. The second line sets the impedance condition on the outer edges, the second boundary part of the domain `'Lshape'`, as Neumann sides: the side domain `dsn` of the forms. Two uniform refinements of the coarse mesh, six triangles, give the first mesh.
p = problemHelmholtzLshapeExact(omega,2/3);
p.geom = Geometry(struct('name','Lshape','bdry','DN'));
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform',2);
plot(mesh);
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the loop stops when the number of degrees of freedom exceeds `minNrDoF`.
nrDoF4lvl = [];
L2Error4lvl = [];
H1Error4lvl = [];
eta4lvl = [];
while true
    % finite element space
    fe = Pk(degree,mesh);
    nrDoF4lvl(end+1) = fe.nrDoF;
    % assemble
    a = BilinearForm('grad(u)*grad(v) dx - omega2(x)*u*v dx + iomega(x)*u*v dsn');
    A = a.assemble(fe,'omega2',omega^2,'iomega',1i*omega);
    l = LinearForm('f(x) v dx + g(x) v dsn');
    b = l.assemble(fe,'f',@p.f,'g',@(x,n)p.g4Rb(x,n));
    % solve
    x = zeros(fe.nrDoF,1);
    x(fe.fixedNodes) = p.u4Db(fe.c4dof(fe.fixedNodes,:));
    b = b - A*x;
    free = fe.freeNodes;
    x(free) = A(free,free)\b(free);
    % estimate
    eta = EstimatorForm(['h^2*(f(x)+div(grad(u))+omega2(x)*u)^2 dx', ...
        '+ h^1*([grad(u)*n])^2 dsi', ...
        '+ h^1*(g(x)-grad(u)*n-iomega(x)*u)^2 dsn']);
    eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g4Rb(x,n),'omega2',omega^2,'iomega',1i*omega);
    eta4lvl(end+1) = norm(eta4e);
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(end+1) = l2.assemble(x,fe,2*(degree+3),'exact',@p.exactU);
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(end+1) = h1.assemble(x,fe,2*(degree+3),'exact',@p.exactGradU);
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
%[text] The errors are integrated with a higher quadrature degree, $2(k+3)$, since the solution oscillates. The last column is the ratio of the energy error and the estimator.
table(nrDoF4lvl',L2Error4lvl',H1Error4lvl',eta4lvl',(H1Error4lvl./eta4lvl)', ...
    VariableNames=["N" "L2 error" "H1 error" "eta" "error/eta"])
%[text] On the first meshes the error is about as large as the solution itself. Once the waves are resolved, the error drops quickly and the ratio of error and estimator settles at a constant.
%%
%[text] ### Adaptive mesh
clf
plot(mesh);
%[text] The mesh is refined in the whole domain, to resolve the waves, and strongly at the reentrant corner.
%%
%[text] ### Discrete solution
%[text] The real part of $u\_h$; the third argument of `plot` refines every triangle twice uniformly for the plot, so that the quadratic polynomials are drawn smoothly.
clf
plot(fe,real(x),2);
view([40,35]);
colorbar;
shading interp;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(u-u_h)||_{L^2}');
plotConvergence(nrDoF4lvl,eta4lvl,'\eta');
hold off
%[text] The legend slopes average over the pre-asymptotic phase and the asymptotic range. In the asymptotic range, the last meshes of the table, the energy error and the estimator decrease with the optimal rate $k/2 = 1$.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkUniform;` to refine uniformly: after the pre-asymptotic phase the rate of the energy error decreases towards $1/3$.
%[text] - Set `degree = 1`: the pollution effect is much stronger, the error even grows on the first meshes, and the pre-asymptotic phase lasts much longer.
%[text] - Double the frequency, `omega = 40`, and raise `minNrDoF`: the pre-asymptotic phase needs many more degrees of freedom.
%[text] - Remove the line `p.geom = ...` to impose the Dirichlet condition on all edges. Then $\\omega^2 = 400$ lies between the Dirichlet eigenvalues $40\\pi^2\\approx 394.8$ and $41\\pi^2\\approx 404.7$ of the L-shaped domain, the problem is close to resonance, the matrix is real, and on the first mesh the error is about 15 times larger than the solution, on the next meshes still a multiple of it. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
