%[text] # Hellan-Herrmann-Johnson mixed finite elements for the biharmonic problem
%[text] Solves the biharmonic equation of a clamped plate in a mixed formulation with the Hessian as a second unknown. The Hellan-Herrmann-Johnson (HHJ) elements approximate the Hessian by symmetric matrix fields with continuous normal-normal components and the deflection by continuous Lagrange elements, so that no $C^1$ elements as in `exBiharmonicArgyris` are needed.
%[text] **Builds on:** `exPoissonRTAdaptive`, `exBiharmonicArgyris`. **Next:** `exStokesTaylorHood`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] As in `exBiharmonicArgyris`, the deflection $u$ of a clamped plate under the load $f$ solves
%[text]{"align":"center"} $ \\Delta^2 u = f \\textrm{ in } \\Omega, \\qquad u = 0 \\textrm{ and } \\nabla u\\cdot n = 0 \\textrm{ on } \\partial\\Omega, $
%[text] and the problem object `problemBiharmonicSquareClampedExact` provides the load of the exact solution $u = x\_1^2(1-x\_1)^2\\,x\_2^2(1-x\_2)^2$ on the unit square. With the Hessian $\\sigma = D^2 u$, the curvature tensor of the plate (up to constants the bending moments for the Poisson ratio 0), the equation splits into the system
%[text]{"align":"center"} $ \\sigma - D^2 u = 0, \\qquad \\mathrm{div}\\,\\mathrm{div}\\,\\sigma = f \\quad \\textrm{in } \\Omega. $
%[text] ### Mixed formulation
%[text] Let $\\mathcal{E}$ be the set of all edges of the triangulation and $\\mathcal{E}(\\Omega)$ that of the interior edges. Let $\\tau$ be a symmetric matrix field that is a polynomial on every triangle and whose normal-normal component $n\_E\\cdot\\tau n\_E$ is continuous across every interior edge $E$ with normal $n\_E$ and tangent $t\_E$. Integration by parts on every triangle $T$ gives
%[text]{"align":"center"} $ \\int\_T D^2 u:\\tau\\,dx = -\\int\_T \\nabla u\\cdot\\mathrm{div}\\,\\tau\\,dx + \\int\_{\\partial T} \\left(\\partial\_n u\\,(n\\cdot\\tau n) + \\partial\_t u\\,(t\\cdot\\tau n)\\right) ds $
%[text] with the normal and tangent of $\\partial T$. In the sum over all triangles the normal parts cancel on the interior edges, since $n\\cdot\\tau n$ is continuous and $\\partial\_n u$ changes its sign with $n$, and vanish on the boundary, since $\\partial\_n u = 0$. The tangential parts leave the jumps $\[t\_E\\cdot\\tau n\_E\]$ of the normal-tangential component; they vanish on the boundary as well, since $u = 0$ there. With the bilinear form
%[text]{"align":"center"} $ b(\\tau,v) = \\sum\_{T\\in\\mathcal{T}} \\int\_T \\mathrm{div}\\,\\tau\\cdot\\nabla v\\,dx - \\sum\_{E\\in\\mathcal{E}} \\int\_E \[t\_E\\cdot\\tau n\_E\]\\,\\partial\_{t}v\\,ds $
%[text] the first equation becomes $\\int\_\\Omega \\sigma:\\tau\\,dx + b(\\tau,u) = 0$. For continuous $v$ with $v = 0$ on $\\partial\\Omega$, $b(\\sigma,v)$ is the integral of $-\\mathrm{div}\\,\\mathrm{div}\\,\\sigma$ against $v$ in the distributional sense, so the second equation becomes $b(\\sigma,v) = -\\int\_\\Omega f v\\,dx$. The mixed problem seeks $\\sigma$ and $u\\in H^1\_0(\\Omega)$ with
%[text]{"align":"center"} $ \\int\_\\Omega \\sigma:\\tau\\,dx + b(\\tau,u) = 0 \\quad \\textrm{for all } \\tau, \\qquad b(\\sigma,v) = -\\int\_\\Omega f v\\,dx \\quad \\textrm{for all } v\\in H^1\_0(\\Omega). $
%[text] Only first derivatives of $u$ and $v$ appear. The clamped condition $u = 0$ is essential and belongs to $H^1\_0(\\Omega)$, while $\\partial\_n u = 0$ is natural: it was used in the derivation and imposes no condition on $\\sigma$.
%[text] ### Discretization
%[text] The Hellan-Herrmann-Johnson space of degree $k\\ge 0$ (`HHJk(degree,mesh)`)
%[text]{"align":"center"} $ \\Sigma\_h = \\{\\tau : \\tau|\_T\\in P\_k(T;\\mathbb{S}) \\textrm{ for all } T\\in\\mathcal{T},\\ n\_E\\cdot\\tau n\_E \\textrm{ continuous across all } E\\in\\mathcal{E}(\\Omega)\\} $
%[text] consists of the symmetric matrix fields of degree $k$ with continuous normal-normal components, stored as $\[\\sigma\_{11}, \\sigma\_{22}, \\sigma\_{12}\]$. Its degrees of freedom are the coefficients of a hierarchical basis with $k+1$ functions per edge, which carry the normal-normal component, and $3k(k+1)/2$ functions per triangle. The deflection is approximated in $V\_h = P\_{k+1}(\\mathcal{T})\\cap H^1\_0(\\Omega)$ (`Pk(degree+1,mesh)`, free nodes only). For the continuous $v\_h\\in V\_h$ the tangential derivative is single valued on the edges and vanishes on the boundary, which the form expresses by the average $\\{\\nabla v\\cdot t\\}$ over all edges `ds`. The coefficient vectors $\\sigma$ and $u$ solve
%[text]{"align":"center"} $ \\left(\\begin{array}{cc} A & B^T \\\\ B & 0 \\end{array}\\right) \\left(\\begin{array}{c} \\sigma \\\\ u \\end{array}\\right) = \\left(\\begin{array}{c} 0 \\\\ -b \\end{array}\\right). $
%[text:table]
%[text] | Matrix | Integral | fem4m form |
%[text] | --- | --- | --- |
%[text] | $A$ | $\\int\_\\Omega \\sigma:\\tau\\,dx$ | `sigma:tau dx` |
%[text] | $B$ | $b(\\sigma,v)$ | `div(sigma)*grad(v) dx - [sigma*nt]*{grad(v)*t} ds` |
%[text] | $b$ | $\\int\_\\Omega f v\\,dx$ | `f(x)*v dx` |
%[text] |  | $\\Vert u-u\_h\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - u )^2 dx` |
%[text] |  | $\\Vert \\nabla(u-u\_h)\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - grad(u) )^2 dx` |
%[text] |  | $\\Vert D^2u-\\sigma\_h\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - sigma )^2 dx` |
%[text:table]
%[text] Here `[sigma*nt]` is the jump of the normal-tangential component $t\\cdot\\sigma n$ across an edge.
%[text] ### Expected convergence
%[text] For $k\\ge 1$ and a smooth solution the method converges with
%[text]{"align":"center"} $ \\|D^2u-\\sigma\_h\\|\_{L^2(\\Omega)} + \\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)} \\le C\\,h^{k+1}, $
%[text] i.e. like $N^{-(k+1)/2}$ on uniformly refined meshes, and the $L^2$ error of the deflection one order faster, like $N^{-(k+2)/2}$. For $k = 0$ the Hessian converges with the order $h$.
%%
%[text] ## Parameters
degree = 1;     % degree k of sigma_h, u_h has degree k+1
nrLevels = 5;   % number of meshes
%%
%[text] ## Problem data and initial mesh
%[text] The coarse mesh of the square has two triangles; one uniform refinement by newest vertex bisection gives the first mesh.
p = problemBiharmonicSquareClampedExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform','nvb');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
nrDoF4lvl = zeros(1,nrLevels);
L2Error4lvl = zeros(1,nrLevels);
H1Error4lvl = zeros(1,nrLevels);
sigmaError4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element spaces: normal-normal continuous symmetric matrices and Lagrange elements
    feSigma = HHJk(degree,mesh);
    feU = Pk(degree+1,mesh);
    nrDoF4lvl(lvl) = feSigma.nrDoF + feU.nrDoF;
    % assemble
    a = BilinearForm('sigma:tau dx',{'sigma','tau'});
    A = a.assemble(feSigma);
    bForm = BilinearForm('div(sigma)*grad(v) dx - [sigma*nt]*{grad(v)*t} ds',{'sigma','v'});
    B = bForm.assemble(feSigma,feU);
    l = LinearForm('f(x)*v dx');
    b = l.assemble(feU,12,'f',@p.f);
    % solve: all sigma dofs are free, u = 0 on the boundary
    S = [A, B'; B, sparse(feU.nrDoF,feU.nrDoF)];
    rhs = [zeros(feSigma.nrDoF,1); -b];
    free = [1:feSigma.nrDoF, feSigma.nrDoF+feU.freeNodes];
    x = zeros(size(S,1),1);
    x(free) = S(free,free)\rhs(free);
    sigma = x(1:feSigma.nrDoF);
    u = x(feSigma.nrDoF+1:end);
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(lvl) = l2.assemble(u,feU,12,'exact',@p.exactU);
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(lvl) = h1.assemble(u,feU,12,'exact',@p.exactGradU);
    l2Sigma = ErrorForm('( exact(x) - sigma )^2 dx','sigma');
    sigmaError4lvl(lvl) = l2Sigma.assemble(sigma,feSigma,12,'exact',@p.exactD2);
    % refine
    mesh = mesh.refine('uniform','nvb');
end
%%
%[text] ## Results
%[text] ### Convergence table
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
H1Rate = [NaN, -diff(log(H1Error4lvl))./diff(log(nrDoF4lvl))];
sigmaRate = [NaN, -diff(log(sigmaError4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',H1Error4lvl',H1Rate',sigmaError4lvl',sigmaRate', ...
    VariableNames=["N" "L2 error" "L2 rate" "H1 error" "H1 rate" "sigma error" "sigma rate"])
%[text] For $k = 1$ the rates approach $(k+1)/2 = 1$ for the gradient of the deflection and for the Hessian, and $(k+2)/2 = 3/2$ for the $L^2$ error of the deflection.
%%
%[text] ### Discrete solution
clf
plot(feU,u);
view([40,35]);
colorbar;
shading interp;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(u-u_h)||_{L^2}');
plotConvergence(nrDoF4lvl,sigmaError4lvl,'||D^2u-\sigma_h||_{L^2}');
hold off
%[text] For the default degree the rates are lower than for the quintic Argyris element of `exBiharmonicArgyris`, but the elements are simpler: continuous Lagrange elements and normal-normal continuous matrix fields of low degree. The rates grow with the degree.
%%
%[text] ## Possibilities for extensions
%[text] - Set `degree = 2`: the rates grow to $N^{-3/2}$ for the gradient and the Hessian and to $N^{-2}$ for the deflection.
%[text] - Set `degree = 0`, the lowest-order method with piecewise constant $\\sigma\_h$ and $P\_1$ deflections: the Hessian converges like $N^{-1/2}$, the gradient approaches the same rate only slowly, and the deflection converges like $N^{-1}$.
%[text] - Compare $\\sigma\_h$ with the piecewise Hessian of $u\_h$, `ErrorForm('( exact(x) - D2(u) )^2 dx')` assembled for `u` and `feU`: for $k = 1$ it converges only like $N^{-1/2}$, one order in $h$ less than $\\sigma\_h$.
%[text] - Use `problemBiharmonicLshapeClampedExact`, whose solution behaves like $r^{1.5445}$ at the reentrant corner of the L-shaped domain: the Hessian is unbounded and uniform refinement converges much slower. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
