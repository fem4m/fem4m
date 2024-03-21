%[text] # Hybridizable discontinuous Galerkin method for the Poisson problem
%[text] Discretizes the Poisson problem in mixed form with discontinuous approximations of the flux and of the solution, which are coupled only through a single-valued approximation of the trace on the edges: the hybridizable discontinuous Galerkin (HDG) method. The element unknowns are eliminated triangle by triangle, so that only the trace unknowns remain in the global system.
%[text] **Builds on:** `exPoissonDG`, `exPoissonMixedDG`. **Next:** `exStokesDG`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] On the unit square $\\Omega = (0,1)^2$ the Poisson problem $-\\Delta u = f$ with $u = 0$ on $\\partial\\Omega$ reads in mixed form, with the flux $\\sigma = -\\nabla u$,
%[text]{"align":"center"} $ \\sigma + \\nabla u = 0, \\qquad \\mathrm{div}\\,\\sigma = f \\quad \\textrm{in } \\Omega, \\qquad u = 0 \\quad \\textrm{on } \\partial\\Omega. $
%[text] The problem object `problemScalarEllipticSquareExact` provides the source of the exact solution $u = x\_1(1-x\_1)\\,x\_2(1-x\_2)$.
%[text] ### Discretization
%[text] On every triangle $T\\in\\mathcal{T}$ the method approximates $\\sigma$ by $\\sigma\_h|\_T\\in P\_k(T)^2$ and $u$ by $u\_h|\_T\\in P\_l(T)$, both discontinuous across the edges (`DGDkVector(degree,mesh)` and `DGDk(degreeU,mesh)`). A third unknown $\\lambda\_h$ approximates the trace of $u$ on the interior edges $\\mathcal{E}(\\Omega)$ by piecewise polynomials of degree $m = \\max(k,l)$ (`DGDkEdge`); on $\\partial\\Omega$ the Dirichlet condition sets $\\lambda\_h = 0$. Integration by parts on every triangle, with the trace $\\lambda\_h$ and the numerical flux
%[text]{"align":"center"} $ \\hat\\sigma\_h\\cdot n\_T = \\sigma\_h\\cdot n\_T + \\tau\\,(u\_h - \\lambda\_h) \\quad \\textrm{on } \\partial T $
%[text] with the stabilization parameter $\\tau\\ge 0$, gives the discrete problem: find $(\\sigma\_h, u\_h, \\lambda\_h)$ with
%[text]{"align":"center"} $ (\\sigma\_h, q)\_{\\mathcal{T}} - (u\_h, \\mathrm{div}\\,q)\_{\\mathcal{T}} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} (\\lambda\_h, \[q\]\\cdot n\_E)\_E = 0, $
%[text]{"align":"center"} $ (\\mathrm{div}\\,\\sigma\_h, v)\_{\\mathcal{T}} + \\sum\_{T\\in\\mathcal{T}} (\\tau\\,(u\_h-\\lambda\_h), v)\_{\\partial T} = (f, v), $
%[text]{"align":"center"} $ \\sum\_{E\\in\\mathcal{E}(\\Omega)} (\[\\sigma\_h\]\\cdot n\_E + 2\\tau\\,(\\{u\_h\\}-\\lambda\_h), \\mu)\_E = 0 $
%[text] for all test functions $q$, $v$ and $\\mu$ of the three spaces. Here $(\\cdot,\\cdot)\_{\\mathcal{T}}$ sums the integrals over the triangles, and $\[\\cdot\]$ and $\\{\\cdot\\}$ are the jump and the average across an edge. The first two equations are the mixed equations on every triangle, the third one makes the normal component of the numerical flux continuous across the edges, so the method is locally conservative.
%[text] ### Linear system
%[text] On an interior edge the stabilization term of the second equation splits into $2\\tau(\\{u\_h\\},\\{v\\}) + \\frac{\\tau}{2}(\[u\_h\],\[v\]) - 2\\tau(\\lambda\_h,\\{v\\})$, on a boundary edge it is $\\tau(u\_h,v)$. With the matrices of the forms below and the coefficient vectors $\\sigma$, $u$, $\\lambda$ the method reads
%[text]{"align":"center"} $ \\left(\\begin{array}{ccc} A & -B^T & C^T \\\\ B & 2D\_1 + \\frac{1}{2} D\_2 + D\_\\partial & -2E^T \\\\ C & 2E & -2F \\end{array}\\right) \\left(\\begin{array}{c} \\sigma \\\\ u \\\\ \\lambda \\end{array}\\right) = \\left(\\begin{array}{c} 0 \\\\ b \\\\ 0 \\end{array}\\right). $
%[text:table]
%[text] | Matrix | Integral | fem4m form |
%[text] | --- | --- | --- |
%[text] | $A$ | $(\\sigma, q)\_{\\mathcal{T}}$ | `sigma*q dx` |
%[text] | $B$ | $(\\mathrm{div}\\,\\sigma, v)\_{\\mathcal{T}}$ | `div(sigma) v dx` |
%[text] | $C$ | $\\sum\_{E\\in\\mathcal{E}(\\Omega)} (\[\\sigma\]\\cdot n\_E, \\mu)\_E$ | `[sigma]*n {mu} dsi` |
%[text] | $D\_1$ | $\\sum\_{E\\in\\mathcal{E}(\\Omega)} (\\tau\\{u\\}, \\{v\\})\_E$ | `t(x)*{u}{v} dsi` |
%[text] | $D\_2$ | $\\sum\_{E\\in\\mathcal{E}(\\Omega)} (\\tau\[u\], \[v\])\_E$ | `t(x)*[u][v] dsi` |
%[text] | $D\_\\partial$ | $\\sum\_{E\\subset\\Gamma\_D} (\\tau u, v)\_E$ | `t(x)*{u}{v} dsd` |
%[text] | $E$ | $\\sum\_{E\\in\\mathcal{E}(\\Omega)} (\\tau\\{u\\}, \\mu)\_E$ | `t(x)*{u}{mu} dsi` |
%[text] | $F$ | $\\sum\_{E\\in\\mathcal{E}(\\Omega)} (\\tau\\lambda, \\mu)\_E$ | `t(x)*{lambda}{mu} dsi` |
%[text] | $b$ | $(f, v)$ | `f(x) v dx` |
%[text:table]
%[text] ### Static condensation
%[text] The upper left block $A\_{11}$ of the first two block rows does not couple different triangles: $A$ and $B$ act on one triangle, and the combination $2D\_1 + \\frac{1}{2} D\_2 + D\_\\partial$ is the elementwise sum $\\sum\_T (\\tau u\_h, v)\_{\\partial T}$. Eliminating $\\sigma$ and $u$ triangle by triangle leaves the Schur complement system for $\\lambda$ alone, which couples only the edges of a common triangle:
%[text]{"align":"center"} $ \\left(A\_{21}A\_{11}^{-1}A\_{12} - A\_{22}\\right)\\lambda = A\_{21}A\_{11}^{-1}r, \\qquad \\left(\\begin{array}{c} \\sigma \\\\ u \\end{array}\\right) = A\_{11}^{-1}\\left(r - A\_{12}\\lambda\\right), $
%[text] with the blocks $A\_{12}$, $A\_{21}$, $A\_{22} = -2F$ of the system above and $r = (0, b)$.
%[text] ### Expected convergence
%[text] For $k = l$ and $\\tau$ of order one the flux and the solution converge with the order $k+1$:
%[text]{"align":"center"} $ \\|\\sigma - \\sigma\_h\\|\_{L^2(\\Omega)} + \\|u - u\_h\\|\_{L^2(\\Omega)} \\le C\\,h^{k+1}, $
%[text] i.e. like $N^{-(k+1)/2}$ in the number of unknowns $N$.
%%
%[text] ## Parameters
degree = 0;    % degree k of the flux sigma_h
degreeU = 0;   % degree l of u_h
tau = 1;       % stabilization parameter
% degree = 1; degreeU = 1; tau = 1;  % k = l = 1: order 2 for both errors
% degree = 1; degreeU = 0; tau = 0;  % no stabilization: the hybridized BDM_1 x P_0 mixed method
nrLevels = 6;  % number of meshes
%%
%[text] ## Problem data and initial mesh
p = problemScalarEllipticSquareExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
nrDoF4lvl = zeros(1,nrLevels);
L2Error4lvl = zeros(1,nrLevels);
sigmaError4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element spaces
    feSigma = DGDkVector(degree,mesh);
    feU = DGDk(degreeU,mesh);
    feLambda = DGDkEdge(max(degree,degreeU),mesh);
    free = feLambda.freeNodes;  % the interior edges, lambda = 0 on the Dirichlet boundary
    nrDoF4lvl(lvl) = feSigma.nrDoF + feU.nrDoF + numel(free);
    % assemble
    a = BilinearForm('sigma*q dx',{'sigma','q'});
    A = a.assemble(feSigma);
    bDiv = BilinearForm('div(sigma) v dx',{'sigma','v'});
    B = bDiv.assemble(feSigma,feU);
    c = BilinearForm('[sigma]*n {mu} dsi',{'sigma','mu'});
    C = c.assemble(feSigma,feLambda);
    d1 = BilinearForm('t(x)*{u}{v} dsi');
    D1 = d1.assemble(feU,feU,'t',tau);
    d2 = BilinearForm('t(x)*[u][v] dsi');
    D2 = d2.assemble(feU,feU,'t',tau);
    dBd = BilinearForm('t(x)*{u}{v} dsd');
    DBd = dBd.assemble(feU,feU,'t',tau);
    e = BilinearForm('t(x)*{u}{mu} dsi',{'u','mu'});
    E = e.assemble(feU,feLambda,'t',tau);
    f = BilinearForm('t(x)*{lambda}{mu} dsi',{'lambda','mu'});
    F = f.assemble(feLambda,feLambda,'t',tau);
    l = LinearForm('f(x) v dx');
    b = l.assemble(feU,'f',@p.f);
    % solve: static condensation, then the element unknowns
    A11 = [A, -B'; B, 2*D1 + D2/2 + DBd];
    A12 = [C(free,:)'; -2*E(free,:)'];
    A21 = [C(free,:), 2*E(free,:)];
    A22 = -2*F(free,free);
    r = [zeros(feSigma.nrDoF,1); b];
    lambda = (A21*(A11\A12) - A22)\(A21*(A11\r));
    x = A11\(r - A12*lambda);
    sigma = x(1:feSigma.nrDoF);
    u = x(feSigma.nrDoF+1:end);
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(lvl) = l2.assemble(u,feU,2*degreeU+2,'exact',@p.exactU);
    l2Sigma = ErrorForm('( exact(x) + sigma )^2 dx','sigma');
    sigmaError4lvl(lvl) = l2Sigma.assemble(sigma,feSigma,2*degree+2,'exact',@p.exactGradU);
    % refine
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Convergence table
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
sigmaRate = [NaN, -diff(log(sigmaError4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',sigmaError4lvl',sigmaRate', ...
    VariableNames=["N" "u error" "u rate" "sigma error" "sigma rate"])
%[text] Both rates approach $(k+1)/2$ with respect to $N$, i.e. the order $k+1$ in $h$.
%%
%[text] ### Discrete solution
clf
plot(feU,u);
view([40,35]);
colorbar;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,sigmaError4lvl,'||\sigma-\sigma_h||_{L^2}');
hold off
%%
%[text] ### Condensed system
%[text] The Schur complement of the finest mesh couples only the edges of a common triangle.
clf
spy(A21*(A11\A12) - A22);
%%
%[text] ## Possibilities for extensions
%[text] - Use `degree = 1; degreeU = 1;`: both errors converge with the order 2 in $h$. Raising only `degreeU` gains nothing here, the flux approximation $\\sigma\_h\\in P\_0$ limits the accuracy.
%[text] - Use `degree = 1; degreeU = 0; tau = 0;`: without stabilization the method is the hybridized Brezzi-Douglas-Marini mixed method $BDM\_1\\times P\_0$; compare it with `exPoissonRTAdaptive`.
%[text] - Vary `tau` over several orders of magnitude and observe its influence on the two errors.
%[text] - Solve the whole system at once, `[A11, A12; A21, A22]\[r; zeros(numel(free),1)]`, and compare its size and the solution with the condensed one.
%[text] - `exPoissonLDG` and `exPoissonMixedDG` discretize the same mixed formulation with other numerical fluxes and without the trace unknown. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
