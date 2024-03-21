%[text] # Adaptive Raviart-Thomas mixed finite elements for the Poisson problem
%[text] Solves the Poisson problem on the L-shaped domain in mixed form, with the flux $\\sigma = \\nabla u$ as an unknown of its own, by Raviart-Thomas elements for the flux and discontinuous piecewise polynomials for $u$. In the mixed form the Neumann condition becomes an essential and the Dirichlet condition a natural boundary condition. An error estimator for the flux steers the adaptive refinement towards the reentrant corner.
%[text] **Builds on:** `exPoissonMixedP1P0Unstable`, `exPoissonPkAdaptive`. **Next:** `exElasticityAW`, `exBiharmonicHHJ`, `exStokesTaylorHood`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] The L-shaped domain $\\Omega = (-1,1)^2\\setminus(\[0,1\]\\times\[-1,0\])$ has a reentrant corner at the origin. Its boundary consists of the two edges $\\Gamma\_D$ at the reentrant corner and the remaining edges $\\Gamma\_N$. With the flux $\\sigma = \\nabla u$ the Poisson problem $-\\Delta u = f$ becomes the first-order system
%[text]{"align":"center"} $ \\sigma - \\nabla u = 0, \\qquad \\mathrm{div}\\,\\sigma = -f \\quad \\textrm{in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\Gamma\_D, \\qquad \\sigma\\cdot n = g \\textrm{ on } \\Gamma\_N. $
%[text] The problem object `problemScalarEllipticLshapeExact` provides $f = 0$ and the boundary data of the exact solution, in polar coordinates $(r,\\varphi)$ with $\\varphi\\in\[0,2\\pi)$,
%[text]{"align":"center"} $ u(r,\\varphi) = r^{2/3}\\sin(2\\varphi/3), $
%[text] which vanishes on $\\Gamma\_D$, so $u\_D = 0$, and $g = \\nabla u\\cdot n$. The flux $\\sigma = \\nabla u$ behaves like $r^{-1/3}$ and is unbounded at the origin; it belongs to $H^s(\\Omega)^2$ only for $s \< 2/3$.
%[text] ### Mixed formulation
%[text] The flux lives in $H(\\mathrm{div},\\Omega)$, the vector fields in $L^2(\\Omega)^2$ with a divergence in $L^2(\\Omega)$, whose normal components on the boundary are well defined. Multiplying $\\sigma = \\nabla u$ by a test function $\\tau$ with $\\tau\\cdot n = 0$ on $\\Gamma\_N$ and integrating by parts gives $\\int\_\\Omega \\nabla u\\cdot\\tau\\,dx = -\\int\_\\Omega u\\,\\mathrm{div}\\,\\tau\\,dx + \\int\_{\\Gamma\_D} u\_D\\,\\tau\\cdot n\\,ds$. Find $\\sigma\\in H(\\mathrm{div},\\Omega)$ with $\\sigma\\cdot n = g$ on $\\Gamma\_N$ and $u\\in L^2(\\Omega)$ with
%[text]{"align":"center"} $ \\int\_\\Omega \\sigma\\cdot\\tau\\,dx + \\int\_\\Omega u\\,\\mathrm{div}\\,\\tau\\,dx = \\int\_{\\Gamma\_D} u\_D\\,\\tau\\cdot n\\,ds \\quad \\textrm{for all } \\tau\\in H\_N(\\mathrm{div},\\Omega) = \\{\\tau\\in H(\\mathrm{div},\\Omega) : \\tau\\cdot n = 0 \\textrm{ on } \\Gamma\_N\\}, $
%[text]{"align":"center"} $ \\int\_\\Omega \\mathrm{div}\\,\\sigma\\,v\\,dx = -\\int\_\\Omega f v\\,dx \\quad \\textrm{for all } v\\in L^2(\\Omega). $
%[text] The roles of the boundary conditions are exchanged with respect to the formulation of `exPoissonPkAdaptive`: the Dirichlet data enter the right-hand side (natural boundary condition), the Neumann data restrict the space of the flux (essential boundary condition). The system is a saddle point problem; it has a unique solution because the divergence maps $H\_N(\\mathrm{div},\\Omega)$ onto $L^2(\\Omega)$ (inf-sup condition).
%[text] ### Discretization
%[text] On a triangle $T$ the Raviart-Thomas space $RT\_{k-1}(T) = P\_{k-1}(T)^2 + x\\,P\_{k-1}(T)$ consists of polynomial vector fields of degree $k$ whose normal components on the edges have degree $k-1$ and whose divergence lies in $P\_{k-1}(T)$. The global space $RT\_{k-1}(\\mathcal{T})$ of the vector fields with $\\sigma\_h|\_T\\in RT\_{k-1}(T)$ and continuous normal components across the edges is a subspace of $H(\\mathrm{div},\\Omega)$. fem4m counts the polynomial degree $k$: `RTk(k,mesh)` is $RT\_{k-1}(\\mathcal{T})$ with $k$ degrees of freedom per edge, and `RTk(1,mesh)` is the lowest-order space $RT\_0$. The discrete $u\_h$ lies in $P\_{k-1}(\\mathcal{T})$, the discontinuous piecewise polynomials of degree $k-1$ (`DGDk(k-1,mesh)`). Since $\\mathrm{div}\\,RT\_{k-1}(\\mathcal{T}) = P\_{k-1}(\\mathcal{T})$, the discrete inf-sup condition holds with a constant independent of the mesh; `exPoissonMixedP1P0Unstable` shows what happens without it.
%[text] The essential boundary condition fixes the degrees of freedom of the Neumann edges (`feSigma.fixedNodes`) to the $L^2$ projection of $g$ onto the normal components, `feSigma.projectNbBoundary(@p.g)`. With the coefficient vectors $\\sigma$ and $u$ the discrete problem is the linear system
%[text]{"align":"center"} $ \\left(\\begin{array}{cc} A & B^T \\\\ B & 0 \\end{array}\\right) \\left(\\begin{array}{c} \\sigma \\\\ u \\end{array}\\right) = \\left(\\begin{array}{c} b\_D \\\\ -b \\end{array}\\right) $
%[text] in the free degrees of freedom, after the fixed ones have been moved to the right-hand side.
%[text:table]
%[text] | Matrix | Integral | fem4m form |
%[text] | --- | --- | --- |
%[text] | $A$ | $\\int\_\\Omega \\sigma\\cdot\\tau\\,dx$ | `sigma*tau dx` |
%[text] | $B$ | $\\int\_\\Omega \\mathrm{div}\\,\\sigma\\,v\\,dx$ | `div(sigma)*v dx` |
%[text] | $b\_D$ | $\\int\_{\\Gamma\_D} u\_D\\,\\tau\\cdot n\\,ds$ | `uD(x) tau*n dsd` |
%[text] | $b$ | $\\int\_\\Omega f v\\,dx$ | `f(x) v dx` |
%[text:table]
%[text] ### Error estimator
%[text] The exact flux is a gradient: $\\mathrm{rot}\\,\\sigma = \\partial\_1\\sigma\_2 - \\partial\_2\\sigma\_1 = 0$, its tangential component $\\sigma\\cdot t$ is continuous across the edges and equals the tangential derivative $\\partial u\_D/\\partial s$ on $\\Gamma\_D$, and $\\mathrm{div}\\,\\sigma = -f$. The discrete flux satisfies only $\\mathrm{div}\\,\\sigma\_h = -\\Pi\_{k-1}f$ with the $L^2$ projection $\\Pi\_{k-1}$ onto $P\_{k-1}(\\mathcal{T})$. The estimator measures the violations of the other properties:
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{T\\in\\mathcal{T}} h\_T^2\\left(\\|f + \\mathrm{div}\\,\\sigma\_h\\|^2\_{L^2(T)} + \\|\\mathrm{rot}\\,\\sigma\_h\\|^2\_{L^2(T)}\\right) + \\sum\_{E\\in\\mathcal{E}(\\Omega)} h\_E \\|\[\\sigma\_h\\cdot t\_E\]\\|^2\_{L^2(E)} + \\sum\_{E\\in\\mathcal{E}(\\Gamma\_D)} h\_E \\|\\partial u\_D/\\partial s - \\sigma\_h\\cdot t\\|^2\_{L^2(E)}. $
%[text] The first term is the oscillation $f - \\Pi\_{k-1}f$ of the source, zero for $f = 0$; in the last term $\\partial u\_D/\\partial s = 0$ because $u\_D = 0$.
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2\\Vert f+\\mathrm{div}\\,\\sigma\_h\\Vert^2\_{L^2(T)}$ | `h^2*(f(x)+div(sigma))^2 dx` |
%[text] | $h\_T^2\\Vert \\mathrm{rot}\\,\\sigma\_h\\Vert^2\_{L^2(T)}$ | `h^2*(rot(sigma))^2 dx` |
%[text] | $h\_E\\Vert\[\\sigma\_h\\cdot t\_E\]\\Vert^2\_{L^2(E)}$ | `h^1*([sigma*t])^2 dsi` |
%[text] | $h\_E\\Vert\\sigma\_h\\cdot t\\Vert^2\_{L^2(E)}$, $E\\subset\\Gamma\_D$ | `h^1*(sigma*t)^2 dsd` |
%[text:table]
%[text] `RTk` evaluates the values and the divergence of its functions, but not the derivatives in `rot(sigma)`. The estimator is therefore assembled for the $L^2$ projection of $\\sigma\_h$ onto `DGPkVector(k,mesh)`, the discontinuous functions in $P\_k(T)^2$ on every triangle; since $RT\_{k-1}(T)\\subset P\_k(T)^2$, the projection reproduces $\\sigma\_h$ exactly. The estimator is reliable and efficient up to oscillations of the data: $\\|\\sigma-\\sigma\_h\\|\_{L^2(\\Omega)} \\le C\_{\\textrm{rel}}\\,\\eta$ and $\\eta \\le C\_{\\textrm{eff}}\\,\\|\\sigma-\\sigma\_h\\|\_{L^2(\\Omega)}$. The adaptive algorithm of `exPoissonPkAdaptive` uses it with bulk marking.
%[text] ### Expected convergence
%[text] For a smooth solution the flux and the solution converge with the order $k$: $\\|\\sigma-\\sigma\_h\\|\_{L^2(\\Omega)} + \\|u-u\_h\\|\_{L^2(\\Omega)} = O(h^k) = O(N^{-k/2})$. On uniformly refined meshes the corner singularity limits the flux error to $O(h^{2/3}) = O(N^{-1/3})$ for every $k$; adaptive refinement restores $N^{-k/2}$ for both errors.
%%
%[text] ## Parameters
degree = 3;       % polynomial degree k: flux in RTk(k), u in DGDk(k-1)
minNrDoF = 1e5;   % stop as soon as the number of degrees of freedom exceeds minNrDoF
mark = MarkBulk;  % bulk (Doerfler) marking with theta = 1/2
% mark = MarkUniform;  % uniform refinement, for comparison
%%
%[text] ## Problem data and initial mesh
p = problemScalarEllipticLshapeExact();
mesh = Mesh(p.geom);
plot(mesh);
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the loop stops when the number of degrees of freedom of both spaces exceeds `minNrDoF`.
nrDoF4lvl = [];
L2Error4lvl = [];
sigmaError4lvl = [];
eta4lvl = [];
while true
    % finite element spaces
    feSigma = RTk(degree,mesh);
    feU = DGDk(degree-1,mesh);
    nrDoF4lvl(end+1) = feSigma.nrDoF + feU.nrDoF;
    % assemble
    a = BilinearForm('sigma*tau dx',{'sigma','tau'});
    A = a.assemble(feSigma);
    bDiv = BilinearForm('div(sigma)*v dx',{'sigma','v'});
    B = bDiv.assemble(feSigma,feU);
    lD = LinearForm('uD(x) tau*n dsd',{'sigma','tau'});
    bD = lD.assemble(feSigma,'uD',@p.u4Db);
    l = LinearForm('f(x) v dx');
    b = l.assemble(feU,'f',@p.f);
    % solve
    x = zeros(feSigma.nrDoF + feU.nrDoF,1);
    x(1:feSigma.nrDoF) = feSigma.projectNbBoundary(@p.g);
    r = [bD; -b] - [A; B]*x(1:feSigma.nrDoF);
    S = [A, B'; B, sparse(feU.nrDoF,feU.nrDoF)];
    free = [feSigma.freeNodes, feSigma.nrDoF + feU.freeNodes];
    x(free) = S(free,free)\r(free);
    sigma = x(1:feSigma.nrDoF);
    u = x(feSigma.nrDoF+1:end);
    % estimate
    feSigmaPk = DGPkVector(degree,mesh);
    sigmaPk = feSigmaPk.project(@(pts,pts_ref,parts) feSigma.getValue(pts,pts_ref,parts,sigma));
    eta = EstimatorForm(['h^2*(f(x)+div(sigma))^2 dx + h^2*(rot(sigma))^2 dx', ...
        '+ h^1*([sigma*t])^2 dsi + h^1*(sigma*t)^2 dsd'],'sigma');
    eta4e = eta.assemble(sigmaPk,feSigmaPk,'f',@p.f);
    eta4lvl(end+1) = norm(eta4e);
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(end+1) = l2.assemble(u,feU,max(12,2*(degree+1)),'exact',@p.exactU);
    l2Sigma = ErrorForm('( exact(x) - sigma )^2 dx','sigma');
    sigmaError4lvl(end+1) = l2Sigma.assemble(sigma,feSigma,max(12,2*(degree+1)),'exact',@p.exactSigma);
    % stop
    if nrDoF4lvl(end) >= minNrDoF
        break
    end
    % mark
    mark = mark.mark(feSigma,eta4e);
    % refine
    mesh = mesh.refine(mark.marked4ed);
end
%%
%[text] ## Results
%[text] ### Convergence table
%[text] The rates are the slopes of the errors over $N$ in a log-log scale between two consecutive meshes; the last column is the ratio of the flux error and the estimator.
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
sigmaRate = [NaN, -diff(log(sigmaError4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',sigmaError4lvl',sigmaRate',eta4lvl',(sigmaError4lvl./eta4lvl)', ...
    VariableNames=["N" "u error" "u rate" "sigma error" "sigma rate" "eta" "error/eta"])
%[text] The ratio of the flux error and the estimator stays bounded: the estimator overestimates the error by a roughly constant factor.
%%
%[text] ### Adaptive mesh
clf
plot(mesh);
%[text] The refinement concentrates at the reentrant corner, where the flux is singular.
%%
%[text] ### Discrete solution
%[text] The discrete $u\_h$ is discontinuous, a polynomial of degree $k-1$ on every triangle.
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
plotConvergence(nrDoF4lvl,sigmaError4lvl,'||\sigma-\sigma_h||_{L^2}');
plotConvergence(nrDoF4lvl,eta4lvl,'\eta');
hold off
%[text] Both errors and the estimator decrease with the rate $k/2$ of a smooth solution, despite the singularity of the flux: the least-squares slopes in the legend are close to $k/2$, and the rates in the table approach it after a preasymptotic phase.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkUniform;`: the flux error decays only like $N^{-1/3}$.
%[text] - Set `degree = 1` for the lowest-order pair $RT\_0\\times P\_0$: the adaptive rate is $N^{-1/2}$.
%[text] - Use `problemScalarEllipticSquareExact`, the unit square with $u = 0$ on the whole boundary and a smooth solution: there are no Neumann edges, and for `degree` at most 2 the oscillation term of the quadratic source $f$ no longer vanishes.
%[text] - Replace `RTk(degree,mesh)` by the Brezzi-Douglas-Marini space `BDMk(degree,mesh)` of all vector fields in $P\_k(T)^2$ with continuous normal components: the flux error improves to $O(h^{k+1})$, while $u\_h\\in P\_{k-1}(\\mathcal{T})$ still converges with the order $k$.
%[text] - Raise the polynomial degree instead of refining the mesh: replace the refine step by `degree = degree + 1;` and stop at a moderate degree. The flux error decreases, but the estimator does not, because its weights $h\_T^2$ and $h\_E$ ignore the degree. `exPoissonHPAdaptive` uses an estimator weighted with the degrees and decides on every triangle between both kinds of refinement. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
