%[text] # Local discontinuous Galerkin method for the Poisson problem
%[text] Discretizes the Poisson problem with the local discontinuous Galerkin (LDG) method: the gradient of the discrete solution is replaced by a discrete gradient that includes a lifting of its jumps. The resulting method looks like the conforming one with an additional jump penalty, and in contrast to the interior penalty method of `exPoissonDG` it is stable for every penalty parameter.
%[text] **Builds on:** `exPoissonDG`. **Next:** `exPoissonMixedDG`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] On the unit square $\\Omega = (0,1)^2$ find $u$ with $-\\Delta u = f$ in $\\Omega$ and $u = 0$ on $\\partial\\Omega$. The problem object `problemScalarEllipticSquareExact` provides the source $f = 2x\_1(1-x\_1) + 2x\_2(1-x\_2)$ of the exact solution $u = x\_1(1-x\_1)\\,x\_2(1-x\_2)$. With the flux $\\sigma = -\\nabla u$ the problem is the first-order system
%[text]{"align":"center"} $ \\sigma + \\nabla u = 0, \\qquad \\mathrm{div}\\,\\sigma = f \\quad \\textrm{in } \\Omega, \\qquad u = 0 \\quad \\textrm{on } \\partial\\Omega. $
%[text] ### Numerical fluxes
%[text] The notation for jumps $\[\\cdot\]$, averages $\\{\\cdot\\}$, normals $n\_E$ and the edges $\\mathcal{E}$, $\\mathcal{E}(\\Omega)$ is that of `exPoissonDG`. Let $V\_h = P\_k(\\mathcal{T})$ and $\\Sigma\_h = P\_{k-1}(\\mathcal{T})^2$ be discontinuous piecewise polynomials (`DGDk(degree,mesh)` and `DGDkVector(degree-1,mesh)`). Multiplying the two equations by $\\tau\\in\\Sigma\_h$ and $v\\in V\_h$, integrating by parts on every triangle $T$ and replacing the traces of $u$ and $\\sigma$ on $\\partial T$ by single-valued numerical fluxes $\\hat u\_h$ and $\\hat\\sigma\_h$ gives
%[text]{"align":"center"} $ (\\sigma\_h,\\tau)\_T - (u\_h,\\mathrm{div}\\,\\tau)\_T + (\\hat u\_h, \\tau\\cdot n\_T)\_{\\partial T} = 0, \\qquad -(\\sigma\_h,\\nabla v)\_T + (\\hat\\sigma\_h\\cdot n\_T, v)\_{\\partial T} = (f,v)\_T. $
%[text] The LDG method takes the average $\\hat u\_h = \\{u\_h\\}$ on the interior edges and the boundary value $\\hat u\_h = 0$ on $\\partial\\Omega$, and $\\hat\\sigma\_h = \\{\\sigma\_h\\} + \\frac{\\alpha}{h\_E}\[u\_h\]\\,n\_E$ on all edges with a penalty parameter $\\alpha \> 0$ (on a boundary edge $\\{\\sigma\_h\\} = \\sigma\_h$ and $\[u\_h\] = u\_h$). Both fluxes are consistent: they equal $u$ and $\\sigma$ for the exact solution.
%[text] ### Lifting and discrete gradient
%[text] Since $\\hat u\_h$ does not depend on $\\sigma\_h$, the first equation determines $\\sigma\_h$ triangle by triangle from $u\_h$; this local elimination gives the method its name. Integrating the first equation by parts once more shows $\\sigma\_h = -G\_h(u\_h)$ with the discrete gradient
%[text]{"align":"center"} $ G\_h(v) = \\nabla\_h v - R\_h(\[v\]), \\qquad (R\_h(\\varphi), \\tau)\_{\\mathcal{T}} = \\sum\_{E\\in\\mathcal{E}} (\\varphi, \\{\\tau\\}\\cdot n\_E)\_E \\quad \\textrm{for all } \\tau\\in\\Sigma\_h. $
%[text] The lifting operator $R\_h$ maps the jumps on the edges to a piecewise polynomial vector field in $\\Sigma\_h$; on a triangle $T$, $R\_h(\[v\])$ depends only on the jumps on the three edges of $T$. For a continuous $v$ that vanishes on $\\partial\\Omega$ all jumps vanish and $G\_h(v) = \\nabla v$. Inserting $\\sigma\_h = -G\_h(u\_h)$ into the second equation, whose edge terms combine to $(\\sigma\_h, R\_h(\[v\]))\_{\\mathcal{T}}$ plus the penalty, leaves the primal problem: find $u\_h\\in V\_h$ with
%[text]{"align":"center"} $ (G\_h(u\_h), G\_h(v\_h))\_{\\mathcal{T}} + \\sum\_{E\\in\\mathcal{E}} \\frac{\\alpha}{h\_E}(\[u\_h\],\[v\_h\])\_E = (f, v\_h) \\quad \\textrm{for all } v\_h\\in V\_h. $
%[text] Expanding the first term gives the bilinear form of the SIPG method of `exPoissonDG` plus the term $(R\_h(\[u\_h\]), R\_h(\[v\_h\]))\_{\\mathcal{T}}$. The left-hand side is positive for every $v\_h \\ne 0$ and every $\\alpha \> 0$: if all jumps of $v\_h$ vanish, then $v\_h$ is continuous, zero on $\\partial\\Omega$ and $G\_h(v\_h) = \\nabla v\_h$. No lower bound on the penalty is needed, unlike for SIPG.
%[text] ### Linear system
%[text] With the mass matrix $M$ of $\\Sigma\_h$, which is block diagonal with one block per triangle, and the matrix $B$ of $(G\_h(u),\\tau)\_{\\mathcal{T}}$, the coefficients of the discrete gradient are $G = M^{-1}B$, and the coefficient vector $x$ of $u\_h$ solves
%[text]{"align":"center"} $ (G^T M G + S)\\,x = b. $
%[text:table]
%[text] | Matrix | Integral | fem4m form |
%[text] | --- | --- | --- |
%[text] | $M$ | $(\\sigma,\\tau)\_{\\mathcal{T}}$ | `sigma*tau dx` |
%[text] | $B$ | $(\\nabla\_h u,\\tau)\_{\\mathcal{T}} - \\sum\_{E\\in\\mathcal{E}} (\[u\], \\{\\tau\\}\\cdot n\_E)\_E$ | `grad(u)*tau dx - [u]*n{tau} ds` |
%[text] | $S$ | $\\sum\_{E\\in\\mathcal{E}} \\frac{\\alpha}{h\_E}(\[u\],\[v\])\_E$ | `h^-1*alpha(x)*[u][v] ds` |
%[text] | $b$ | $(f,v)$ | `f(x) v dx` |
%[text:table]
%[text] Besides $\\|u-u\_h\\|\_{L^2(\\Omega)}$ and the DG norm $\\|u-u\_h\\|\_{DG}$ of `exPoissonDG`, the example measures the error of the discrete gradient, $\\|\\nabla u - G\_h(u\_h)\\|\_{L^2(\\Omega)} = \\|\\sigma-\\sigma\_h\\|\_{L^2(\\Omega)}$.
%[text] ### Expected convergence
%[text] For a smooth solution and a penalty of the form $\\alpha/h\_E$ the a priori analysis of the LDG method gives the following orders, and the computations below confirm them for the flux space $\\Sigma\_h = P\_{k-1}(\\mathcal{T})^2$ of this example:
%[text]{"align":"center"} $ \\|u-u\_h\\|\_{DG} + \\|\\nabla u - G\_h(u\_h)\\|\_{L^2(\\Omega)} = O(h^k), \\qquad \\|u-u\_h\\|\_{L^2(\\Omega)} = O(h^{k+1}), $
%[text] i.e. like $N^{-k/2}$ and $N^{-(k+1)/2}$ in the number $N$ of unknowns of $u\_h$.
%%
%[text] ## Parameters
degree = 1;    % polynomial degree k >= 1 of u_h
penalty = 1;   % penalty parameter alpha > 0
nrLevels = 5;  % number of meshes
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
DGError4lvl = zeros(1,nrLevels);
gradError4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element spaces
    fe = DGDk(degree,mesh);
    feSigma = DGDkVector(degree-1,mesh);
    nrDoF4lvl(lvl) = fe.nrDoF;
    % assemble: the discrete gradient G = M\B
    m = BilinearForm('sigma*tau dx',{'sigma','tau'});
    M = m.assemble(feSigma);
    bGrad = BilinearForm('grad(u)*tau dx - [u]*n{tau} ds',{'u','tau'});
    B = bGrad.assemble(fe,feSigma);
    G = M\B;
    s = BilinearForm('h^-1*alpha(x)*[u][v] ds');
    S = s.assemble(fe,'alpha',penalty);
    l = LinearForm('f(x) v dx');
    b = l.assemble(fe,'f',@p.f);
    % solve
    x = (G'*M*G + S)\b;
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(lvl) = l2.assemble(x,fe,2*degree+2,'exact',@p.exactU);
    dgNorm = ErrorForm('( exact(x) - grad(u) )^2 dx + h^-1*([u])^2 ds');
    DGError4lvl(lvl) = dgNorm.assemble(x,fe,2*degree+2,'exact',@p.exactGradU);
    l2Grad = ErrorForm('( exact(x) - sigma )^2 dx','sigma');
    gradError4lvl(lvl) = l2Grad.assemble(G*x,feSigma,2*degree+2,'exact',@p.exactGradU);
    % refine
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Convergence table
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
DGRate = [NaN, -diff(log(DGError4lvl))./diff(log(nrDoF4lvl))];
gradRate = [NaN, -diff(log(gradError4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',DGError4lvl',DGRate',gradError4lvl',gradRate', ...
    VariableNames=["N" "L2 error" "L2 rate" "DG error" "DG rate" "G error" "G rate"])
%[text] The $L^2$ rates approach $(k+1)/2$, the rates of the DG norm and of the discrete gradient $k/2$.
%%
%[text] ### Discrete solution
clf
plot(fe,x);
view([40,35]);
colorbar;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,DGError4lvl,'||u-u_h||_{DG}');
plotConvergence(nrDoF4lvl,gradError4lvl,'||\nabla u-G_h(u_h)||_{L^2}');
hold off
%[text] The error of the discrete gradient and the error in the DG norm are of the same size and decay with the same rate.
%%
%[text] ## Possibilities for extensions
%[text] - Make the penalty small, e.g. `penalty = 0.01`: the method stays stable and the error of the discrete gradient hardly changes, but the jumps of $u\_h$ grow, and with them the errors of $u\_h$ in the DG and in the $L^2$ norm. The SIPG method of `exPoissonDG` loses its coercivity for such a penalty.
%[text] - Set `degree = 2`: the $L^2$ rate becomes $3/2$, the rates of the DG norm and of the discrete gradient $1$.
%[text] - Use the flux space of the original LDG method, `feSigma = DGDkVector(degree,mesh)`, i.e. $\\Sigma\_h = P\_k(\\mathcal{T})^2$: the discrete gradient changes, the rates stay the same.
%[text] - Compare `nnz(G'*M*G + S)` with the number of nonzero entries of the SIPG matrix of `exPoissonDG` on the same mesh: the discrete gradient on a triangle involves its neighbours, so the LDG matrix couples a triangle also with the neighbours of its neighbours and has more nonzero entries. This wider stencil is the price for the stability for every penalty.
%[text] - `exPoissonMixedDG` keeps $\\sigma\_h\\in P\_k(\\mathcal{T})^2$ as an unknown and varies both numerical fluxes; the LDG fluxes are its special case $C\_{11} = \\alpha/h\_E$, $C\_{12} = 0$, $C\_{22} = 0$. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
