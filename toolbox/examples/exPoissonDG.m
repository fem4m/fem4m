%[text] # Symmetric interior penalty DG method for the Poisson problem
%[text] Discretizes the Poisson problem with piecewise polynomials that may jump across the edges of the triangulation: the symmetric interior penalty discontinuous Galerkin (SIPG) method couples the triangles by averages and jumps on the edges and imposes the boundary condition weakly. The first example of the DG family introduces jumps and averages, the edge domains of the form strings and the DG norm.
%[text] **Builds on:** `exPoissonP1`. **Next:** `exPoissonDGAdaptive`, `exPoissonLDG`, `exPoissonMixedDG`, `exStokesDG`, `exTransportDG`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] Let $\\Omega = (0,1)^2$ be the unit square. Given the source $f$, find $u$ with
%[text]{"align":"center"} $ -\\Delta u = f \\quad \\textrm{in } \\Omega, \\qquad u = 0 \\quad \\textrm{on } \\partial\\Omega. $
%[text] The problem object `problemScalarEllipticSquareExact` provides the source $f = 2x\_1(1-x\_1) + 2x\_2(1-x\_2)$ of the exact solution $u = x\_1(1-x\_1)\\,x\_2(1-x\_2)$.
%[text] ### Jumps and averages
%[text] Let $\\mathcal{T}$ be a regular triangulation of $\\Omega$ with the set of edges $\\mathcal{E}$, the interior edges $\\mathcal{E}(\\Omega)$ and the boundary edges $\\mathcal{E}(\\partial\\Omega)$, and let $h\_E$ be the length of the edge $E$. Every edge has a fixed unit normal $n\_E$: on an interior edge it points from the first neighbouring triangle $T\_+$ to the second one $T\_-$, on a boundary edge it is the outer normal. For a piecewise smooth function $v$ the jump and the average on an interior edge are
%[text]{"align":"center"} $ \[v\] = v|\_{T\_+} - v|\_{T\_-}, \\qquad \\{v\\} = \\frac{1}{2}\\left(v|\_{T\_+} + v|\_{T\_-}\\right), $
%[text] and on a boundary edge $\[v\] = \\{v\\} = v$; vector fields jump and average componentwise. $\\nabla\_h$ is the gradient on every triangle, $(v,w)\_{\\mathcal{T}} = \\sum\_{T\\in\\mathcal{T}}\\int\_T vw\\,dx$ and $(v,w)\_E = \\int\_E vw\\,ds$. In the form strings `[v]` and `{v}` are the jump and the average, `n` is the normal $n\_E$ (`mesh.normal4ed`), and `ds`, `dsi`, `dsd`, `dsn` integrate over all edges, the interior edges, the Dirichlet edges and the Neumann edges.
%[text] ### Discretization
%[text] The discrete space $V\_h = P\_k(\\mathcal{T})$ consists of all piecewise polynomials of degree $k$ without any continuity (`DGPk(degree,mesh)`): every Lagrange basis function lives on one triangle only. Multiplying $-\\Delta u = f$ by $v\_h\\in V\_h$ and integrating by parts on every triangle gives
%[text]{"align":"center"} $ (\\nabla u, \\nabla\_h v\_h)\_{\\mathcal{T}} - \\sum\_{T\\in\\mathcal{T}} \\int\_{\\partial T} \\nabla u\\cdot n\_T\\, v\_h\\,ds = (f, v\_h) $
%[text] with the outer normal $n\_T$ of $T$. The normal derivative of $u$ is continuous, so the integrals over the boundaries of the triangles collect to $\\sum\_{E\\in\\mathcal{E}} (\\nabla u\\cdot n\_E, \[v\_h\])\_E$. Replacing $u$ by $u\_h$ and $\\nabla u$ on the edges by the average $\\{\\nabla\_h u\_h\\}$ gives a consistent method: the exact solution satisfies its equations. Two terms that vanish for the exact solution complete the method, because $\[u\] = 0$ on the interior edges and $u = 0$ on $\\partial\\Omega$: the transposed consistency term $-(\[u\_h\], \\{\\nabla\_h v\_h\\}\\cdot n\_E)\_E$ makes the bilinear form symmetric, and the penalty term $\\frac{\\alpha}{h\_E}(\[u\_h\],\[v\_h\])\_E$ with the penalty parameter $\\alpha \> 0$ makes it coercive. The SIPG method seeks $u\_h\\in V\_h$ with $a\_h(u\_h,v\_h) = (f,v\_h)$ for all $v\_h\\in V\_h$, where
%[text]{"align":"center"} $ a\_h(u\_h,v\_h) = (\\nabla\_h u\_h, \\nabla\_h v\_h)\_{\\mathcal{T}} - \\sum\_{E\\in\\mathcal{E}} (\\{\\nabla\_h u\_h\\}\\cdot n\_E, \[v\_h\])\_E - \\sum\_{E\\in\\mathcal{E}} (\[u\_h\], \\{\\nabla\_h v\_h\\}\\cdot n\_E)\_E + \\sum\_{E\\in\\mathcal{E}} \\frac{\\alpha}{h\_E}(\[u\_h\],\[v\_h\])\_E. $
%[text] On the boundary edges the three edge terms impose $u = 0$ weakly (Nitsche's method), so no degree of freedom is fixed and the linear system $Ax = b$ involves all of them. `alpha(x)` is a coefficient of the form; `assemble` receives its value `penalty`.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $(\\nabla\_h u, \\nabla\_h v)\_{\\mathcal{T}}$ | `grad(u)*grad(v) dx` |
%[text] | $-\\sum\_{E\\in\\mathcal{E}} (\[u\], \\{\\nabla\_h v\\}\\cdot n\_E)\_E$ | `- [u]*n.{grad(v)} ds` |
%[text] | $-\\sum\_{E\\in\\mathcal{E}} (\\{\\nabla\_h u\\}\\cdot n\_E, \[v\])\_E$ | `- {grad(u)}.[v]*n ds` |
%[text] | $\\sum\_{E\\in\\mathcal{E}} \\frac{\\alpha}{h\_E}(\[u\],\[v\])\_E$ | `h^-1*alpha(x)*[u][v] ds` |
%[text] | $(f, v)$ | `f(x) v dx` |
%[text:table]
%[text] ### Penalty parameter and DG norm
%[text] The error is measured in the DG norm
%[text]{"align":"center"} $ \\|v\\|\_{DG}^2 = \\|\\nabla\_h v\\|^2\_{L^2(\\Omega)} + \\sum\_{E\\in\\mathcal{E}} h\_E^{-1}\\|\[v\]\\|^2\_{L^2(E)}, $
%[text] the energy norm of the continuous problem plus the jumps of $v$, which vanish for the exact solution; for $v\_h\\in V\_h$ it is a norm, since $\\|v\_h\\|\_{DG} = 0$ forces $v\_h$ to be continuous, piecewise constant and zero on $\\partial\\Omega$. The consistency terms of $a\_h$ have no sign. Polynomials $w$ of degree $k-1$ satisfy the discrete trace inequality $h\_E\\|w\\|^2\_{L^2(E)} \\le C\_T\\,\\frac{k(k+1)}{2}\\|w\\|^2\_{L^2(T)}$ on every edge $E$ of a triangle $T$, with a constant $C\_T$ that depends only on the shape of $T$. Applied to $w = \\nabla v\_h$, together with Young's inequality, it shows that $a\_h$ is coercive in the DG norm if $\\alpha$ exceeds a constant times $k(k+1)$ that depends on the shape regularity of $\\mathcal{T}$. The example uses $\\alpha = 2k(k+1)$. A too small penalty destroys the coercivity, a too large one increases the condition number of $A$.
%[text] The meshes are refined uniformly by newest vertex bisection, `mesh.refine('uniform','nvb')`, instead of the default red-green-blue refinement: every triangle is divided into four by bisecting each of its edges. Both refinements keep the meshes shape regular; newest vertex bisection is also the refinement of the adaptive version `exPoissonDGAdaptive`.
%[text] ### Expected convergence
%[text] Consistency, coercivity and boundedness give the quasi-optimal error estimate in the DG norm, and the symmetry of $a\_h$ allows a duality argument on the convex domain: for $u\\in H^{k+1}(\\Omega)$
%[text]{"align":"center"} $ \\|u-u\_h\\|\_{DG} \\le C\\,h^k\\,|u|\_{H^{k+1}(\\Omega)}, \\qquad \\|u-u\_h\\|\_{L^2(\\Omega)} \\le C\\,h^{k+1}\\,|u|\_{H^{k+1}(\\Omega)}. $
%[text] Every triangle carries $(k+1)(k+2)/2$ unknowns, so $N\\approx h^{-2}$ and the errors decay like $N^{-k/2}$ and $N^{-(k+1)/2}$.
%%
%[text] ## Parameters
degree = 1;                     % polynomial degree k
penalty = 2*degree*(degree+1);  % penalty parameter alpha
nrLevels = 6;                   % number of meshes
%%
%[text] ## Problem data and initial mesh
%[text] The coarse mesh of the square has two triangles; one refinement gives the first mesh.
p = problemScalarEllipticSquareExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform','nvb');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
nrDoF4lvl = zeros(1,nrLevels);
L2Error4lvl = zeros(1,nrLevels);
DGError4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element space
    fe = DGPk(degree,mesh);
    nrDoF4lvl(lvl) = fe.nrDoF;
    % assemble
    a = BilinearForm(['grad(u)*grad(v) dx', ...
        ' - [u]*n.{grad(v)} ds - {grad(u)}.[v]*n ds', ...
        ' + h^-1*alpha(x)*[u][v] ds']);
    A = a.assemble(fe,'alpha',penalty);
    l = LinearForm('f(x) v dx');
    b = l.assemble(fe,'f',@p.f);
    % solve
    x = A\b;
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(lvl) = l2.assemble(x,fe,2*degree+2,'exact',@p.exactU);
    dgNorm = ErrorForm('( exact(x) - grad(u) )^2 dx + h^-1*([u])^2 ds');
    DGError4lvl(lvl) = dgNorm.assemble(x,fe,2*degree+2,'exact',@p.exactGradU);
    % refine
    mesh = mesh.refine('uniform','nvb');
end
%%
%[text] ## Results
%[text] ### Convergence table
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
DGRate = [NaN, -diff(log(DGError4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',DGError4lvl',DGRate', ...
    VariableNames=["N" "L2 error" "L2 rate" "DG error" "DG rate"])
%[text] After a preasymptotic phase on the coarsest meshes, the $L^2$ rates approach $(k+1)/2$ and the rates in the DG norm $k/2$, as predicted.
%%
%[text] ### Discrete solution
clf
plot(fe,x);
view([40,35]);
colorbar;
%[text] The jumps of $u\_h$ across the edges are of the size of the error and not visible on the finest mesh.
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,DGError4lvl,'||u-u_h||_{DG}');
hold off
%[text] The least-squares slopes in the legend include the coarsest meshes and are therefore close to, but not exactly at, the asymptotic rates.
%%
%[text] ## Possibilities for extensions
%[text] - Set `degree = 2` or `3`: the errors decay like $N^{-k/2}$ in the DG norm and $N^{-(k+1)/2}$ in the $L^2$ norm. For $k\\ge 4$ the polynomial solution lies in $V\_h$ and the errors drop to rounding level.
%[text] - Lower the penalty, e.g. `penalty = 1` for `degree = 1`: on these meshes $A$ is positive definite only if $\\alpha$ exceeds a threshold between about $1.1\\,k(k+1)$ and $1.4\\,k(k+1)$ (about $2.8$ for $k = 1$); below it `min(eig(full(A)))` is negative and the errors no longer converge. A large penalty, e.g. `penalty = 20*degree*(degree+1)`, keeps the rates, but enlarges the errors and the condition number of $A$.
%[text] - Flip the sign of the transposed consistency term, `' + [u]*n.{grad(v)} ds'`: the nonsymmetric interior penalty method (NIPG) is coercive for every $\\alpha \> 0$ and converges with the same rate in the DG norm, but it is not adjoint consistent, and for even $k$ its $L^2$ error loses one order.
%[text] - Refine with `mesh.refine('uniform')`, the red refinement: the meshes have the same nodes but other diagonals, and the rates are the same.
%[text] - `exPoissonDGAdaptive` solves the L-shaped domain problem adaptively, `exPoissonLDG` is a DG method that is stable for every penalty $\\alpha \> 0$. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
