%[text] # Interior penalty DG method for the Stokes problem
%[text] Discretizes the Stokes equations with discontinuous velocities and pressures of the same polynomial degree: the velocity part is the symmetric interior penalty method of `exPoissonDG` for every component, the inhomogeneous Dirichlet data enter by Nitsche's method, a penalty on the pressure jumps stabilizes the equal-order pair, and a Lagrange multiplier fixes the mean value of the pressure.
%[text] **Builds on:** `exPoissonDG`, `exStokesTaylorHood`. **Next:** `exTransportDG`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] On the unit square $\\Omega = (0,1)^2$, given the body force $f$ and the boundary velocity $u\_D$, find the velocity $u$ and the pressure $p$ with
%[text]{"align":"center"} $ -\\Delta u + \\nabla p = f, \\qquad \\mathrm{div}\\,u = 0 \\quad \\textrm{in } \\Omega, \\qquad u = u\_D \\quad \\textrm{on } \\partial\\Omega. $
%[text] The pressure is determined only up to a constant; the mean value condition $\\int\_\\Omega p\\,dx = 0$ fixes it. The problem object `problemStokesExactSquareSinus` with the viscosity $\\nu = 1$ provides $f$ and $u\_D$ of the exact solution
%[text]{"align":"center"} $ u = \\left(\\begin{array}{c} s \\\\ s \\end{array}\\right), \\quad s = \\sin(\\pi(x\_1-x\_2)), \\qquad p = -2\\pi\\cos(\\pi(x\_1-x\_2)) + \\frac{8}{\\pi}, $
%[text] whose velocity does not vanish on the boundary and whose pressure has the mean value zero.
%[text] ### Discretization
%[text] The notation for jumps $\[\\cdot\]$, averages $\\{\\cdot\\}$, normals $n\_E$ and the edges $\\mathcal{E}$, $\\mathcal{E}(\\Omega)$ is that of `exPoissonDG`; all boundary edges are Dirichlet edges. Velocity and pressure are discontinuous piecewise polynomials of the same degree $k$, $u\_h\\in V\_h = P\_k(\\mathcal{T})^2$ (`DGPkVector(degree,mesh)`) and $p\_h\\in Q\_h = P\_k(\\mathcal{T})$ (`DGPk(degree,mesh)`). The discrete problem seeks $(u\_h,p\_h)$ and a multiplier $\\lambda\\in\\mathbb{R}$ with
%[text]{"align":"center"} $ a\_h(u\_h,v) + b\_h(v,p\_h) = \\ell\_h(v), $
%[text]{"align":"center"} $ b\_h(u\_h,q) - s\_h(p\_h,q) - \\lambda\\,(1,q) = \\sum\_{E\\in\\mathcal{E}(\\partial\\Omega)} (u\_D\\cdot n, q)\_E, $
%[text]{"align":"center"} $ (p\_h, 1) = 0 $
%[text] for all $v\\in V\_h$ and $q\\in Q\_h$. The parts are:
%[text] - $a\_h$ is the SIPG bilinear form of `exPoissonDG` applied to both velocity components, with the matrix-valued gradient and the penalty parameter $\\alpha$: $a\_h(u,v) = (\\nabla\_h u : \\nabla\_h v)\_{\\mathcal{T}} - \\sum\_{E\\in\\mathcal{E}} \\left( (\\{\\nabla\_h u\\}n\_E, \[v\])\_E + (\[u\], \\{\\nabla\_h v\\}n\_E)\_E - \\frac{\\alpha}{h\_E}(\[u\],\[v\])\_E \\right)$.
%[text] - $b\_h(v,q) = -(\\mathrm{div}\_h v, q)\_{\\mathcal{T}} + \\sum\_{E\\in\\mathcal{E}} (\[v\]\\cdot n\_E, \\{q\\})\_E$ is the DG divergence: elementwise integration by parts of $(\\nabla p, v)$ for a continuous pressure gives exactly these terms. For the exact velocity, $\\mathrm{div}\\,u = 0$ and $\[u\] = 0$ on the interior edges leave the boundary term $(u\_D\\cdot n, q)\_E$ of the right-hand side.
%[text] - Nitsche's method imposes $u = u\_D$: on the boundary edges $a\_h$ contains the jump terms with $\[u\_h\] = u\_h$, and the right-hand side $\\ell\_h(v) = (f,v) + \\sum\_{E\\in\\mathcal{E}(\\partial\\Omega)} \\left(\\frac{\\alpha}{h\_E}(u\_D,v)\_E - (u\_D, \\nabla\_h v\\,n)\_E\\right)$ adds the same terms with $u\_D$ in place of $u\_h$, so that the exact solution satisfies the discrete equations.
%[text] - The pressure jump stabilization $s\_h(p,q) = \\sum\_{E\\in\\mathcal{E}(\\Omega)} h\_E(\[p\],\[q\])\_E$ vanishes for a continuous pressure. The divergence term alone does not control the pressure of the equal-order pair uniformly in $h$; together with $s\_h$ the method is stable.
%[text] - $b\_h(v,1) = 0$ for every $v\\in V\_h$ (the edge terms cancel the divergence), so constant pressures are not seen by the velocity equation. The Lagrange multiplier $\\lambda$ enforces $(p\_h,1) = 0$; testing the second equation with $q = 1$ shows that $\\lambda$ vanishes up to the quadrature error of $\\int\_{\\partial\\Omega} u\_D\\cdot n\\,ds = 0$. \
%[text] With the vector $c$ of the integrals $(1,q)$ of the pressure basis functions the linear system reads
%[text]{"align":"center"} $ \\left(\\begin{array}{ccc} A & B^T & 0 \\\\ -B & S & c \\\\ 0 & c^T & 0 \\end{array}\\right) \\left(\\begin{array}{c} u \\\\ p \\\\ \\lambda \\end{array}\\right) = \\left(\\begin{array}{c} b\_1 \\\\ -b\_2 \\\\ 0 \\end{array}\\right). $
%[text:table]
%[text] | Matrix | Integral | fem4m form |
%[text] | --- | --- | --- |
%[text] | $A$ | $a\_h(u,v)$ | `grad(u):grad(v) dx - [uxn]:{grad(v)} ds - {grad(u)}:[vxn] ds + h^-1*alpha(x)*[uxn]:[vxn] ds` |
%[text] | $B$ | $b\_h(v,q)$ | `-div(v)*q dx + [v*n]{q} ds` |
%[text] | $S$ | $s\_h(p,q)$ | `h^1*[p][q] dsi` |
%[text] | $c$ | $(1,q)$ | `v dx` |
%[text] | $b\_1$ | $\\ell\_h(v)$ | `f(x)*v dx` and `h^-1*alphaud(x)*v dsd - ud(x)*grad(v)n dsd` |
%[text] | $b\_2$ | $\\sum\_{E\\in\\mathcal{E}(\\partial\\Omega)} (u\_D\\cdot n, q)\_E$ | `ud(x)*vn dsd` |
%[text:table]
%[text] For a vector field, `[uxn]` is the matrix $\[u\]\\otimes n\_E$, so that `[uxn]:{grad(v)}` is $\[u\]\\cdot\\{\\nabla\_h v\\}n\_E$ and `[uxn]:[vxn]` is $\[u\]\\cdot\[v\]$; `[v*n]` is the jump of the normal component. The function `alphaud(x)` is $\\alpha u\_D$, `ud(x)` is $u\_D$.
%[text] ### Expected convergence
%[text] For a smooth solution the equal-order method with pressure stabilization converges with
%[text]{"align":"center"} $ \\|u-u\_h\\|\_{DG} + \\|p-p\_h\\|\_{L^2(\\Omega)} = O(h^k), \\qquad \\|u-u\_h\\|\_{L^2(\\Omega)} = O(h^{k+1}), $
%[text] with the DG norm $\\|v\\|\_{DG}^2 = \\|\\nabla\_h v\\|^2\_{L^2(\\Omega)} + \\sum\_{E\\in\\mathcal{E}} h\_E^{-1}\\|\[v\]\\|^2\_{L^2(E)}$ of `exPoissonDG`, where the jump on a boundary edge is $u\_D - u\_h$. In the number $N$ of unknowns the rates are $k/2$ and $(k+1)/2$.
%%
%[text] ## Parameters
degree = 1;                     % polynomial degree k of velocity and pressure
penalty = 2*degree*(degree+1);  % penalty parameter alpha
nrLevels = 5;                   % number of meshes
%%
%[text] ## Problem data and initial mesh
p = problemStokesExactSquareSinus();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
nrDoF4lvl = zeros(1,nrLevels);
uL2Error4lvl = zeros(1,nrLevels);
uDGError4lvl = zeros(1,nrLevels);
pL2Error4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element spaces
    feU = DGPkVector(degree,mesh);
    feP = DGPk(degree,mesh);
    nrDoF4lvl(lvl) = feU.nrDoF + feP.nrDoF;
    % assemble
    a = BilinearForm(['grad(u):grad(v) dx', ...
        ' - [uxn]:{grad(v)} ds - {grad(u)}:[vxn] ds', ...
        ' + h^-1*alpha(x)*[uxn]:[vxn] ds']);
    A = a.assemble(feU,'alpha',penalty);
    bDiv = BilinearForm('-div(v)*q dx + [v*n]{q} ds',{'v','q'});
    B = bDiv.assemble(feU,feP);
    s = BilinearForm('h^1*[p][q] dsi',{'p','q'});
    S = s.assemble(feP);
    c = LinearForm('v dx');
    C = c.assemble(feP);
    O = sparse(1,feU.nrDoF);
    l = LinearForm('f(x)*v dx');
    b1 = l.assemble(feU,'f',@p.f);
    lNitsche = LinearForm('h^-1*alphaud(x)*v dsd - ud(x)*grad(v)n dsd');
    b1 = b1 + lNitsche.assemble(feU,2*degree+2,'alphaud',@(x)penalty*p.u4Db(x),'ud',@p.u4Db);
    lDiv = LinearForm('ud(x)*vn dsd');
    b2 = lDiv.assemble(feP,2*degree+2,'ud',@p.u4Db);
    % solve
    x = [A, B', O'; -B, S, C; O, C', 0]\[b1; -b2; 0];
    xU = x(1:feU.nrDoF);
    xP = x(feU.nrDoF+(1:feP.nrDoF));
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    uL2Error4lvl(lvl) = l2.assemble(xU,feU,6,'exact',@p.exactU);
    pL2Error4lvl(lvl) = l2.assemble(xP,feP,6,'exact',@p.exactP);
    dgNorm = ErrorForm('( exact(x) - grad(u) )^2 dx + h^-1*([u])^2 dsi + h^-1*( ud(x) - u )^2 dsd');
    uDGError4lvl(lvl) = dgNorm.assemble(xU,feU,6,'exact',@p.exactGradU,'ud',@p.u4Db);
    % refine
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Convergence table
uL2Rate = [NaN, -diff(log(uL2Error4lvl))./diff(log(nrDoF4lvl))];
uDGRate = [NaN, -diff(log(uDGError4lvl))./diff(log(nrDoF4lvl))];
pL2Rate = [NaN, -diff(log(pL2Error4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',uL2Error4lvl',uL2Rate',uDGError4lvl',uDGRate',pL2Error4lvl',pL2Rate', ...
    VariableNames=["N" "u L2 error" "u L2 rate" "u DG error" "u DG rate" "p L2 error" "p L2 rate"])
%[text] The velocity converges with the rates $(k+1)/2$ in the $L^2$ norm and $k/2$ in the DG norm. The pressure converges faster than the guaranteed rate $k/2$ on the coarse meshes; its rate decreases towards $k/2$ under refinement.
%%
%[text] ### Discrete pressure
clf
plot(feP,xP);
view([40,35]);
colorbar;
%%
%[text] ### Streamlines
clf
plotStreamLines(feU,xU);
%[text] The flow runs along the lines $x\_1 - x\_2 = \\textrm{const}$, in opposite directions on the two sides of the diagonal $x\_1 = x\_2$, where $s$ changes its sign.
%%
%[text] ### Velocity field
clf
plotVelocity(feU,xU);
%[text] The arrows show the direction of $u\_h$, normalized to unit length. Away from the diagonal $x\_1 = x\_2$, where $u$ vanishes, they point along $\\pm(1,1)$.
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,uL2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,uDGError4lvl,'||u-u_h||_{DG}');
plotConvergence(nrDoF4lvl,pL2Error4lvl,'||p-p_h||_{L^2}');
hold off
%%
%[text] ## Possibilities for extensions
%[text] - Set `degree = 2`: the velocity converges with the rates $3/2$ and $1$, the pressure at least with the rate $1$.
%[text] - Drop the stabilization, `S = 0*S;` after its assembly: on these meshes the system stays solvable, but the pressure error is considerably larger.
%[text] - Use the pressure space of one degree less, `feP = DGDk(degree-1,mesh);`, and drop the stabilization: the pair $P\_k(\\mathcal{T})^2\\times P\_{k-1}(\\mathcal{T})$ is stable without $s\_h$, and velocity and pressure converge with the same rates as above.
%[text] - `exStokesTaylorHood` solves the Stokes problem with continuous Taylor-Hood elements, `exPoissonHDG` shows the hybridizable DG approach, which also exists for the Stokes equations. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
