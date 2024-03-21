%[text] # Adaptive interior penalty DG method for the Poisson problem
%[text] Solves the Poisson problem with mixed boundary conditions on the L-shaped domain by the symmetric interior penalty DG method of `exPoissonDG` on adaptively refined meshes. The residual-based error estimator contains, besides the residuals of conforming methods, the jumps of the discrete solution, which measure its distance to the continuous functions.
%[text] **Builds on:** `exPoissonDG`, `exPoissonPkAdaptive`. **Next:** `exPoissonLDG`, `exPoissonMixedDG`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] The L-shaped domain $\\Omega = (-1,1)^2\\setminus(\[0,1\]\\times\[-1,0\])$ has a reentrant corner at the origin. Its boundary consists of the two edges $\\Gamma\_D$ at the reentrant corner and the remaining edges $\\Gamma\_N$. Given $f$, $u\_D$ and $g$, find $u$ with
%[text]{"align":"center"} $ -\\Delta u = f \\textrm{ in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\Gamma\_D, \\qquad \\nabla u\\cdot n = g \\textrm{ on } \\Gamma\_N. $
%[text] The problem object `problemScalarEllipticLshapeExact` provides $f = 0$ and the Neumann data $g = \\nabla u\\cdot n$ of the exact solution, in polar coordinates $(r,\\varphi)$ with $\\varphi\\in\[0,2\\pi)$,
%[text]{"align":"center"} $ u(r,\\varphi) = r^{2/3}\\sin(2\\varphi/3). $
%[text] It vanishes on $\\Gamma\_D$, so $u\_D = 0$. Its gradient behaves like $r^{-1/3}$ and is unbounded at the origin: $u$ belongs to $H^{1+s}(\\Omega)$ only for $s \< 2/3$.
%[text] ### Discretization
%[text] The notation for jumps $\[\\cdot\]$, averages $\\{\\cdot\\}$ and normals $n\_E$ is that of `exPoissonDG`. The edges split into the interior edges $\\mathcal{E}(\\Omega)$, the Dirichlet edges $\\mathcal{E}(\\Gamma\_D)$ and the Neumann edges $\\mathcal{E}(\\Gamma\_N)$; let $\\mathcal{E}\_D = \\mathcal{E}(\\Omega)\\cup\\mathcal{E}(\\Gamma\_D)$. Integration by parts on every triangle gives edge integrals of $\\nabla u\\cdot n\_E$ on all edges; on the Neumann edges $\\nabla u\\cdot n = g$ is known and moves to the right-hand side, on the other edges the SIPG method of `exPoissonDG` replaces it by the average and adds the symmetric and the penalty terms. The method seeks $u\_h\\in V\_h = P\_k(\\mathcal{T})$ (`DGPk(degree,mesh)`) with
%[text]{"align":"center"} $ a\_h(u\_h,v\_h) = (f,v\_h) + \\sum\_{E\\in\\mathcal{E}(\\Gamma\_N)} (g,v\_h)\_E \\quad \\textrm{for all } v\_h\\in V\_h, $
%[text]{"align":"center"} $ a\_h(u\_h,v\_h) = (\\nabla\_h u\_h, \\nabla\_h v\_h)\_{\\mathcal{T}} - \\sum\_{E\\in\\mathcal{E}\_D} \\left( (\\{\\nabla\_h u\_h\\}\\cdot n\_E, \[v\_h\])\_E + (\[u\_h\], \\{\\nabla\_h v\_h\\}\\cdot n\_E)\_E - \\frac{\\alpha}{h\_E}(\[u\_h\],\[v\_h\])\_E \\right). $
%[text] The Dirichlet edges impose $u = u\_D = 0$ weakly; for $u\_D \\ne 0$ the right-hand side would contain the terms $\\frac{\\alpha}{h\_E}(u\_D, v\_h)\_E - (u\_D, \\nabla v\_h\\cdot n)\_E$ on these edges (see `exStokesDG`). The form strings repeat the terms of `exPoissonDG` with the domains `dsi` and `dsd`:
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $(\\nabla\_h u, \\nabla\_h v)\_{\\mathcal{T}}$ | `grad(u)*grad(v) dx` |
%[text] | $-\\sum\_{E\\in\\mathcal{E}(\\Omega)} \\left((\[u\], \\{\\nabla\_h v\\}\\cdot n\_E)\_E + (\\{\\nabla\_h u\\}\\cdot n\_E, \[v\])\_E - \\frac{\\alpha}{h\_E}(\[u\],\[v\])\_E\\right)$ | `- [u]*n.{grad(v)} dsi - {grad(u)}.[v]*n dsi + h^-1*alpha(x)*[u][v] dsi` |
%[text] | the same sum over $\\mathcal{E}(\\Gamma\_D)$ | the same terms with `dsd` |
%[text] | $(f,v) + \\sum\_{E\\in\\mathcal{E}(\\Gamma\_N)} (g,v)\_E$ | `f(x) v dx + g(x) v dsn` |
%[text:table]
%[text] The error is measured in the DG norm $\\|v\\|\_{DG}^2 = \\|\\nabla\_h v\\|^2\_{L^2(\\Omega)} + \\sum\_{E\\in\\mathcal{E}\_D} h\_E^{-1}\\|\[v\]\\|^2\_{L^2(E)}$; the jumps of $u$ vanish on $\\mathcal{E}\_D$.
%[text] ### Error estimator
%[text] The residual-based error estimator $\\eta^2 = \\sum\_{T\\in\\mathcal{T}}\\eta\_T^2$ consists of the residuals of the conforming method of `exPoissonPkAdaptive`, i.e. the residual of the equation on the triangles, the jumps of the normal derivatives across the interior edges and the residual of the Neumann condition, and of the penalty-weighted jumps of $u\_h$ on $\\mathcal{E}\_D$:
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{T\\in\\mathcal{T}} h\_T^2\\|f+\\Delta u\_h\\|^2\_{L^2(T)} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} h\_E\\|\[\\nabla u\_h\\cdot n\_E\]\\|^2\_{L^2(E)} + \\sum\_{E\\in\\mathcal{E}(\\Gamma\_N)} h\_E\\|g-\\nabla u\_h\\cdot n\\|^2\_{L^2(E)} + \\sum\_{E\\in\\mathcal{E}\_D} \\frac{\\alpha}{h\_E}\\|\[u\_h\]\\|^2\_{L^2(E)}. $
%[text] The last sum is new: a discontinuous $u\_h$ is not in $H^1(\\Omega)$, and its jumps on the interior edges and its values on $\\Gamma\_D$ (where $u\_D = 0$) bound its distance to the continuous functions with the boundary values $u\_D$, the nonconforming part of the error. `EstimatorForm` squares the whole term $(\\dots)^2$, so the coefficient `alpha(x)` inside the jump receives $\\sqrt{\\alpha}$; an edge term counts one half for each of the two neighbouring triangles.
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2\\Vert f+\\Delta u\_h\\Vert^2\_{L^2(T)}$ | `h^2*(f(x)+div(grad(u)))^2 dx` |
%[text] | $h\_E\\Vert\[\\nabla u\_h\\cdot n\_E\]\\Vert^2\_{L^2(E)}$ | `h^1*([grad(u)*n])^2 dsi` |
%[text] | $h\_E\\Vert g-\\nabla u\_h\\cdot n\\Vert^2\_{L^2(E)}$ | `h^1*(g(x)-grad(u)*n)^2 dsn` |
%[text] | $\\frac{\\alpha}{h\_E}\\Vert\[u\_h\]\\Vert^2\_{L^2(E)}$ on $\\mathcal{E}(\\Omega)$ and $\\mathcal{E}(\\Gamma\_D)$ | `h^-1*([alpha(x)*u])^2 dsi` and `h^-1*([alpha(x)*u])^2 dsd` |
%[text:table]
%[text] The estimator is reliable, $\\|u-u\_h\\|\_{DG} \\le C\_{\\textrm{rel}}\\,\\eta$, and efficient up to oscillations of the data. The adaptive loop marks the triangles by bulk (Dörfler) marking, `MarkBulk` with $\\theta = 1/2$, and refines the marked edges by newest vertex bisection, `mesh.refine(mark.marked4ed,'nvb')`.
%[text] ### Expected convergence
%[text] Uniform refinement converges only like $\\|u-u\_h\\|\_{DG} = O(h^{2/3}) = O(N^{-1/3})$ for every degree $k$ because of the corner singularity. The adaptive algorithm recovers the rate $N^{-k/2}$ of a smooth solution, and the estimator decreases with the same rate. The default $k = 3$ makes the difference large: $N^{-3/2}$ instead of $N^{-1/3}$.
%%
%[text] ## Parameters
degree = 3;                     % polynomial degree k
penalty = 2*degree*(degree+1);  % penalty parameter alpha
minNrDoF = 1e4;                 % stop as soon as the number of degrees of freedom reaches minNrDoF
mark = MarkBulk;                % bulk (Doerfler) marking with theta = 1/2
% mark = MarkUniform;           % uniform refinement, for comparison
%%
%[text] ## Problem data and initial mesh
p = problemScalarEllipticLshapeExact();
mesh = Mesh(p.geom);
plot(mesh);
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the loop stops when the number of degrees of freedom reaches `minNrDoF`.
nrDoF4lvl = [];
DGError4lvl = [];
eta4lvl = [];
while true
    % finite element space
    fe = DGPk(degree,mesh);
    nrDoF4lvl(end+1) = fe.nrDoF;
    % assemble
    a = BilinearForm(['grad(u)*grad(v) dx', ...
        ' - [u]*n.{grad(v)} dsi - {grad(u)}.[v]*n dsi + h^-1*alpha(x)*[u][v] dsi', ...
        ' - [u]*n.{grad(v)} dsd - {grad(u)}.[v]*n dsd + h^-1*alpha(x)*[u][v] dsd']);
    A = a.assemble(fe,'alpha',penalty);
    l = LinearForm('f(x) v dx + g(x) v dsn');
    b = l.assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
    % solve
    x = A\b;
    % estimate
    eta = EstimatorForm(['h^2*(f(x)+div(grad(u)))^2 dx', ...
        ' + h^1*([grad(u)*n])^2 dsi + h^-1*([alpha(x)*u])^2 dsi', ...
        ' + h^1*(g(x)-grad(u)*n)^2 dsn + h^-1*([alpha(x)*u])^2 dsd']);
    eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g(x,n),'alpha',sqrt(penalty));
    eta4lvl(end+1) = norm(eta4e);
    % errors
    dgNorm = ErrorForm('( exact(x) - grad(u) )^2 dx + h^-1*([u])^2 dsi + h^-1*([u])^2 dsd');
    DGError4lvl(end+1) = dgNorm.assemble(x,fe,2*degree+2,'exact',@p.exactGradU);
    % stop
    if fe.nrDoF >= minNrDoF
        break
    end
    % mark
    mark = mark.mark(fe,eta4e);
    % refine
    mesh = mesh.refine(mark.marked4ed,'nvb');
end
%%
%[text] ## Results
%[text] ### Convergence table
%[text] The last column is the ratio of the error and the estimator.
table(nrDoF4lvl',DGError4lvl',eta4lvl',(DGError4lvl./eta4lvl)', ...
    VariableNames=["N" "DG error" "eta" "error/eta"])
%[text] The ratio stays bounded: up to the unknown constant $C\_{\\textrm{rel}}$ the estimator reflects the error on every mesh.
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
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,DGError4lvl,'||u-u_h||_{DG}');
hold on
plotConvergence(nrDoF4lvl,eta4lvl,'\eta');
hold off
%[text] After a preasymptotic phase on the coarse meshes, error and estimator decrease with the rate $k/2$ of a smooth solution, far above the rate $1/3$ of uniform refinement. The least-squares slopes in the legend include the preasymptotic phase and are therefore smaller than $k/2$.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkUniform;` to refine uniformly: the rate drops to $N^{-1/3}$.
%[text] - Set `degree = 1`: the adaptive rate is $N^{-1/2}$; compare the number of degrees of freedom for a given accuracy with that of `degree = 3`.
%[text] - Measure only the broken gradient, `ErrorForm('( exact(x) - grad(u) )^2 dx')`: on the adaptive meshes it decays with the same rate as the DG norm.
%[text] - Mark with `MarkMaximum` instead of `MarkBulk` and compare the meshes.
%[text] - `exPoissonPkAdaptive` solves the same problem with conforming elements and the estimator without the jump terms of $u\_h$. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
