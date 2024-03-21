%[text] # hp-adaptive finite elements for the Poisson problem
%[text] Solves the Poisson problem on the slit domain, whose solution has a strong singularity at the tip of the slit, with conforming finite elements whose polynomial degree varies from triangle to triangle. On every marked triangle an hp-marking strategy decides whether to bisect the triangle (h-refinement) or to raise its degree (p-refinement); the error then decays exponentially in the cube root of the number of degrees of freedom.
%[text] **Builds on:** `exPoissonPkAdaptive`. **Next:** `exPoisson3dPkAdaptive`, `exEigenvaluePkAdaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] The slit domain $\\Omega = (-1,1)^2\\setminus(\[0,1\]\\times\\{0\\})$ is the square cut along the segment from the origin to the right edge. In polar coordinates $(r,\\varphi)$ with $\\varphi\\in\[0,2\\pi)$ the upper side of the slit is $\\varphi = 0$ and the lower side $\\varphi = 2\\pi$. The Neumann boundary $\\Gamma\_N$ is the lower side of the slit, the Dirichlet boundary $\\Gamma\_D$ consists of the upper side and the outer boundary. Given $f$, $u\_D$ and $g$, find $u$ with
%[text]{"align":"center"} $ -\\Delta u = f \\textrm{ in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\Gamma\_D, \\qquad \\nabla u\\cdot n = g \\textrm{ on } \\Gamma\_N. $
%[text] The problem object `problemScalarEllipticSlitExact` provides $f = 0$ and the boundary data of the exact solution
%[text]{"align":"center"} $ u(r,\\varphi) = r^{1/4}\\sin(\\varphi/4). $
%[text] It vanishes on the upper side of the slit, its normal derivative vanishes on the lower side, so $g = 0$, and $u\_D = u$ on the outer boundary. Its gradient behaves like $r^{-3/4}$ at the tip of the slit: $u$ belongs to $H^{1+s}(\\Omega)$ only for $s \< 1/4$, and away from the origin $u$ is analytic.
%[text] ### Weak formulation
%[text] Find $u\\in H^1(\\Omega)$ with $u = u\_D$ on $\\Gamma\_D$ and
%[text]{"align":"center"} $ \\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx = \\int\_\\Omega f v\\,dx + \\int\_{\\Gamma\_N} g v\\,ds \\quad \\textrm{for all } v\\in H^1\_D(\\Omega) = \\{v\\in H^1(\\Omega) : v = 0 \\textrm{ on } \\Gamma\_D\\}. $
%[text] ### Discretization
%[text] Every triangle $T\\in\\mathcal{T}$ carries its own polynomial degree $p\_T\\ge 1$, the entry of the vector `degree4e`. The discrete space
%[text]{"align":"center"} $ V\_{hp} = \\{v\\in C(\\overline\\Omega) : v|\_T\\in P\_{p\_T}(T) \\textrm{ for all } T\\in\\mathcal{T}\\} $
%[text] is `HPk(degree4e,mesh)`. Its hierarchical basis consists of the $P\_1$ nodal functions, integrated Legendre polynomials on the edges and element bubbles. An edge $E$ gets the smaller degree $p\_E$ of its two triangles (`fe.degree4ed`), which keeps the functions continuous; basis functions above the degree of their triangle or edge are not used. `fe.nrDoF` counts the used basis functions, while coefficient vectors and matrices have `fe.allnrDoF` entries or rows, and `fe.freeNodes` lists the used ones off $\\Gamma\_D$. The Dirichlet data are not polynomial; `fe.projectDbBoundary` computes their $L^2$ projection onto the traces of $V\_{hp}$ on $\\Gamma\_D$. The forms are those of `exPoissonPkAdaptive`:
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx$ | `grad(u) * grad(v) dx` |
%[text] | $\\int\_\\Omega f v\\,dx + \\int\_{\\Gamma\_N} g v\\,ds$ | `f(x) v dx + g(x) v dsn` |
%[text:table]
%[text] ### Error estimator
%[text] The residual-based estimator weights the terms of the estimator of `exPoissonPkAdaptive` with the local degrees: the residual on $T$ with $h\_T^2/p\_T^2$, the residuals on an edge with $h\_E/p\_E$,
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{T\\in\\mathcal{T}} \\frac{h\_T^2}{p\_T^2} \\|f + \\Delta u\_{hp}\\|^2\_{L^2(T)} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} \\frac{h\_E}{p\_E} \\|\[\\nabla u\_{hp}\\cdot n\_E\]\\|^2\_{L^2(E)} + \\sum\_{E\\in\\mathcal{E}(\\Gamma\_N)} \\frac{h\_E}{p\_E} \\|g - \\nabla u\_{hp}\\cdot n\\|^2\_{L^2(E)}. $
%[text] In the form strings `p` is $p\_T$ in the volume terms (`dx`) and $p\_E$ in the edge terms (`dsi`, `dsn`):
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2p\_T^{-2}\\Vert f+\\Delta u\_{hp}\\Vert^2\_{L^2(T)}$ | `h^2*p^-2*(f(x)+div(grad(u)))^2 dx` |
%[text] | $h\_Ep\_E^{-1}\\Vert\[\\nabla u\_{hp}\\cdot n\_E\]\\Vert^2\_{L^2(E)}$ | `h^1*p^-1*([grad(u)*n])^2 dsi` |
%[text] | $h\_Ep\_E^{-1}\\Vert g-\\nabla u\_{hp}\\cdot n\\Vert^2\_{L^2(E)}$ | `h^1*p^-1*(g(x)-grad(u)*n)^2 dsn` |
%[text:table]
%[text] The estimator is reliable, $\\|\\nabla(u-u\_{hp})\\|\_{L^2(\\Omega)} \\le C\\,\\eta$, with a constant $C$ independent of the mesh size and of the degrees; the efficiency constant grows with the degree.
%[text] ### hp-marking
%[text] `MarkHP` follows the predicted error reduction strategy. It keeps a prediction $\\eta\_T^{\\textrm{pred}}$ of the squared indicator of every triangle (`mark.predict4e`) and
%[text] 1. selects the triangles with $\\eta\_T^2 \\ge \\frac{3}{4}\\cdot\\frac{1}{|\\mathcal{T}|}\\sum\_{T'\\in\\mathcal{T}}\\eta\_{T'}^2$, at least three quarters of the mean,
%[text] 2. bisects a selected triangle with $\\eta\_T^2 \> \\eta\_T^{\\textrm{pred}}$ (h-refinement by newest vertex bisection of its refinement edge) and raises the degree of the other selected triangles by one (p-refinement),
%[text] 3. updates the prediction: $\\eta\_T^{\\textrm{pred}} = \\frac{\\gamma\_H}{2}(\\frac{1}{2})^{p\_T}\\eta\_T^2$ after h-refinement, inherited by both children, $\\eta\_T^{\\textrm{pred}} = \\gamma\_P\\,\\eta\_T^2$ after p-refinement, and $\\eta\_T^{\\textrm{pred}} \\leftarrow \\gamma\_N\\,\\eta\_T^{\\textrm{pred}}$ on the triangles not selected. \
%[text] The predictions are the indicators to expect if $u$ is smooth on $T$: a bisection halves the area, so $h\_T^{2p\_T}$ shrinks by $(1/2)^{p\_T}$, and the factor $1/2$ splits $\\eta\_T^2$ among the two children; a higher degree reduces the indicator of an analytic function by a fixed factor $\\gamma\_P \< 1$. An indicator above its prediction means that the last refinement of $T$ fell short of this reduction, so $u$ is not smooth near $T$ and $T$ is bisected; otherwise its degree is raised. The safety factors $\\gamma\_H = 4$, $\\gamma\_P = 0.4$ and $\\gamma\_N = 1$ are the defaults of `Mark`; the initial prediction $0$ (`mark.initialPredict`) bisects every selected triangle in the first step.
%[text] ### Expected convergence
%[text] On uniform meshes the singularity limits the convergence to $\\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)} = O(N^{-1/8})$ for every fixed degree; h-adaptive meshes with a fixed degree $p$ achieve $O(N^{-p/2})$, as in `exPoissonPkAdaptive`. For solutions that are analytic up to the corners and singular points, such as this $u$, hp-finite elements on meshes geometrically refined towards the singular points, with degrees that grow linearly away from them, converge exponentially,
%[text]{"align":"center"} $ \\|\\nabla(u-u\_{hp})\\|\_{L^2(\\Omega)} \\le C\\exp(-b N^{1/3}) $
%[text] with constants $b, C \> 0$. The hp-adaptive algorithm generates meshes and degree distributions of this kind; in a semilogarithmic plot over $N^{1/3}$ the error then lies on a straight line.
%%
%[text] ## Parameters
degree = 2;              % initial polynomial degree of all triangles
minNrDoF = 1e4;          % stop as soon as the number of degrees of freedom exceeds minNrDoF
mark = MarkHP;           % hp-marking by predicted error reduction
mark.gammaH = 4;         % prediction factor after h-refinement
mark.gammaP = 0.4;       % prediction factor after p-refinement
mark.gammaN = 1;         % prediction factor for triangles that are not refined
mark.initialPredict = 0; % first step: bisect every selected triangle
%%
%[text] ## Problem data and initial mesh
%[text] The coarse mesh of the slit domain has eight triangles; one uniform refinement by newest vertex bisection gives the first mesh, on which all triangles have the degree `degree`.
p = problemScalarEllipticSlitExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform','nvb');
degree4e = degree*ones(mesh.nrElems,1);
plot(mesh);
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the refinement first raises the degrees of the p-refined triangles and then bisects the h-refined ones, whose children inherit the degree of their parent `mesh.parent4e`.
nrDoF4lvl = [];
maxDegree4lvl = [];
H1Error4lvl = [];
eta4lvl = [];
while true
    % finite element space
    fe = HPk(degree4e,mesh);
    nrDoF4lvl(end+1) = fe.nrDoF;
    maxDegree4lvl(end+1) = fe.degree;
    % assemble
    a = BilinearForm('grad(u) * grad(v) dx');
    A = a.assemble(fe);
    l = LinearForm('f(x) v dx + g(x) v dsn');
    b = l.assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
    % solve
    x = fe.projectDbBoundary(@p.u4Db);
    b = b - A*x;
    free = fe.freeNodes;
    x(free) = A(free,free)\b(free);
    % estimate
    eta = EstimatorForm(['h^2*p^-2*(f(x)+div(grad(u)))^2 dx', ...
        '+ h^1*p^-1*([grad(u)*n])^2 dsi', ...
        '+ h^1*p^-1*(g(x)-grad(u)*n)^2 dsn']);
    eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g(x,n));
    eta4lvl(end+1) = norm(eta4e);
    % errors
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(end+1) = h1.assemble(x,fe,12,'exact',@p.exactGradU);
    % stop
    if fe.nrDoF >= minNrDoF
        break
    end
    % mark
    mark = mark.mark(fe,eta4e);
    % refine
    degree4e = degree4e + mark.markedP4e;
    if any(mark.marked4ed)
        mesh = mesh.refine('nvb',mark.marked4ed);
        degree4e = degree4e(mesh.parent4e);
    end
end
%%
%[text] ## Results
%[text] ### Convergence table
%[text] Every fourth level and the last one; the last column is the ratio of the error and the estimator.
rows = unique([1:4:numel(nrDoF4lvl), numel(nrDoF4lvl)]);
table(nrDoF4lvl(rows)',maxDegree4lvl(rows)',H1Error4lvl(rows)',eta4lvl(rows)', ...
    (H1Error4lvl(rows)./eta4lvl(rows))', ...
    VariableNames=["N" "max degree" "H1 error" "eta" "error/eta"])
%[text] The ratio stays bounded: the estimator reflects the error on every mesh, although the degrees grow.
%%
%[text] ### Polynomial degrees
%[text] The color of a triangle is its degree $p\_T$.
coordX = reshape(mesh.c4n(mesh.n4e,1),[],3);
coordY = reshape(mesh.c4n(mesh.n4e,2),[],3);
clf
patch(coordX',coordY',degree4e','EdgeColor','none');
axis equal
colorbar
%[text] The algorithm bisects towards the tip of the slit, where $u$ is singular, and keeps low degrees on the tiny triangles there, too small to be seen in this plot; the highest degrees sit on larger triangles at some distance from the tip. Near the outer corners, where the error is small, the triangles are seldom selected and their degrees stay lower. Around the tip this resembles the geometric meshes with degrees growing away from the singularity of the a priori theory.
%%
%[text] ### Discrete solution
clf
plot(fe,x);
view([40,35]);
colorbar;
shading interp;
%%
%[text] ### Convergence
%[text] Error and estimator over $N^{1/3}$ in a semilogarithmic scale.
clf
semilogy(nrDoF4lvl.^(1/3),H1Error4lvl,'-o','DisplayName','||\nabla(u-u_{hp})||_{L^2}');
hold on
semilogy(nrDoF4lvl.^(1/3),eta4lvl,'-o','DisplayName','\eta');
hold off
xlabel('N^{1/3}');
legend show;
%[text] Error and estimator follow nearly straight lines, which even steepen slightly: they decrease like $\\exp(-bN^{1/3})$, faster than any algebraic rate $N^{-s}$, which would bend upwards in this plot.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkBulk;`: the degrees stay at `degree` and the loop is h-adaptive with bulk marking; the error decays only algebraically, like $N^{-p/2}$.
%[text] - Set `mark = MarkHP3Edges; mark.initialPredict = 0;`: an h-refined triangle is split into four children by bisecting all its edges. Compare the errors at the same number of degrees of freedom.
%[text] - Vary the safety factors: a smaller `mark.gammaP` makes the next refinement of a p-refined triangle more likely an h-refinement, a larger `mark.gammaH` makes the children of a bisected triangle more likely p-refined. `mark.initialPredict = inf;`, the default of `MarkHP`, raises the degree of every selected triangle in the first step.
%[text] - Use `problemScalarEllipticLshapeExact` with the weaker corner singularity $r^{2/3}$: the degrees grow higher and the error decays much faster.
%[text] - `exPoissonPkAdaptive` adapts only the mesh; `exPoisson3dPkAdaptive` is the h-adaptive algorithm in 3D. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
