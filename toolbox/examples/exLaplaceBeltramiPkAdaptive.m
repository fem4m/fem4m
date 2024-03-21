%[text] # Adaptive surface finite elements for the Laplace-Beltrami problem
%[text] Solves the Laplace-Beltrami equation on a three-quarter sphere, whose solution has singular derivatives at the two poles, with $P\_k$ finite elements on an adaptively refined triangulation of the surface. New are the surface mesh `MeshSurface`, whose vertices are projected onto the surface, the lift of the data onto the discrete surface by `SurfaceTransformation`, and a geometric error estimator that, together with the residual estimator, steers the refinement.
%[text] **Builds on:** `exPoissonPkAdaptive`. **Next:** `exEigenvaluePkAdaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] Let $\\Gamma\\subset\\mathbb{R}^3$ be a smooth surface, the zero level set of a function $\\phi$ with the unit normal $\\nu = \\nabla\\phi/|\\nabla\\phi|$. The tangential gradient $\\nabla\_\\Gamma u = (I - \\nu\\nu^T)\\nabla\\tilde u$ of a function $u$ on $\\Gamma$, with any extension $\\tilde u$ to a neighbourhood, and the surface divergence define the Laplace-Beltrami operator $\\Delta\_\\Gamma u = \\mathrm{div}\_\\Gamma\\nabla\_\\Gamma u$. Given $f$ and $u\_D$, find $u$ with
%[text]{"align":"center"} $ -\\Delta\_\\Gamma u = f \\textrm{ on } \\Gamma, \\qquad u = u\_D \\textrm{ on } \\partial\\Gamma. $
%[text] The problem object `problemSurfaceSphereLExact` takes the unit sphere $|x| = 1$ without the wedge $\\{x\_1 \> 0, x\_2 \< 0\\}$. The boundary of this three-quarter sphere consists of two half great circles, which meet at the north and the south pole at the interior angle $3\\pi/2$, like the reentrant corner of the L-shaped domain. With the polar angle $\\theta$ and the azimuth $\\varphi\\in\[0,3\\pi/2\]$ the exact solution and the source are
%[text]{"align":"center"} $ u = \\sin(\\theta)^{2/3}\\sin(2\\varphi/3), \\qquad f = \\frac{10}{9}\\,u, $
%[text] since $-\\Delta\_\\Gamma u = \\alpha(\\alpha+1)u$ for $\\alpha = 2/3$. The solution vanishes on $\\partial\\Gamma$, $u\_D = 0$. Near a pole $\\sin\\theta$ agrees with the geodesic distance $r$ to the pole up to $O(r^3)$, so $u$ behaves like $r^{2/3}\\sin(2\\varphi/3)$ in the plane and its gradient is unbounded at both poles.
%[text] ### Weak formulation
%[text] Find $u\\in H^1(\\Gamma)$ with $u = u\_D$ on $\\partial\\Gamma$ and
%[text]{"align":"center"} $ \\int\_\\Gamma \\nabla\_\\Gamma u\\cdot\\nabla\_\\Gamma v\\,ds = \\int\_\\Gamma f v\\,ds \\quad \\textrm{for all } v\\in H^1\_0(\\Gamma). $
%[text] On a closed surface, $\\partial\\Gamma = \\emptyset$, the solution is unique only up to a constant and requires $\\int\_\\Gamma f\\,ds = 0$; the side condition $\\int\_\\Gamma u\\,ds = 0$ fixes the constant.
%[text] ### Discretization
%[text] The surface finite element method works on a discrete surface $\\Gamma\_h$ close to $\\Gamma$. `MeshSurface(p.geom,@p.surfaceFunction,@p.gradSurfaceFunction,surfaceDegree)` projects the vertices of a coarse triangulation onto $\\Gamma$ (closest point projection along the normal), and every refinement projects the new vertices onto $\\Gamma$ again. For `surfaceDegree` $r = 1$ the triangles of $\\Gamma\_h$ are flat; for $r \> 1$ they are curved isoparametric triangles, the images of $P\_r$ maps whose Lagrange nodes lie on $\\Gamma$. On $\\Gamma\_h$ the space $V\_h$ of continuous functions that are polynomials of degree $k$ on the reference triangle (`Pk(degree,mesh)`; for $r = 1$ polynomials on the flat triangles) gives the discrete problem: find $u\_h\\in V\_h$ with $u\_h = u\_D$ at the boundary nodes and
%[text]{"align":"center"} $ \\int\_{\\Gamma\_h} \\nabla\_{\\Gamma\_h} u\_h\\cdot\\nabla\_{\\Gamma\_h} v\_h\\,ds\_h = \\int\_{\\Gamma\_h} f\_h v\_h\\,ds\_h \\quad \\textrm{for all } v\_h\\in V\_h \\textrm{ with } v\_h = 0 \\textrm{ on } \\partial\\Gamma\_h. $
%[text] The forms evaluate `grad(u)` as the tangential gradient on $\\Gamma\_h$ and integrate `dx` over $\\Gamma\_h$. The closest point projection $a(x)$ onto $\\Gamma$ lifts the data: $f\_h = \\delta\_h\\, f\\circ a$ with the ratio $\\delta\_h$ of the area elements, $ds = \\delta\_h\\,ds\_h$, so that $\\int\_{\\Gamma\_h} f\_h v\_h\\,ds\_h$ is the integral of $f$ times the lifted test function over $\\Gamma$. `q = SurfaceTransformation(mesh,p)` provides $f\_h$ as `q.f`. A function handle with the three arguments `(pts,ptsRef,parts)`, the points, the reference point and the element indices, receives the reference point that `q.f` needs for the normal of $\\Gamma\_h$.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_{\\Gamma\_h} \\nabla\_{\\Gamma\_h} u\\cdot\\nabla\_{\\Gamma\_h} v\\,ds\_h$ | `grad(u) * grad(v) dx` |
%[text] | $\\int\_{\\Gamma\_h} f\_h v\\,ds\_h$ | `f(x) v dx` |
%[text] | $\\int\_{\\Gamma\_h} v\\,ds\_h$ | `v dx` |
%[text:table]
%[text] On a closed surface the mesh has no boundary edges (`mesh.n4Db` is empty) and the coefficient vector $x$ of $u\_h$ and a Lagrange multiplier $\\rho$ for the side condition $\\int\_{\\Gamma\_h} u\_h\\,ds\_h = 0$ solve the saddle point system
%[text]{"align":"center"} $ \\left(\\begin{array}{cc} A & c \\\\ c^T & 0 \\end{array}\\right) \\left(\\begin{array}{c} x \\\\ \\rho \\end{array}\\right) = \\left(\\begin{array}{c} b \\\\ 0 \\end{array}\\right), \\qquad c\_j = \\int\_{\\Gamma\_h} \\varphi\_j\\,ds\_h. $
%[text] The multiplier also absorbs the small violation of $\\int\_{\\Gamma\_h} f\_h\\,ds\_h = 0$ by the discrete data.
%[text] ### Error estimator
%[text] The error has two sources: the finite element approximation on $\\Gamma\_h$ and the distance of $\\Gamma\_h$ from $\\Gamma$. The residual estimator measures the first,
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{T\\in\\mathcal{T}} h\_T^2 \\|f\_h + \\Delta\_{\\Gamma\_h} u\_h\\|^2\_{L^2(T)} + \\sum\_{E\\in\\mathcal{E}(\\Gamma\_h)} h\_E \\|\[\\nabla\_{\\Gamma\_h} u\_h\\cdot n\_E\]\\|^2\_{L^2(E)}, $
%[text] where $\\mathcal{E}(\\Gamma\_h)$ are the interior edges of $\\Gamma\_h$. The two triangles at an edge do not lie in one plane, so the jump is the sum of the derivatives $\\nabla\_{\\Gamma\_h} u\_h\\cdot n\_{E,T}$ in the directions of the outward co-normals $n\_{E,T}$ of both triangles $T$, the unit vectors tangent to $T$ and normal to $E$. The geometric estimator measures the second,
%[text]{"align":"center"} $ \\mu^2 = \\sum\_{T\\in\\mathcal{T}} \\|B\\,\\nabla\_{\\Gamma\_h} u\_h\\|^2\_{L^2(T)}, $
%[text] with the coefficient `q.B` of `SurfaceTransformation`. $B$ combines the signed distance $d$ of $\\Gamma\_h$ from $\\Gamma$, the principal curvatures $\\kappa\_1,\\kappa\_2$ of the parallel surface through the point and the normals $\\nu$ of $\\Gamma$ and $\\nu\_h$ of $\\Gamma\_h$, $B = \\big(|d(\\kappa\_1-\\kappa\_2)| + |1-\\nu\\cdot\\nu\_h|\\,(1+4\\max\_i(1-d\\kappa\_i))\\big)/\\delta\_h$. It vanishes for $\\Gamma\_h = \\Gamma$ and is of the order $h^{r+1}$, because $d = O(h^{r+1})$ and $1-\\nu\\cdot\\nu\_h = |\\nu-\\nu\_h|^2/2 = O(h^{2r})$ with $|\\nu-\\nu\_h| = O(h^r)$ for the surface degree $r$.
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2\\Vert f\_h+\\Delta\_{\\Gamma\_h} u\_h\\Vert^2\_{L^2(T)}$ | `h^2*(f(x)+div(grad(u)))^2 dx` |
%[text] | $h\_E\\Vert\[\\nabla\_{\\Gamma\_h} u\_h\\cdot n\_E\]\\Vert^2\_{L^2(E)}$ | `h^1*([grad(u)*n])^2 dsi` |
%[text] | $\\Vert B\\nabla\_{\\Gamma\_h} u\_h\\Vert^2\_{L^2(T)}$ | `(B(x)*grad(u))^2 dx` |
%[text:table]
%[text] Up to a constant, $\\eta + \\mu$ bounds the error of the lifted discrete solution on $\\Gamma$.
%[text] ### Adaptive algorithm
%[text] Starting from the coarse mesh, the adaptive finite element method repeats
%[text] 1. SOLVE: compute $u\_h$ on the current discrete surface,
%[text] 2. ESTIMATE: compute $\\eta\_T$ and $\\mu\_T$ for all $T\\in\\mathcal{T}$,
%[text] 3. MARK: `mark.mark(fe,eta4e,mu4e)` applies the bulk criterion to all indicators $\\eta\_T$ and $\\mu\_T$ together: it selects the largest of them until their squares sum to at least $\\theta(\\eta^2+\\mu^2)$ (`MarkBulk` with $\\theta = 1/2$), and marks the triangles of the selected indicators,
%[text] 4. REFINE: refine the marked triangles by red-green-blue refinement and project the new vertices onto $\\Gamma$ \
%[text] until the number of degrees of freedom exceeds `minNrDoF`.
%[text] ### Expected convergence
%[text] The surface is two-dimensional, so $N\\approx h^{-2}$ on uniform meshes as in the plane. The singularities at the poles limit uniform refinement to $\\|\\nabla\_\\Gamma(u-u\_h)\\|\_{L^2} = O(N^{-1/3})$ and $\\|u-u\_h\\|\_{L^2} = O(N^{-2/3})$. Adaptive refinement restores the rates $N^{-k/2}$ and $N^{-(k+1)/2}$ of a smooth solution, provided the surface degree satisfies $r\\ge k$: the geometric error is of the order $h^{r+1}$, and the errors below are measured on $\\Gamma\_h$ with the exact solution and its tangential gradient at the closest points on $\\Gamma$, which adds the angle $|\\nu-\\nu\_h| = O(h^r)$ between the tangent planes to the gradient error. For $k = r = 1$ the geometric estimator, asymptotically $\\mu = O(N^{-1})$, is of higher order than $\\eta = O(N^{-1/2})$.
%%
%[text] ## Parameters
degree = 1;         % polynomial degree k
surfaceDegree = 1;  % degree r of the discrete surface, 1 = flat triangles; choose r >= k
minNrDoF = 1e4;     % stop as soon as the number of degrees of freedom exceeds minNrDoF
mark = MarkBulk;    % bulk (Doerfler) marking with theta = 1/2
% mark = MarkUniform;  % uniform refinement, for comparison
%%
%[text] ## Problem data and initial mesh
%[text] One uniform refinement of the coarse triangulation of the three-quarter sphere, with all vertices on $\\Gamma$, gives the first mesh. The plots look into the removed wedge, at the two boundary half circles and the poles where they meet.
p = problemSurfaceSphereLExact();
mesh = MeshSurface(p.geom,@p.surfaceFunction,@p.gradSurfaceFunction,surfaceDegree);
mesh = mesh.refine('uniform');
plot(mesh);
axis equal;
view([45,20]);
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the loop stops when the number of degrees of freedom exceeds `minNrDoF`.
nrDoF4lvl = [];
L2Error4lvl = [];
H1Error4lvl = [];
eta4lvl = [];
mu4lvl = [];
while true
    % finite element space
    fe = Pk(degree,mesh);
    nrDoF4lvl(end+1) = fe.nrDoF;
    % assemble
    q = SurfaceTransformation(mesh,p);
    fh = @(pts,ptsRef,parts)(q.f(pts,ptsRef,parts));  % lifted source f_h
    a = BilinearForm('grad(u) * grad(v) dx');
    A = a.assemble(fe);
    l = LinearForm('f(x) v dx');
    b = l.assemble(fe,'f',fh);
    % solve
    if isempty(mesh.n4Db)
        % closed surface: mean value zero by a Lagrange multiplier
        meanValue = LinearForm('v dx');
        c = meanValue.assemble(fe);
        x = [A,c;c',0]\[b;0];
        x = x(1:end-1);
    else
        x = zeros(fe.nrDoF,1);
        x(fe.fixedNodes) = p.u4Db(fe.c4dof(fe.fixedNodes,:));
        b = b - A*x;
        free = fe.freeNodes;
        x(free) = A(free,free)\b(free);
    end
    % estimate
    eta = EstimatorForm(['h^2*(f(x)+div(grad(u)))^2 dx', ...
        '+ h^1*([grad(u)*n])^2 dsi']);
    eta4e = eta.assemble(x,fe,'f',fh);
    eta4lvl(end+1) = norm(eta4e);
    mu = EstimatorForm('(B(x)*grad(u))^2 dx');
    mu4e = mu.assemble(x,fe,'B',@(pts,ptsRef,parts)(q.B(pts,ptsRef,parts)));
    mu4lvl(end+1) = norm(mu4e);
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(end+1) = l2.assemble(x,fe,2*(degree+1),'exact',@p.exactU);
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(end+1) = h1.assemble(x,fe,2*(degree+1),'exact',@p.exactGradU);
    % stop
    if fe.nrDoF >= minNrDoF
        break
    end
    % mark
    mark = mark.mark(fe,eta4e,mu4e);
    % refine
    mesh = mesh.refine(mark.marked4ed);
end
%%
%[text] ## Results
%[text] ### Convergence table
%[text] The last column is the ratio of the error $\\|\\nabla\_\\Gamma(u-u\_h)\\|\_{L^2}$ and the residual estimator.
table(nrDoF4lvl',L2Error4lvl',H1Error4lvl',eta4lvl',mu4lvl',(H1Error4lvl./eta4lvl)', ...
    VariableNames=["N" "L2 error" "H1 error" "eta" "mu" "error/eta"])
%[text] The ratio stays bounded, and the geometric estimator is much smaller than the residual estimator.
%%
%[text] ### Adaptive mesh
clf
plot(mesh);
axis equal;
view([45,20]);
%[text] The refinement concentrates at the two poles, where the boundary curves meet and the gradient of $u$ is singular.
%%
%[text] ### Discrete solution
clf
plot(fe,x);
axis equal;
view([45,20]);
colorbar;
shading interp;
%[text] The discrete solution vanishes on the two boundary half circles and is largest on the far side of the sphere, at the equator point with the azimuth $\\varphi = 3\\pi/4$.
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla_\Gamma(u-u_h)||_{L^2}');
plotConvergence(nrDoF4lvl,eta4lvl,'\eta');
plotConvergence(nrDoF4lvl,mu4lvl,'\mu');
hold off
%[text] The gradient error and $\\eta$ decrease with the same rate, close to $k/2$, the $L^2$ error with a rate close to $(k+1)/2$, and $\\mu$ faster than $\\eta$: for $k = r = 1$ the geometric error does not limit the convergence.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkUniform;` to refine uniformly: the singularities at the poles reduce the rates to $N^{-1/3}$ for the gradient and $N^{-2/3}$ for the $L^2$ error.
%[text] - Set `degree = 2` and `surfaceDegree = 2`: the adaptive rates grow to $N^{-1}$ and $N^{-3/2}$. With `degree = 2` and flat triangles, `surfaceDegree = 1`, the angle between the tangent planes of $\\Gamma$ and $\\Gamma\_h$ dominates the gradient error measured on $\\Gamma\_h$: the rates stay lower and the ratio of the error and $\\eta$ grows from mesh to mesh, while $\\mu$, which bounds the geometric part of the error on $\\Gamma$, stays below $\\eta$.
%[text] - Use a closed surface, `problemSurfaceSphereExact`, `problemSurfaceTorusExact` or `problemSurfaceEllipsoidExact`: the mesh has no boundary edges, and the solve step imposes the mean value zero with the Lagrange multiplier.
%[text] - Call `mesh.plotCoNormal` after the construction of the first mesh: it draws the arrows of the co-normals of both triangles at every edge, which do not point in opposite directions; `mesh.plotSurfaceNormal` draws the normals $\\nu\_h$ of the triangles. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
