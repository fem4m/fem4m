%[text] # Adaptive $P\_k$ finite elements for the Poisson problem in 3D
%[text] Solves the Poisson problem on an L-shaped prism, whose solution has a singular gradient along the reentrant edge, with conforming $P\_k$ finite elements on adaptively refined tetrahedral meshes. The algorithm is that of `exPoissonPkAdaptive`; new are the tetrahedral mesh `Mesh3d`, its refinement by newest vertex bisection and the convergence rates in three dimensions.
%[text] **Builds on:** `exPoissonPkAdaptive`. **Next:** `exLaplaceBeltramiPkAdaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] The L-shaped prism $\\Omega = \\big((-1,1)^2\\setminus(\[0,1\]\\times\[-1,0\])\\big)\\times(0,1)$ has a reentrant edge, the segment $x\_1 = x\_2 = 0$, $0\\le x\_3\\le 1$, at which the two faces $x\_1 = 0$ and $x\_2 = 0$ of the boundary meet at the interior angle $3\\pi/2$. Given $f$, $u\_D$ and $g$, find $u$ with
%[text]{"align":"center"} $ -\\Delta u = f \\textrm{ in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\Gamma\_D, \\qquad \\nabla u\\cdot n = g \\textrm{ on } \\Gamma\_N. $
%[text] The problem object `problemScalarEllipticLbrickExact` provides $f = 0$ and the Dirichlet data on the whole boundary, $\\Gamma\_D = \\partial\\Omega$ and $\\Gamma\_N = \\emptyset$, of the exact solution, in the polar coordinates $(r,\\varphi)$ of $(x\_1,x\_2)$ with $\\varphi\\in\[0,2\\pi)$,
%[text]{"align":"center"} $ u(x) = r^{2/3}\\sin(2\\varphi/3). $
%[text] It does not depend on $x\_3$ and vanishes on the two faces at the reentrant edge. Its gradient behaves like $r^{-1/3}$, the distance to the edge to the power $-1/3$, and is unbounded along the whole edge: $u$ belongs to $H^{1+s}(\\Omega)$ only for $s \< 2/3$. The Neumann terms below vanish for this problem; they keep the code ready for problems with Neumann faces.
%[text] ### Weak formulation
%[text] Find $u\\in H^1(\\Omega)$ with $u = u\_D$ on $\\Gamma\_D$ and
%[text]{"align":"center"} $ \\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx = \\int\_\\Omega f v\\,dx + \\int\_{\\Gamma\_N} g v\\,ds \\quad \\textrm{for all } v\\in H^1\_D(\\Omega) = \\{v\\in H^1(\\Omega) : v = 0 \\textrm{ on } \\Gamma\_D\\}. $
%[text] ### Discretization
%[text] The geometry `p.geom` holds the coarse tetrahedral mesh of the prism from the files `Lshape3d_*.dat`: 36 tetrahedra, twelve around the center of each of the three unit cubes, the Dirichlet faces and the bisection labels; `Mesh3d(p.geom)` builds the mesh from it. On a regular triangulation $\\mathcal{T}$ into tetrahedra the discrete space $V\_h = P\_k(\\mathcal{T})\\cap C(\\overline\\Omega)$ (`Pk3d(degree,mesh)`) has one degree of freedom per Lagrange node: at the vertices, $k-1$ on every edge, $(k-1)(k-2)/2$ on every face and $(k-1)(k-2)(k-3)/6$ inside every tetrahedron. As in 2D, the degrees of freedom on $\\Gamma\_D$ (`fe.fixedNodes`) take the values of $u\_D$, the others (`fe.freeNodes`) are the unknowns. The form strings are those of the 2D problem; the side domains `dsi` and `dsn` now denote the interior faces and the Neumann faces.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx$ | `grad(u) * grad(v) dx` |
%[text] | $\\int\_\\Omega f v\\,dx + \\int\_{\\Gamma\_N} g v\\,ds$ | `f(x) v dx + g(x) v dsn` |
%[text:table]
%[text] ### Error estimator
%[text] With the diameters $h\_T$ of the tetrahedra and $h\_F$ of the faces, the interior faces $\\mathcal{F}(\\Omega)$, the faces $\\mathcal{F}(\\Gamma\_N)$ on $\\Gamma\_N$ and a unit normal $n\_F$ of every face, the residual-based error estimator reads
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{T\\in\\mathcal{T}} h\_T^2 \\|f + \\Delta u\_h\\|^2\_{L^2(T)} + \\sum\_{F\\in\\mathcal{F}(\\Omega)} h\_F \\|\[\\nabla u\_h\\cdot n\_F\]\\|^2\_{L^2(F)} + \\sum\_{F\\in\\mathcal{F}(\\Gamma\_N)} h\_F \\|g - \\nabla u\_h\\cdot n\\|^2\_{L^2(F)}. $
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2\\Vert f+\\Delta u\_h\\Vert^2\_{L^2(T)}$ | `h^2*(f(x)+div(grad(u)))^2 dx` |
%[text] | $h\_F\\Vert\[\\nabla u\_h\\cdot n\_F\]\\Vert^2\_{L^2(F)}$ | `h^1*([grad(u)*n])^2 dsi` |
%[text] | $h\_F\\Vert g-\\nabla u\_h\\cdot n\\Vert^2\_{L^2(F)}$ | `h^1*(g(x)-grad(u)*n)^2 dsn` |
%[text:table]
%[text] It is reliable and efficient as in 2D: $\\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)} \\le C\_{\\textrm{rel}}\\,\\eta$ and, up to oscillations of the data, $\\eta \\le C\_{\\textrm{eff}}\\,\\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)}$.
%[text] ### Refinement by newest vertex bisection
%[text] `Mesh3d` refines by recursive bisection, the three-dimensional newest vertex bisection. Every tetrahedron stores its vertices in an order whose first and last vertex span its refinement edge (`mesh.n4e(:,[1 4])`), and a type. A bisection cuts the tetrahedron at the midpoint of its refinement edge into two children; their vertex order and type, and with them their refinement edges, follow from those of the parent. `mesh.refine(marked4e)` bisects every marked tetrahedron once, unlike the 2D refinement of all edges of a marked triangle, and then further tetrahedra until the mesh has no hanging nodes. The bisections produce only finitely many similarity classes of tetrahedra, so the meshes stay shape regular. Three rounds of bisections of all tetrahedra give eight children each, the uniform refinement `mesh.refine('uniform')`.
%[text] ### Adaptive algorithm
%[text] Starting from the coarse mesh, the adaptive finite element method repeats
%[text] 1. SOLVE: compute $u\_h$ on the current mesh,
%[text] 2. ESTIMATE: compute $\\eta\_T$ for all $T\\in\\mathcal{T}$,
%[text] 3. MARK: select a set $\\mathcal{M}\\subseteq\\mathcal{T}$ of minimal size with $\\theta\\,\\eta^2 \\le \\sum\_{T\\in\\mathcal{M}}\\eta\_T^2$ (`MarkBulk` with $\\theta = 1/2$),
%[text] 4. REFINE: bisect every marked tetrahedron once and complete the mesh \
%[text] until the number of degrees of freedom exceeds `minNrDoF`.
%[text] ### Expected convergence
%[text] In 3D the number of degrees of freedom grows like $N\\approx h^{-3}$, so the rate $h^k$ of a smooth solution is $N^{-k/3}$. On uniformly refined meshes the edge singularity limits the energy error to $O(h^{2/3}) = O(N^{-2/9})$ for every degree $k$. Adaptive refinement recovers $N^{-1/3}$ for $k = 1$. Shape regular tetrahedra cannot do better than $N^{-2/3}$ for this singularity: a layer of tetrahedra of diameter $\\delta$ along the edge of length one contains about $\\delta^{-1}$ of them and carries an error of order $\\delta^{2/3}$. So the default $k = 2$ reaches its optimal rate $N^{-2/3}$ only up to logarithmic factors, and $k\\ge 3$ stays at $N^{-2/3}$; anisotropic tetrahedra, stretched along the edge, are needed for more.
%%
%[text] ## Parameters
degree = 2;       % polynomial degree k
minNrDoF = 2e4;   % stop as soon as the number of degrees of freedom exceeds minNrDoF
mark = MarkBulk;  % bulk (Doerfler) marking with theta = 1/2
% mark = MarkUniform;  % bisect all tetrahedra once, for comparison
%%
%[text] ## Problem data and initial mesh
%[text] One bisection of all 36 tetrahedra of the coarse mesh gives the first mesh. The plot shows the boundary faces, the Dirichlet faces in red.
p = problemScalarEllipticLbrickExact();
mesh = Mesh3d(p.geom);
mesh = mesh.refine(true(mesh.nrElems,1));
plot(mesh,'bd');
view([30,30]);
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the loop stops when the number of degrees of freedom exceeds `minNrDoF`.
nrDoF4lvl = [];
H1Error4lvl = [];
eta4lvl = [];
while true
    % finite element space
    fe = Pk3d(degree,mesh);
    nrDoF4lvl(end+1) = fe.nrDoF;
    % assemble
    a = BilinearForm('grad(u) * grad(v) dx');
    A = a.assemble(fe);
    l = LinearForm('f(x) v dx + g(x) v dsn');
    b = l.assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
    % solve
    x = zeros(fe.nrDoF,1);
    x(fe.fixedNodes) = p.u4Db(fe.c4dof(fe.fixedNodes,:));
    b = b - A*x;
    free = fe.freeNodes;
    x(free) = A(free,free)\b(free);
    % estimate
    eta = EstimatorForm(['h^2*(f(x)+div(grad(u)))^2 dx', ...
        '+ h^1*([grad(u)*n])^2 dsi', ...
        '+ h^1*(g(x)-grad(u)*n)^2 dsn']);
    eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g(x,n));
    eta4lvl(end+1) = norm(eta4e);
    % errors
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(end+1) = h1.assemble(x,fe,'exact',@p.exactGradU);
    % stop
    if fe.nrDoF >= minNrDoF
        break
    end
    % mark
    mark = mark.mark(fe,eta4e);
    % refine
    mesh = mesh.refine(mark.marked4e);
end
%%
%[text] ## Results
%[text] ### Convergence table
%[text] The last column is the ratio of the error and the estimator.
table(nrDoF4lvl',H1Error4lvl',eta4lvl',(H1Error4lvl./eta4lvl)', ...
    VariableNames=["N" "H1 error" "eta" "error/eta"])
%[text] A marked tetrahedron is bisected only once, so the number of degrees of freedom grows slowly from one mesh to the next. The ratio of error and estimator stays almost constant.
%%
%[text] ### Adaptive mesh
%[text] The boundary faces of the final mesh.
clf
plot(mesh);
view([30,30]);
%[text] The refinement concentrates along the reentrant edge, where the gradient of $u$ is singular; since $u$ does not depend on $x\_3$, it spreads along the whole edge.
%%
%[text] ### Discrete solution
%[text] A 3D finite element function is drawn on faces of tetrahedra, colored by its values at their vertices. `plot(fe,x,X0,X1,Y0,Y1,Z0,Z1)` draws it on the boundary of the union of all tetrahedra whose vertices lie in the box $\[X\_0,X\_1\]\\times\[Y\_0,Y\_1\]\\times\[Z\_0,Z\_1\]$; empty or omitted bounds are the bounds of the mesh. `plot(fe,x,[],0)` keeps the tetrahedra with $x\_1\\le 0$, i.e. it cuts the prism along the plane $x\_1 = 0$ and shows the cut face $x\_1 = 0$, $0 \< x\_2 \< 1$, which ends at the reentrant edge.
clf
plot(fe,x,[],0);
view([30,30]);
colorbar;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(u-u_h)||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,eta4lvl,'\eta');
hold off
%[text] Error and estimator decrease with the same rate, somewhat below $2/3$ for $k = 2$, since shape regular meshes reach this rate only up to logarithmic factors, and far above the rate $2/9$ of uniform refinement.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkUniform;` to bisect all tetrahedra: the rate drops to about $N^{-2/9}$.
%[text] - Set `degree = 1`: the adaptive rate is $N^{-1/3}$, the optimal rate of $P\_1$ elements. With `degree = 3` the rate stays close to $N^{-2/3}$, the limit of shape regular meshes for this edge singularity.
%[text] - Use `problemScalarEllipticCubeExact`, the smooth solution $u = x\_1(1-x\_1)\\,x\_2(1-x\_2)\\,x\_3(1-x\_3)$ on the unit cube: uniform and adaptive refinement both converge with the rate $N^{-k/3}$.
%[text] - Use Neumann faces, e.g. `p.geom = Geometry(struct('name','Cube','midpoint',[0.5,0.5,0.5],'length',1,'bdry','DNDNDN'));` for the cube: the terms `g(x) v dsn` and `(g(x)-grad(u)*n)^2 dsn` become active.
%[text] - `exPoissonPkAdaptive` is the same algorithm in 2D; `exLaplaceBeltramiPkAdaptive` solves a Poisson problem on a curved surface in $\\mathbb{R}^3$. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
