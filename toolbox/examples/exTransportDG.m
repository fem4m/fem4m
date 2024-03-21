%[text] # Upwind DG method for the transport equation
%[text] Discretizes the stationary transport equation, a first-order hyperbolic problem, with discontinuous piecewise polynomials and the upwind numerical flux: every triangle receives information only through its inflow edges. The example compares the error with the difference of the solutions on consecutive meshes, which measures the convergence also for problems without exact solution.
%[text] **Builds on:** `exPoissonDG`. **Next:** `exEllipticCRAdaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] Given a velocity field $\\beta$ with $\\mathrm{div}\\,\\beta = 0$, a source $f$ and inflow data $u\_D$, find $u$ with
%[text]{"align":"center"} $ \\beta\\cdot\\nabla u = f \\quad \\textrm{in } \\Omega, \\qquad u = u\_D \\quad \\textrm{on } \\Gamma\_- = \\{x\\in\\partial\\Omega : \\beta(x)\\cdot n(x) \< 0\\}. $
%[text] The solution is transported along the characteristics, the curves tangent to $\\beta$. Information enters only through the inflow boundary $\\Gamma\_-$; on the rest of the boundary no condition is prescribed. The problem object `problemTransportSquareExact` on the unit square $\\Omega = (0,1)^2$ has the constant velocity $\\beta = (1,0)$: the inflow boundary is the left edge $x\_1 = 0$, the outflow boundary the right edge, and the top and bottom edges are parallel to $\\beta$. It provides $f = \\beta\\cdot\\nabla u = (1-2x\_1)\\,x\_2(1-x\_2)$ of the exact solution $u = x\_1(1-x\_1)\\,x\_2(1-x\_2)$, so $u\_D = 0$.
%[text] ### Discretization
%[text] The notation for jumps $\[\\cdot\]$, averages $\\{\\cdot\\}$, normals $n\_E$ and the edges $\\mathcal{E}(\\Omega)$, $\\mathcal{E}(\\partial\\Omega)$ is that of `exPoissonDG`. The discrete space $V\_h = P\_k(\\mathcal{T})$ (`DGPk(degree,mesh)`) has no continuity. Testing the equation with $v\_h\\in V\_h$ on every triangle and coupling the triangles by the jumps of $u\_h$ weighted with the normal velocity gives the method: find $u\_h\\in V\_h$ with $a\_h(u\_h,v\_h) = \\ell\_h(v\_h)$ for all $v\_h\\in V\_h$, where
%[text]{"align":"center"} $ a\_h(u\_h,v\_h) = (\\beta\\cdot\\nabla\_h u\_h, v\_h)\_{\\mathcal{T}} - \\sum\_{E\\in\\mathcal{E}(\\Omega)} ((\\beta\\cdot n\_E)\[u\_h\], \\{v\_h\\})\_E + \\sum\_{E\\in\\mathcal{E}(\\partial\\Omega)} ((\\beta\\cdot n)^- u\_h, v\_h)\_E + \\sum\_{E\\in\\mathcal{E}(\\Omega)} \\left(\\frac{\\gamma}{2}|\\beta\\cdot n\_E|\\,\[u\_h\], \[v\_h\]\\right)\_E, $
%[text]{"align":"center"} $ \\ell\_h(v\_h) = (f, v\_h) + \\sum\_{E\\in\\mathcal{E}(\\partial\\Omega)} ((\\beta\\cdot n)^- u\_D, v\_h)\_E, \\qquad (\\beta\\cdot n)^- = \\max(-\\beta\\cdot n, 0). $
%[text] The factor $(\\beta\\cdot n)^-$ vanishes outside $\\Gamma\_-$, so the boundary terms impose the inflow condition weakly. All edge terms vanish for the exact solution, which is continuous and equals $u\_D$ on $\\Gamma\_-$: the method is consistent for every weight $\\gamma\\ge 0$ (`upwind` in the code).
%[text] - $\\gamma = 0$, the centered flux: the edges see the average $\\{v\_h\\}$ of the test function.
%[text] - $\\gamma = 1$, the upwind flux: on every triangle $T$ the edge terms add up to $(|\\beta\\cdot n\_T|\\,(u\_h - u\_h^{\\textrm{up}}), v\_h)\_{\\partial T\_-}$ on the inflow part $\\partial T\_- = \\{x\\in\\partial T: \\beta\\cdot n\_T \< 0\\}$ of its boundary, where $u\_h^{\\textrm{up}}$ is the trace of $u\_h$ on the upwind neighbour, or $u\_D$ on $\\Gamma\_-$. A triangle receives information only from upstream, as the exact solution does. \
%[text] Integration by parts with $\\mathrm{div}\\,\\beta = 0$ shows the stability of the method,
%[text]{"align":"center"} $ a\_h(v\_h,v\_h) = \\frac{1}{2}\\sum\_{E\\in\\mathcal{E}(\\partial\\Omega)} \\left\\Vert |\\beta\\cdot n|^{1/2}v\_h\\right\\Vert^2\_{L^2(E)} + \\frac{\\gamma}{2}\\sum\_{E\\in\\mathcal{E}(\\Omega)} \\left\\Vert |\\beta\\cdot n\_E|^{1/2}\[v\_h\]\\right\\Vert^2\_{L^2(E)}: $
%[text] the upwind flux controls the jumps of $u\_h$ in addition, the centered flux leaves them free.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $(\\beta\\cdot\\nabla\_h u, v)\_{\\mathcal{T}}$ | `beta(x)*grad(u)v dx` |
%[text] | $-\\sum\_{E\\in\\mathcal{E}(\\Omega)} ((\\beta\\cdot n\_E)\[u\], \\{v\\})\_E$ | `- betan(x)*[u]{v} dsi` |
%[text] | $\\sum\_{E\\in\\mathcal{E}(\\partial\\Omega)} ((\\beta\\cdot n)^- u, v)\_E$ | `betanm(x)*[u]{v} dsd` |
%[text] | $\\sum\_{E\\in\\mathcal{E}(\\Omega)} (\\frac{\\gamma}{2}\\vert\\beta\\cdot n\_E\\vert\[u\], \[v\])\_E$ | `delta(x)*[u][v] dsi` |
%[text] | $(f,v) + \\sum\_{E\\in\\mathcal{E}(\\partial\\Omega)} ((\\beta\\cdot n)^- u\_D, v)\_E$ | `f(x) v dx + inflow(x) v dsd` |
%[text:table]
%[text] All boundary edges of the square are Dirichlet edges (`dsd`); there jump and average are the one-sided values, $\[u\]\\{v\\} = uv$. The coefficients `betan`, `betanm` and `delta` need the normal $n\_E$: a coefficient handle with the three arguments `(x,xRef,parts)` receives the numbers `parts` of the edges and reads their normals from `mesh.normal4ed`, while the function `inflow` of the linear form, a handle `(x,n)`, receives the normals directly.
%[text] ### Measuring the error
%[text] For the default problem the exact solution is known. Without it, the difference of the solutions $u\_\\ell$ and $u\_{\\ell+1}$ on two consecutive meshes measures the error: the meshes are nested, so $u\_\\ell$ is also a function of the finer space, and `fe.prolongate(xOld,feOld)` gives its coefficients there. If the error decreases by a factor $q \< 1$ from one mesh to the next, the triangle inequality gives
%[text]{"align":"center"} $ (1-q)\\,\\|u-u\_\\ell\\|\_{L^2(\\Omega)} \\le \\|u\_{\\ell+1}-u\_\\ell\\|\_{L^2(\\Omega)} \\le (1+q)\\,\\|u-u\_\\ell\\|\_{L^2(\\Omega)}, $
%[text] so the difference decays with the same rate as the error of $u\_\\ell$. The example computes both; the difference belongs to the coarser of the two meshes, so it is missing on the finest mesh.
%[text] ### Expected convergence
%[text] For smooth solutions the centered flux converges with the order $k$, $\\|u-u\_h\\|\_{L^2(\\Omega)} = O(h^k)$, and the upwind flux with the order $k+1/2$. The order $k+1/2$ is sharp on specially constructed meshes; on many meshes, among them the uniform meshes of this example, the upwind method converges even with the optimal order $k+1$. In terms of $N\\approx h^{-2}$ these are the rates $k/2$, $(2k+1)/4$ and $(k+1)/2$.
%%
%[text] ## Parameters
degree = 1;    % polynomial degree k
upwind = 1;    % weight gamma of the jump term: the upwind flux
% upwind = 0;  % centered flux: no control of the jumps, lower order for odd k
nrLevels = 6;  % number of meshes
%%
%[text] ## Problem data and initial mesh
p = problemTransportSquareExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform','nvb');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
%[text] Besides the error, the loop computes the difference to the solution of the previous mesh; the $L^2$ error stays `NaN` for a problem without exact solution.
nrDoF4lvl = zeros(1,nrLevels);
L2Error4lvl = NaN(1,nrLevels);
difference4lvl = NaN(1,nrLevels);
for lvl = 1:nrLevels
    % finite element space
    fe = DGPk(degree,mesh);
    nrDoF4lvl(lvl) = fe.nrDoF;
    % assemble
    normal4ed = mesh.normal4ed;
    betan = @(x,parts) sum(p.beta(x).*normal4ed(parts,:),2);
    a = BilinearForm(['beta(x)*grad(u)v dx - betan(x)*[u]{v} dsi', ...
        ' + betanm(x)*[u]{v} dsd + delta(x)*[u][v] dsi']);
    A = a.assemble(fe,'beta',@p.beta, ...
        'betan',@(x,xRef,parts) betan(x,parts), ...
        'betanm',@(x,xRef,parts) max(-betan(x,parts),0), ...
        'delta',@(x,xRef,parts) upwind/2*abs(betan(x,parts)));
    l = LinearForm('f(x) v dx + inflow(x) v dsd');
    b = l.assemble(fe,'f',@p.f,'inflow',@(x,n) max(-sum(p.beta(x).*n,2),0).*p.u4Db(x));
    % solve
    x = A\b;
    % errors
    if isa(p,'ProblemScalarEllipticTemplateExact')
        l2 = ErrorForm('( exact(x) - u )^2 dx');
        L2Error4lvl(lvl) = l2.assemble(x,fe,12,'exact',@p.exactU);
    end
    if lvl > 1
        l2Difference = ErrorForm('( u - v )^2 dx');
        difference4lvl(lvl-1) = l2Difference.assemble(x,fe,fe.prolongate(xOld,feOld),fe,2*degree+2);
    end
    % refine
    feOld = fe;
    xOld = x;
    mesh = mesh.refine('uniform','nvb');
end
%%
%[text] ## Results
%[text] ### Convergence table
%[text] The column difference holds $\\|u\_{\\ell+1}-u\_\\ell\\|\_{L^2(\\Omega)}$ in the row of the mesh $\\ell$.
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
differenceRate = [NaN, -diff(log(difference4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',difference4lvl',differenceRate', ...
    VariableNames=["N" "L2 error" "L2 rate" "difference" "difference rate"])
%[text] The upwind error converges with the optimal rate $(k+1)/2$ on these meshes, and the difference of consecutive solutions follows the error with the same rate and a ratio close to one.
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
plotConvergence(nrDoF4lvl(1:end-1),difference4lvl(1:end-1),'||u_{h/2}-u_h||_{L^2}');
hold off
%%
%[text] ## Possibilities for extensions
%[text] - Set `upwind = 0`: with the centered flux the error of `degree = 1` decays only like $N^{-1/2}$, the order $k$, and is much larger than with the upwind flux.
%[text] - Set `degree = 2`: the upwind error decays like $N^{-3/2}$; on these meshes the centered flux reaches this rate as well for $k = 2$.
%[text] - Use `p = problemTransportSquare2();`, which has no exact solution: the inflow data jump at two boundary points, and the solution jumps across the characteristics through them. Only the differences measure the convergence, and they decay slowly; the $L^2$ error column and its curve stay empty. Near the discontinuities the upwind solution overshoots the range of the data, the centered one oscillates much more. Its velocity has $\\mathrm{div}\\,\\beta = -11$, and the stability identity gains the positive term $-\\frac{1}{2}(\\mathrm{div}\\,\\beta\\, v\_h, v\_h)$.
%[text] - `exEllipticCRAdaptive` treats the elliptic operator $-\\mathrm{div}(\\alpha\\nabla u) + \\beta\\cdot\\nabla u + \\gamma u$ with diffusion, convection and reaction by Crouzeix-Raviart elements. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
