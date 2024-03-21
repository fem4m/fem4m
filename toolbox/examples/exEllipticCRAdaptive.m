%[text] # Adaptive Crouzeix-Raviart elements for a convection-diffusion-reaction problem
%[text] Solves a convection-diffusion-reaction problem with nonconforming Crouzeix-Raviart elements, piecewise linear functions that are continuous only at the midpoints of the edges, on adaptively refined meshes. The error estimator adds a term for the nonconformity, the jumps of the discrete solution across the edges.
%[text] **Builds on:** `exPoissonPkAdaptive`. **Next:** `exPoissonDG`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] Given a diffusion coefficient $\\alpha$ (a symmetric positive definite matrix at every point), a convection field $\\beta$, a reaction coefficient $\\gamma$, a source $f$, Dirichlet data $u\_D$ on $\\Gamma\_D$ and Neumann data $g$ on $\\Gamma\_N$, find $u$ with
%[text]{"align":"center"} $ -\\mathrm{div}(\\alpha\\nabla u) + \\beta\\cdot\\nabla u + \\gamma u = f \\textrm{ in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\Gamma\_D, \\qquad \\alpha\\nabla u\\cdot n = g \\textrm{ on } \\Gamma\_N. $
%[text] The problem object `problemScalarEllipticSquareFullEllipticExact` provides on the unit square $\\Omega = (0,1)^2$ the coefficients
%[text]{"align":"center"} $ \\alpha = I, \\qquad \\beta(x) = \\left(\\begin{array}{c} 5\\sin(x\_1+x\_2) \\\\ 6\\cos(x\_1+x\_2) \\end{array}\\right), \\qquad \\gamma = 7, $
%[text] the Dirichlet boundary $\\Gamma\_D = \\{0\\}\\times\[0,1\]$ (the left edge), the Neumann boundary $\\Gamma\_N$ of the other three edges, and $f$, $u\_D$, $g$ of the smooth exact solution
%[text]{"align":"center"} $ u(x) = \\sin(x\_1^3)\\cos(x\_2^\\pi) + x\_1^8 - x\_2^9 + x\_1^6x\_2^{10}. $
%[text] The high powers of $x\_1$ and $x\_2$ make its derivatives large near the right edge $x\_1 = 1$ and the upper edge $x\_2 = 1$.
%[text] ### Weak formulation
%[text] Find $u\\in H^1(\\Omega)$ with $u = u\_D$ on $\\Gamma\_D$ and
%[text]{"align":"center"} $ \\int\_\\Omega \\alpha\\nabla u\\cdot\\nabla v + (\\beta\\cdot\\nabla u)\\,v + \\gamma u v\\,dx = \\int\_\\Omega f v\\,dx + \\int\_{\\Gamma\_N} g v\\,ds \\quad \\textrm{for all } v\\in H^1\_D(\\Omega) = \\{v\\in H^1(\\Omega) : v = 0 \\textrm{ on } \\Gamma\_D\\}. $
%[text] ### Discretization
%[text] The Crouzeix-Raviart space consists of the piecewise linear functions that are continuous at the midpoints $\\mathrm{mid}(E)$ of the interior edges,
%[text]{"align":"center"} $ CR^1(\\mathcal{T}) = \\{v\\in P\_1(\\mathcal{T}) : v \\textrm{ is continuous at } \\mathrm{mid}(E) \\textrm{ for all } E\\in\\mathcal{E}(\\Omega)\\}, $
%[text] `CRk(1,mesh)`. Its degrees of freedom are the mean values on the edges, one per edge, which for a linear function are its values at the midpoints. The functions jump across the edges, so $CR^1(\\mathcal{T})\\not\\subset H^1(\\Omega)$: the method is nonconforming. The discrete problem replaces $\\nabla$ by the piecewise gradient $\\nabla\_{NC}$, the gradient on every triangle, which is what `grad(u)` evaluates. The Dirichlet condition fixes the degrees of freedom on $\\Gamma\_D$ (`fe.fixedNodes`) to the mean values of $u\_D$ on the edges, `fe.interpolateBoundary(@p.u4Db)`; the others (`fe.freeNodes`) are the unknowns. Find $u\_{CR}\\in CR^1(\\mathcal{T})$ with these boundary values and
%[text]{"align":"center"} $ \\int\_\\Omega \\alpha\\nabla\_{NC} u\_{CR}\\cdot\\nabla\_{NC} v + (\\beta\\cdot\\nabla\_{NC} u\_{CR})\\,v + \\gamma u\_{CR} v\\,dx = \\int\_\\Omega f v\\,dx + \\int\_{\\Gamma\_N} g v\\,ds $
%[text] for all $v\\in CR^1(\\mathcal{T})$ that vanish at the midpoints of the edges on $\\Gamma\_D$. The coefficients are passed to `assemble` by name, e.g. `'beta',@p.beta`.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega \\alpha\\nabla\_{NC} u\\cdot\\nabla\_{NC} v + (\\beta\\cdot\\nabla\_{NC} u)\\,v + \\gamma u v\\,dx$ | `alpha(x)*grad(u)*grad(v) dx + beta(x)*grad(u) v dx + gamma(x) u v dx` |
%[text] | $\\int\_\\Omega f v\\,dx + \\int\_{\\Gamma\_N} g v\\,ds$ | `f(x) v dx + g(x) v dsn` |
%[text:table]
%[text] ### Error estimator
%[text] The residual-based estimator adds to the residuals of the equation, of the normal fluxes $\\alpha\\nabla\_{NC}u\_{CR}\\cdot n\_E$ and of the Neumann condition a term that measures the nonconformity, the jumps $\[u\_{CR}\]$ across the interior edges, as the residual-based estimators for the Poisson problem do:
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{T\\in\\mathcal{T}} h\_T^2 \\|f + \\mathrm{div}(\\alpha\\nabla u\_{CR}) - \\beta\\cdot\\nabla u\_{CR} - \\gamma u\_{CR}\\|^2\_{L^2(T)} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} \\left(h\_E \\|\[\\alpha\\nabla\_{NC} u\_{CR}\\cdot n\_E\]\\|^2\_{L^2(E)} + h\_E^{-1}\\|\[u\_{CR}\]\\|^2\_{L^2(E)}\\right) + \\sum\_{E\\in\\mathcal{E}(\\Gamma\_N)} h\_E \\|g - \\alpha\\nabla\_{NC} u\_{CR}\\cdot n\\|^2\_{L^2(E)}. $
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2\\Vert f+\\mathrm{div}(\\alpha\\nabla u\_{CR})-\\beta\\cdot\\nabla u\_{CR}-\\gamma u\_{CR}\\Vert^2\_{L^2(T)}$ | `h^2*(f(x)+div(alphac(x)grad(u))-beta(x)*grad(u)-gamma(x)u)^2 dx` |
%[text] | $h\_E\\Vert\[\\alpha\\nabla\_{NC} u\_{CR}\\cdot n\_E\]\\Vert^2\_{L^2(E)}$ | `h^1*([alpha(x)*grad(u)*n])^2 dsi` |
%[text] | $h\_E\\Vert g-\\alpha\\nabla\_{NC} u\_{CR}\\cdot n\\Vert^2\_{L^2(E)}$ | `h^1*(g(x)-alpha(x)*grad(u)*n)^2 dsn` |
%[text] | $h\_E^{-1}\\Vert\[u\_{CR}\]\\Vert^2\_{L^2(E)}$ | `h^-1*([u])^2 dsi` |
%[text:table]
%[text] Inside `div(...)` the product rule $\\mathrm{div}(\\alpha\\nabla u\_{CR}) = \\sum\_{i,j}\\alpha\_{ij}\\partial\_i\\partial\_j u\_{CR} + \\sum\_{i,j}(\\partial\_i\\alpha\_{ij})\\,\\partial\_j u\_{CR}$ needs the derivatives of the coefficient. The pair `'alphac',{@p.alpha,@p.dalpha}` passes the coefficient together with its derivatives $\\partial\_k\\alpha\_{ij}$; for the piecewise linear $u\_{CR}$ only the second sum remains. The coefficient `alpha` of the edge terms needs no derivatives, it is passed alone. Since $\[u\_{CR}\]$ is linear on $E$ and vanishes at the midpoint, the nonconformity term equals $\\frac{h\_E}{12}\\|\[\\nabla\_{NC} u\_{CR}\\cdot t\_E\]\\|^2\_{L^2(E)}$, the jump of the tangential derivative. On $\\Gamma\_D$ the mean values of $u\_{CR}$ and $u\_D$ coincide; the remaining term $h\_E^{-1}\\|u\_D - u\_{CR}\\|^2\_{L^2(E)}$ is negligible for these smooth data and omitted.
%[text] ### Expected convergence
%[text] For smooth solutions the piecewise gradient converges like $\\|\\nabla\_{NC}(u-u\_{CR})\\|\_{L^2(\\Omega)} = O(h) = O(N^{-1/2})$ and the solution like $\\|u-u\_{CR}\\|\_{L^2(\\Omega)} = O(h^2) = O(N^{-1})$. The solution of this problem is smooth, so uniform refinement converges with these rates, too; adaptive refinement concentrates the triangles where the derivatives of $u$ are large and reduces the constants.
%%
%[text] ## Parameters
minNrDoF = 1e4;   % stop as soon as the number of degrees of freedom exceeds minNrDoF
mark = MarkBulk;  % bulk (Doerfler) marking with theta = 1/2
% mark = MarkUniform;  % uniform refinement, for comparison
%%
%[text] ## Problem data and initial mesh
p = problemScalarEllipticSquareFullEllipticExact();
% p = problemScalarEllipticLshapeExact();  % Poisson problem with a corner singularity
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the loop stops when the number of degrees of freedom exceeds `minNrDoF`.
nrDoF4lvl = [];
L2Error4lvl = [];
H1Error4lvl = [];
eta4lvl = [];
while true
    % finite element space
    fe = CRk(1,mesh);
    nrDoF4lvl(end+1) = fe.nrDoF;
    % assemble
    a = BilinearForm('alpha(x)*grad(u)*grad(v) dx + beta(x)*grad(u) v dx + gamma(x) u v dx');
    A = a.assemble(fe,'alpha',@p.alpha,'beta',@p.beta,'gamma',@p.gamma);
    l = LinearForm('f(x) v dx + g(x) v dsn');
    b = l.assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
    % solve
    x = fe.interpolateBoundary(@p.u4Db);
    b = b - A*x;
    free = fe.freeNodes;
    x(free) = A(free,free)\b(free);
    % estimate
    eta = EstimatorForm(['h^2*(f(x)+div(alphac(x)grad(u))-beta(x)*grad(u)-gamma(x)u)^2 dx', ...
        '+ h^1*([alpha(x)*grad(u)*n])^2 dsi', ...
        '+ h^1*(g(x)-alpha(x)*grad(u)*n)^2 dsn', ...
        '+ h^-1*([u])^2 dsi']);
    eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g(x,n),'alphac',{@p.alpha,@p.dalpha}, ...
        'alpha',@p.alpha,'beta',@p.beta,'gamma',@p.gamma);
    eta4lvl(end+1) = norm(eta4e);
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(end+1) = l2.assemble(x,fe,'exact',@p.exactU);
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(end+1) = h1.assemble(x,fe,'exact',@p.exactGradU);
    % stop
    if fe.nrDoF >= minNrDoF
        break
    end
    % mark
    mark = mark.mark(fe,eta4e);
    % refine
    mesh = mesh.refine(mark.marked4ed);
end
%%
%[text] ## Results
%[text] ### Convergence table
%[text] The rates are the slopes of the errors over $N$ in a log-log scale between two consecutive meshes; the last column is the ratio of the error of the piecewise gradient and the estimator.
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
H1Rate = [NaN, -diff(log(H1Error4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',H1Error4lvl',H1Rate',eta4lvl',(H1Error4lvl./eta4lvl)', ...
    VariableNames=["N" "L2 error" "L2 rate" "H1 error" "H1 rate" "eta" "error/eta"])
%[text] The ratio of error and estimator stays bounded: the estimator, including the nonconformity term, follows the error.
%%
%[text] ### Adaptive mesh
clf
plot(mesh);
%[text] The refinement concentrates along the right and the upper edge, where the derivatives of $u$ are large; the lower left part of the square stays coarse.
%%
%[text] ### Discrete solution
%[text] The discrete solution is discontinuous; on fine meshes the jumps are invisible.
clf
plot(fe,x);
view([40,35]);
colorbar;
shading interp;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_{CR}||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla_{NC}(u-u_{CR})||_{L^2}');
plotConvergence(nrDoF4lvl,eta4lvl,'\eta');
hold off
%[text] The least-squares slopes in the legend are close to 1 for the $L^2$ error and to 1/2 for the error of the piecewise gradient and the estimator, the rates $N^{-1}$ and $N^{-1/2}$ predicted for a smooth solution.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkUniform;`: the rates stay the same for this smooth solution, but the errors at the same number of degrees of freedom are larger.
%[text] - Use `p = problemScalarEllipticLshapeExact();`, the Poisson problem ($\\alpha = I$, $\\beta = 0$, $\\gamma = 0$, the defaults of `ProblemScalarEllipticTemplate`) with the corner singularity of `exPoissonPkAdaptive`: uniform refinement converges only like $N^{-1/3}$, adaptive refinement restores $N^{-1/2}$.
%[text] - Compare with conforming elements: use `fe = Pk(1,mesh);`, set the Dirichlet values by `x = zeros(fe.nrDoF,1); x(fe.fixedNodes) = p.u4Db(fe.c4dof(fe.fixedNodes,:));` and drop the jump term `h^-1*([u])^2 dsi`, which vanishes for continuous functions.
%[text] - Write a problem class of your own with other coefficients: a subclass of `ProblemScalarEllipticTemplateExact` overrides `alpha`, `dalpha`, `beta` and `gamma`; `exEigenvalueDiffusionPkAdaptive` uses a variable diffusion matrix.
%[text] - `exPoissonDG` lets the discrete functions jump across all edges and penalizes the jumps instead. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
