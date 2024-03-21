%[text] # Adaptive $P\_k$ finite elements for linear elasticity
%[text] Computes the displacement of a linear elastic body on an L-shaped domain, whose stress is singular at the reentrant corner, with vector-valued conforming $P\_k$ finite elements on adaptively refined meshes. Compared with the Poisson problem, the unknown is a vector field, the bilinear form involves the symmetric gradient and the elasticity tensor, and the error estimator measures the jumps of the tractions.
%[text] **Builds on:** `exPoissonPkAdaptive`. **Next:** `exElasticityAW`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] A linear elastic body occupies the domain $\\Omega\\subset\\mathbb{R}^2$. Its displacement $u:\\Omega\\to\\mathbb{R}^2$ defines the linearized strain $\\varepsilon(u) = (\\nabla u + \\nabla u^T)/2$, and Hooke's law gives the stress $\\sigma = \\mathbb{C}\\varepsilon(u)$ with the elasticity tensor
%[text]{"align":"center"} $ \\mathbb{C}\\varepsilon = 2\\mu\\,\\varepsilon + \\lambda\\,\\mathrm{tr}(\\varepsilon)\\,I, \\qquad \\lambda = \\frac{E\\nu}{(1+\\nu)(1-2\\nu)}, \\qquad \\mu = \\frac{E}{2(1+\\nu)}. $
%[text] The Lamé parameters $\\lambda$ and $\\mu$ follow from Young's modulus $E$ and the Poisson ratio $\\nu$; every elasticity problem of fem4m inherits them from `ElasticityParameters`, which also provides the matrix `p.C` of $\\mathbb{C}$ and `p.CInv` of $\\mathbb{C}^{-1}$ acting on the components $\[\\varepsilon\_{11}, \\varepsilon\_{22}, \\varepsilon\_{12}\]$ of a symmetric matrix. Given the body force $f$, the displacement $u\_D$ on $\\Gamma\_D$ and the traction $g$ on $\\Gamma\_N$, find $u$ with
%[text]{"align":"center"} $ -\\mathrm{div}\\,\\mathbb{C}\\varepsilon(u) = f \\textrm{ in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\Gamma\_D, \\qquad \\mathbb{C}\\varepsilon(u)\\,n = g \\textrm{ on } \\Gamma\_N. $
%[text] The problem object `problemElasticityLshapeExact` uses the L-shaped domain $(-1,1)^2\\setminus(\[0,1\]\\times\[-1,0\])$, rotated by $-3\\pi/4$ about its reentrant corner at the origin, so that the two reentrant edges lie on the rays $\\varphi = \\pm 3\\pi/4$. These edges are traction free ($\\Gamma\_N$ with $g = 0$), the other edges form $\\Gamma\_D$, and $f = 0$. The exact solution is the corner singularity
%[text]{"align":"center"} $ u(r,\\varphi) = \\frac{r^\\alpha}{2\\mu}\\,\\Phi(\\varphi), \\qquad \\alpha \\approx 0.5445, $
%[text] in polar coordinates $(r,\\varphi)$, with an angular function $\\Phi$ that depends on $\\lambda$ and $\\mu$; $u\_D$ is its trace on $\\Gamma\_D$. The exponent $\\alpha$ is the smallest positive root of $\\sin(3\\pi\\alpha/2) = \\alpha$ and does not depend on the material. The stress behaves like $r^{\\alpha-1}$ and is unbounded at the corner: $u$ belongs to $H^{1+s}(\\Omega)^2$ only for $s \< \\alpha$.
%[text] ### Weak formulation
%[text] Multiplication with a test function $v$ that vanishes on $\\Gamma\_D$, integration by parts and the symmetry of $\\sigma$, $\\sigma:\\nabla v = \\sigma:\\varepsilon(v)$, give: find $u\\in H^1(\\Omega)^2$ with $u = u\_D$ on $\\Gamma\_D$ and
%[text]{"align":"center"} $ \\int\_\\Omega \\mathbb{C}\\varepsilon(u):\\varepsilon(v)\\,dx = \\int\_\\Omega f\\cdot v\\,dx + \\int\_{\\Gamma\_N} g\\cdot v\\,ds \\quad \\textrm{for all } v\\in H^1\_D(\\Omega)^2. $
%[text] Korn's inequality bounds $\\|\\nabla v\\|\_{L^2(\\Omega)}$ by $\\|\\varepsilon(v)\\|\_{L^2(\\Omega)}$ for all $v\\in H^1\_D(\\Omega)^2$ when $\\Gamma\_D$ has positive length, so the bilinear form is coercive and the Lax-Milgram lemma gives a unique solution.
%[text] ### Discretization
%[text] The discrete space $V\_h = P\_k(\\mathcal{T})^2\\cap C(\\overline\\Omega)^2$ (`PkVector(degree,mesh)`) has two copies of the Lagrange space `Pk(degree,mesh)` of `exPoissonPkAdaptive`, one per component, numbered component after component. The degrees of freedom on $\\Gamma\_D$ (`fe.fixedNodes`) take the values of the nodal interpolation `fe.interpolate(@p.u4Db)` of $u\_D$, the others (`fe.freeNodes`) are the unknowns. The coefficient `C(x)` is the constant matrix `p.C`.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega \\mathbb{C}\\varepsilon(u):\\varepsilon(v)\\,dx$ | `C(x)*eps(u):eps(v) dx` |
%[text] | $\\int\_\\Omega f\\cdot v\\,dx + \\int\_{\\Gamma\_N} g\\cdot v\\,ds$ | `f(x)*v dx + g(x)*v dsn` |
%[text] | $\\Vert u-u\_h\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - u )^2 dx` |
%[text] | $\\Vert \\sigma-\\sigma\_h\\Vert\_{\\mathbb{C}^{-1}}$ | `( exact(x) - C(x)*eps(u) )^2 Cin(x) dx` |
%[text:table]
%[text] The last line measures the error of the discrete stress $\\sigma\_h = \\mathbb{C}\\varepsilon(u\_h)$ in the norm $\\|\\tau\\|\_{\\mathbb{C}^{-1}}^2 = \\int\_\\Omega \\tau:\\mathbb{C}^{-1}\\tau\\,dx$ (the weight `Cin(x)` is `p.CInv`). Since $\\sigma-\\sigma\_h = \\mathbb{C}\\varepsilon(u-u\_h)$, it equals the energy norm of the displacement error, $\\|\\sigma-\\sigma\_h\\|\_{\\mathbb{C}^{-1}}^2 = \\int\_\\Omega \\mathbb{C}\\varepsilon(u-u\_h):\\varepsilon(u-u\_h)\\,dx$.
%[text] ### Error estimator
%[text] The residual-based error estimator of `exPoissonPkAdaptive` carries over with the stress $\\mathbb{C}\\varepsilon(u\_h)$ in place of the gradient: it collects the residuals of the equilibrium equation on the triangles, the jumps of the tractions $\\mathbb{C}\\varepsilon(u\_h)\\,n\_E$ across the interior edges $\\mathcal{E}(\\Omega)$ and the residuals of the traction condition on the edges $\\mathcal{E}(\\Gamma\_N)$:
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{T\\in\\mathcal{T}} h\_T^2 \\|f + \\mathrm{div}\\,\\mathbb{C}\\varepsilon(u\_h)\\|^2\_{L^2(T)} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} h\_E \\|\[\\mathbb{C}\\varepsilon(u\_h)\\,n\_E\]\\|^2\_{L^2(E)} + \\sum\_{E\\in\\mathcal{E}(\\Gamma\_N)} h\_E \\|g - \\mathbb{C}\\varepsilon(u\_h)\\,n\\|^2\_{L^2(E)}. $
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2\\Vert f+\\mathrm{div}\\,\\mathbb{C}\\varepsilon(u\_h)\\Vert^2\_{L^2(T)}$ | `h^2*(f(x)+div(C(x)*eps(u)))^2 dx` |
%[text] | $h\_E\\Vert\[\\mathbb{C}\\varepsilon(u\_h)\\,n\_E\]\\Vert^2\_{L^2(E)}$ | `h^1*([C(x)*eps(u)*n])^2 dsi` |
%[text] | $h\_E\\Vert g-\\mathbb{C}\\varepsilon(u\_h)\\,n\\Vert^2\_{L^2(E)}$ | `h^1*(g(x)-C(x)*eps(u)*n)^2 dsn` |
%[text:table]
%[text] The estimator is reliable and efficient up to data oscillations, with constants that depend on $\\mu$ and $\\lambda$. The adaptive loop SOLVE, ESTIMATE, MARK (bulk marking `MarkBulk` with $\\theta = 1/2$) and REFINE is that of `exPoissonPkAdaptive`.
%[text] ### Expected convergence
%[text] On uniformly refined meshes the singularity limits the energy error to $O(h^\\alpha) = O(N^{-\\alpha/2})\\approx O(N^{-0.27})$ for every degree $k$. Adaptive refinement concentrates the triangles at the corner and restores the rate $N^{-k/2}$ of a smooth solution; the estimator decreases with the same rate and the $L^2$ error of the displacement faster.
%[text] ### Locking
%[text] For nearly incompressible material, $\\nu\\to 1/2$, the parameter $\\lambda$ tends to infinity while $\\mu$ stays bounded. The term $\\lambda\\|\\mathrm{div}\\,u\_h\\|^2\_{L^2(\\Omega)}$ of the energy then forces $u\_h$ to be almost divergence free, which leaves too few functions in low-order conforming spaces: the constant of the error estimate grows with $\\lambda$, and for $P\_1$ the discrete displacement is far too small on all practical meshes (volumetric locking). Mixed methods such as `exElasticityAW` are robust in this limit.
%%
%[text] ## Parameters
degree = 2;       % polynomial degree k
minNrDoF = 1e4;   % stop as soon as the number of degrees of freedom exceeds minNrDoF
mark = MarkBulk;  % bulk (Doerfler) marking with theta = 1/2
% mark = MarkUniform;  % uniform refinement, for comparison
E = 1e5;          % Young's modulus
nu = 0.3;         % Poisson ratio
%%
%[text] ## Problem data and initial mesh
p = problemElasticityLshapeExact('E',E,'nu',nu);
mesh = Mesh(p.geom);
plot(mesh);
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the loop stops when the number of degrees of freedom exceeds `minNrDoF`.
nrDoF4lvl = [];
L2Error4lvl = [];
energyError4lvl = [];
eta4lvl = [];
while true
    % finite element space
    fe = PkVector(degree,mesh);
    nrDoF4lvl(end+1) = fe.nrDoF;
    % assemble
    a = BilinearForm('C(x)*eps(u):eps(v) dx');
    A = a.assemble(fe,'C',p.C);
    l = LinearForm('f(x)*v dx + g(x)*v dsn');
    b = l.assemble(fe,'f',@p.f,'g',@(x,n)p.g(x,n));
    % solve
    x = zeros(fe.nrDoF,1);
    u4Db = fe.interpolate(@p.u4Db);
    x(fe.fixedNodes) = u4Db(fe.fixedNodes);
    b = b - A*x;
    free = fe.freeNodes;
    x(free) = A(free,free)\b(free);
    % estimate
    eta = EstimatorForm(['h^2*(f(x)+div(C(x)*eps(u)))^2 dx', ...
        '+ h^1*([C(x)*eps(u)*n])^2 dsi', ...
        '+ h^1*(g(x)-C(x)*eps(u)*n)^2 dsn']);
    eta4e = eta.assemble(x,fe,'f',@p.f,'g',@(x,n)p.g(x,n),'C',p.C);
    eta4lvl(end+1) = norm(eta4e);
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(end+1) = l2.assemble(x,fe,'exact',@p.exactU);
    energy = ErrorForm('( exact(x) - C(x)*eps(u) )^2 Cin(x) dx');
    energyError4lvl(end+1) = energy.assemble(x,fe,'C',p.C,'exact',@p.exactSigma,'Cin',p.CInv);
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
%[text] The last column is the ratio of the energy error and the estimator.
table(nrDoF4lvl',L2Error4lvl',energyError4lvl',eta4lvl',(energyError4lvl./eta4lvl)', ...
    VariableNames=["N" "L2 error" "energy error" "eta" "error/eta"])
%[text] The ratio stays roughly constant: the estimator reflects the error on every mesh. It is small because $\\eta$ measures residuals of the stress, while the energy norm $\\|\\sigma-\\sigma\_h\\|\_{\\mathbb{C}^{-1}}$ is of the size $\\|\\sigma-\\sigma\_h\\|\_{L^2(\\Omega)}/\\sqrt{2\\mu}$ and $\\mu$ is large for $E = 10^5$.
%%
%[text] ### Adaptive mesh
clf
plot(mesh);
%[text] The refinement concentrates at the reentrant corner, where the stress is singular.
%%
%[text] ### Discrete solution
%[text] The color shows the von Mises stress.
clf
feSigma = DGPkSymMatrix(fe.degree,mesh);
sigma = feSigma.interpolate(@(pts,pts_ref,parts)matMul(permute(p.C,[3 1 2]),fe.getSymGradValue(pts,pts_ref,parts,x(1:fe.nrDoF))));
fe.plotDisplacement(x,5e3,sigma,feSigma);
clim([0,8]);
colorbar;
shading interp;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,energyError4lvl,'||\sigma-\sigma_h||_{C^{-1}}');
plotConvergence(nrDoF4lvl,eta4lvl,'\eta');
hold off
%[text] The energy error and the estimator decrease with the same rate. On the finer meshes it approaches $k/2$, the rate of a smooth solution; the least-squares slopes in the legend include the coarse meshes and are somewhat smaller, but far above the rate $\\alpha/2\\approx 0.27$ of uniform refinement. The $L^2$ error of the displacement decreases faster, close to $(k+1)/2$ on the finer meshes.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkUniform;` to refine uniformly: the rate of the energy error drops to $N^{-\\alpha/2}$.
%[text] - Set `degree = 1` or `3`: the adaptive rate of the energy error is $N^{-k/2}$; for $k = 3$ it shows only on finer meshes, so raise `minNrDoF`.
%[text] - Observe locking: use `p = problemElasticitySquareExactLocking('E',E,'nu',nu);`, whose stress and body force do not depend on $\\lambda$, with `degree = 1`, `mark = MarkUniform;` and `nu = 0.4999`. The relative $L^2$ error of the displacement stays close to one up to a few thousand degrees of freedom and decreases only slowly afterwards (to about 0.4 on the last mesh), while it decreases with the rate $N^{-1}$ to about $10^{-3}$ for `nu = 0.3`. For `degree = 2` the errors converge, but with a lower rate and more than an order of magnitude larger than for `nu = 0.3`. `exElasticityAW` solves the same problem with a mixed method that is robust in $\\lambda$.
%[text] - Compute the deformation of Cook's membrane with `p = problemElasticityCooks('E',E,'nu',nu);`. This problem has no exact solution: remove the errors step, plot the estimator alone and show the deformed mesh with `fe.plotDisplacement(x,p.scale);`.
%[text] - Use the smooth solution of `problemElasticitySquareExact('E',E,'nu',nu)` on the unit square: uniform refinement then converges with the rate $N^{-k/2}$ of the energy error, and adaptive refinement gains little. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
