%[text] # Taylor-Hood finite elements for the Stokes problem
%[text] Computes the velocity and the pressure of a Stokes flow on the unit square with the Taylor-Hood elements, continuous piecewise polynomials of degree $k$ for the velocity and $k-1$ for the pressure. The example introduces the saddle-point structure of incompressible flow, the inf-sup condition that a stable pair of spaces must satisfy, and the mean-value constraint that fixes the pressure; the Mini element, the nonconforming Crouzeix-Raviart element and the unstable $P\_1\\times P\_0$ pair are alternatives in the parameters.
%[text] **Builds on:** `exPoissonMixedP1P0Unstable`, `exPoissonRTAdaptive`. **Next:** `exStokesDG`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] The slow, viscous flow of an incompressible fluid with viscosity $\\nu = 1$ through $\\Omega = (0,1)^2$ is described by the velocity $u:\\Omega\\to\\mathbb{R}^2$ and the pressure $p:\\Omega\\to\\mathbb{R}$. Given the body force $f$ and the velocity $u\_D$ on the boundary, they satisfy the Stokes equations
%[text]{"align":"center"} $ -\\Delta u + \\nabla p = f \\textrm{ in } \\Omega, \\qquad \\mathrm{div}\\,u = 0 \\textrm{ in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\partial\\Omega. $
%[text] The pressure acts as the Lagrange multiplier of the incompressibility constraint $\\mathrm{div}\\,u = 0$. Only its gradient appears, so it is unique only up to a constant, which is fixed by the condition $\\int\_\\Omega p\\,dx = 0$. The divergence theorem requires the compatibility condition $\\int\_{\\partial\\Omega} u\_D\\cdot n\\,ds = 0$.
%[text] The problem object `p = problemStokesExactSquareSinus` provides the data of the exact solution with $s = \\sin(\\pi(x\_1-x\_2))$,
%[text]{"align":"center"} $ u = (s, s), \\qquad p = -2\\pi\\cos(\\pi(x\_1-x\_2)) + 8/\\pi, $
%[text] whose velocity does not vanish on the boundary, so $u\_D = u|\_{\\partial\\Omega}$ are inhomogeneous Dirichlet data; the pressure has mean zero. In the code `p` is the problem object, the pressure coefficients are called `pressure`.
%[text] ### Weak formulation
%[text] With $L^2\_0(\\Omega) = \\{q\\in L^2(\\Omega) : \\int\_\\Omega q\\,dx = 0\\}$ the weak formulation reads: find $u\\in H^1(\\Omega)^2$ with $u = u\_D$ on $\\partial\\Omega$ and $p\\in L^2\_0(\\Omega)$ with
%[text]{"align":"center"} $ \\int\_\\Omega \\nabla u:\\nabla v\\,dx - \\int\_\\Omega p\\,\\mathrm{div}\\,v\\,dx = \\int\_\\Omega f\\cdot v\\,dx \\quad \\textrm{for all } v\\in H^1\_0(\\Omega)^2, $
%[text]{"align":"center"} $ -\\int\_\\Omega q\\,\\mathrm{div}\\,u\\,dx = 0 \\quad \\textrm{for all } q\\in L^2\_0(\\Omega). $
%[text] This is a saddle-point problem. It is well posed because of the inf-sup (Ladyzhenskaya-Babuška-Brezzi) condition: there is a $\\beta \> 0$ with
%[text]{"align":"center"} $ \\sup\_{v\\in H^1\_0(\\Omega)^2\\setminus\\{0\\}} \\frac{\\int\_\\Omega q\\,\\mathrm{div}\\,v\\,dx}{\\|\\nabla v\\|\_{L^2(\\Omega)}} \\ge \\beta\\,\\|q\\|\_{L^2(\\Omega)} \\quad \\textrm{for all } q\\in L^2\_0(\\Omega). $
%[text] ### Discretization
%[text] A pair of discrete spaces $V\_h\\times Q\_h$ is stable if it satisfies the inf-sup condition with a constant $\\beta$ independent of the mesh size. The velocity space must then be rich enough compared with the pressure space: if $Q\_h$ is too large, the discrete divergence-free velocities are too few to approximate $u$ (locking), or spurious pressure modes $q\_h\\ne 0$ with $\\int\_\\Omega q\_h\\,\\mathrm{div}\\,v\_h\\,dx = 0$ for all $v\_h$ make the pressure non-unique. The example offers four pairs:
%[text:table]
%[text] | Pair | $V\_h$ | $Q\_h$ | fem4m spaces | stable |
%[text] | --- | --- | --- | --- | --- |
%[text] | Taylor-Hood | $P\_k(\\mathcal{T})^2\\cap C(\\overline\\Omega)^2$ | $P\_{k-1}(\\mathcal{T})\\cap C(\\overline\\Omega)$ | `PkVector(degree,mesh)`, `Pk(degree-1,mesh)` | yes, $k\\ge 2$ |
%[text] | Mini | $P\_1$ plus cubic bubbles | $P\_1(\\mathcal{T})\\cap C(\\overline\\Omega)$ | `PkVector(1,mesh) + FEVector(Bk(3,mesh),Bk(3,mesh))`, `Pk(1,mesh)` | yes |
%[text] | Crouzeix-Raviart | $CR\_1(\\mathcal{T})^2$ | $P\_0(\\mathcal{T})$ | `FEVector(CRk(1,mesh),CRk(1,mesh))`, `DGDk(0,mesh)` | yes |
%[text] | $P\_1\\times P\_0$ | $P\_1(\\mathcal{T})^2\\cap C(\\overline\\Omega)^2$ | $P\_0(\\mathcal{T})$ | `PkVector(1,mesh)`, `DGDk(0,mesh)` | no |
%[text:table]
%[text] The Crouzeix-Raviart functions are piecewise linear and continuous only at the edge midpoints; they are nonconforming, and the gradients in the bilinear form are taken on every triangle. The Taylor-Hood pairs satisfy the inf-sup condition on every triangulation with at least three triangles \[1, Theorem 8.8.1\]. The example replaces the default mesh of the square, two triangles, by the criss-cross mesh, four triangles that meet at the midpoint; on its refinements the unstable pair $P\_1\\times P\_0$ shows the typical spurious pressure modes (see the extensions), while the Taylor-Hood pairs converge with the same rates on both meshes.
%[text] The velocity degrees of freedom on $\\partial\\Omega$ (`feU.fixedNodes`) take the values of $u\_D$ at the Lagrange nodes or, for Crouzeix-Raviart, at the edge midpoints. The mean-value condition $\\int\_\\Omega p\_h\\,dx = 0$ enters with a scalar Lagrange multiplier $\\theta$, which keeps the matrix symmetric. The coefficient vectors $u$ of the free velocity degrees of freedom, $p$ of the pressure and $\\theta$ solve the linear system
%[text]{"align":"center"} $ \\left(\\begin{array}{ccc} A & B^T & 0 \\\\ B & 0 & c \\\\ 0 & c^T & 0 \\end{array}\\right) \\left(\\begin{array}{c} u \\\\ p \\\\ \\theta \\end{array}\\right) = \\left(\\begin{array}{c} b \\\\ 0 \\\\ 0 \\end{array}\\right), $
%[text] where the contributions of the Dirichlet values are moved to the right-hand side (`rhs = ... - S*x` in the code).
%[text] The last row is $\\int\_\\Omega p\_h\\,dx = 0$. The test function $q = 1$ in the second row gives $\\theta\\,|\\Omega| = \\int\_\\Omega \\mathrm{div}\\,u\_h\\,dx$, so $\\theta$ vanishes when the discrete boundary data have zero flux.
%[text:table]
%[text] | Matrix | Integral | fem4m form |
%[text] | --- | --- | --- |
%[text] | $A$ | $\\int\_\\Omega \\nabla u:\\nabla v\\,dx$ | `grad(u)*grad(v) dx` |
%[text] | $B$ | $-\\int\_\\Omega q\\,\\mathrm{div}\\,u\\,dx$ | `-div(u)*v dx` |
%[text] | $c$ | $\\int\_\\Omega q\\,dx$ | `v dx` |
%[text] | $b$ | $\\int\_\\Omega f\\cdot v\\,dx$ | `f(x)*v dx` |
%[text] |  | $\\Vert u-u\_h\\Vert\_{L^2(\\Omega)}$, $\\Vert p-p\_h\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - u )^2 dx` |
%[text] |  | $\\Vert \\nabla(u-u\_h)\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - grad(u) )^2 dx` |
%[text:table]
%[text] In the form for $B$ the test function of the pressure space is called `v`.
%[text] ### Expected convergence
%[text] For a stable conforming pair the discrete solution is quasi-optimal,
%[text]{"align":"center"} $ \\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)} + \\|p-p\_h\\|\_{L^2(\\Omega)}\\le C\\inf\_{v\_h\\in V\_h,\\, q\_h\\in Q\_h}\\left(\\|\\nabla(u-v\_h)\\|\_{L^2(\\Omega)} + \\|p-q\_h\\|\_{L^2(\\Omega)}\\right), $
%[text] and a duality argument gives one more order for $\\|u-u\_h\\|\_{L^2(\\Omega)}$ on the convex square. For the nonconforming Crouzeix-Raviart pair the gradient is the piecewise one and a consistency error adds to the right-hand side, of the same order. For the smooth solution of this example and $N\\approx h^{-2}$:
%[text:table]
%[text] | Pair | $\\Vert\\nabla(u-u\_h)\\Vert$ and $\\Vert p-p\_h\\Vert$ | $\\Vert u-u\_h\\Vert$ |
%[text] | --- | --- | --- |
%[text] | Taylor-Hood, degree $k$ | $N^{-k/2}$ | $N^{-(k+1)/2}$ |
%[text] | Mini | $N^{-1/2}$ | $N^{-1}$ |
%[text] | Crouzeix-Raviart | $N^{-1/2}$ | $N^{-1}$ |
%[text:table]
%%
%[text] ## Parameters
%[text] The two function handles build the velocity and the pressure space on a mesh; the commented lines choose another pair.
degree = 2;                                    % Taylor-Hood degree k >= 2
velocitySpace = @(mesh) PkVector(degree,mesh);   % continuous P_k velocity
pressureSpace = @(mesh) Pk(degree-1,mesh);       % continuous P_{k-1} pressure
% velocitySpace = @(mesh) PkVector(1,mesh) + FEVector(Bk(3,mesh),Bk(3,mesh)); pressureSpace = @(mesh) Pk(1,mesh);  % Mini element
% velocitySpace = @(mesh) FEVector(CRk(1,mesh),CRk(1,mesh)); pressureSpace = @(mesh) DGDk(0,mesh);  % nonconforming Crouzeix-Raviart element
% velocitySpace = @(mesh) PkVector(1,mesh); pressureSpace = @(mesh) DGDk(0,mesh);  % unstable P_1 x P_0
nrLevels = 5;                                  % number of meshes
%%
%[text] ## Problem data and initial mesh
%[text] The criss-cross mesh of the square replaces the coarse mesh of the problem object.
p = problemStokesExactSquareSinus;
p.geom = Geometry(struct('name','Square','reftype','crisscross'));
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
nrDoF4lvl = zeros(1,nrLevels);
L2Error4lvl = zeros(1,nrLevels);
H1Error4lvl = zeros(1,nrLevels);
pressureError4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element spaces
    feU = velocitySpace(mesh);
    feP = pressureSpace(mesh);
    nrDoF4lvl(lvl) = feU.nrDoF + feP.nrDoF;
    % assemble
    a = BilinearForm('grad(u)*grad(v) dx');
    A = a.assemble(feU);
    bDiv = BilinearForm('-div(u)*v dx');
    B = bDiv.assemble(feU,feP);
    lMean = LinearForm('v dx');
    c = lMean.assemble(feP);
    l = LinearForm('f(x)*v dx');
    b = l.assemble(feU,'f',@p.f);
    % solve: Dirichlet values at the Lagrange nodes or edge midpoints, mean-value constraint by a multiplier
    x = zeros(feU.nrDoF+feP.nrDoF+1,1);
    if isprop(feU.fe{1},'c4dof')
        x(feU.fixedNodes) = reshape(p.u4Db(feU.fe{1}.c4dof(feU.fe{1}.fixedNodes,:)),[],1);
    else
        mid4ed = mesh.point4ed(1/2);
        x(feU.fixedNodes) = reshape(p.u4Db(mid4ed(feU.fe{1}.fixedNodes,:)),[],1);
    end
    S = [A, B', sparse(feU.nrDoF,1); B, sparse(feP.nrDoF,feP.nrDoF), c; sparse(1,feU.nrDoF), c', 0];
    rhs = [b; zeros(feP.nrDoF+1,1)] - S*x;
    free = [feU.freeNodes(:); feU.nrDoF+(1:feP.nrDoF+1)'];
    x(free) = S(free,free)\rhs(free);
    u = x(1:feU.nrDoF);
    pressure = x(feU.nrDoF+(1:feP.nrDoF));
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(lvl) = l2.assemble(u,feU,6,'exact',@p.exactU);
    pressureError4lvl(lvl) = l2.assemble(pressure,feP,6,'exact',@p.exactP);
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(lvl) = h1.assemble(u,feU,6,'exact',@p.exactGradU);
    % refine
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Convergence table
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
H1Rate = [NaN, -diff(log(H1Error4lvl))./diff(log(nrDoF4lvl))];
pressureRate = [NaN, -diff(log(pressureError4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',H1Error4lvl',H1Rate',pressureError4lvl',pressureRate', ...
    VariableNames=["N" "u L2 error" "u L2 rate" "u H1 error" "u H1 rate" "p error" "p rate"])
%[text] For the Taylor-Hood pair with $k = 2$ the rates approach $3/2$ for the $L^2$ error of the velocity and $1$ for the gradient of the velocity and for the pressure.
%%
%[text] ### Discrete pressure
clf
plot(feP,pressure);
view([40,35]);
colorbar;
%%
%[text] ### Stream lines
%[text] The stream lines of the discrete velocity on the finest mesh, interpolated from its values at the vertices.
clf
feU.plotStreamLines(u);
%%
%[text] ### Direction of the flow
%[text] Unit arrows in the direction of the discrete velocity.
clf
feU.plotVelocity(u);
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(u-u_h)||_{L^2}');
plotConvergence(nrDoF4lvl,pressureError4lvl,'||p-p_h||_{L^2}');
hold off
%[text] The slopes in the legend agree with the expected rates of the chosen pair.
%%
%[text] ## Possibilities for extensions
%[text] - Set `degree = 3` for the Taylor-Hood pair $P\_3\\times P\_2$: the rates grow to $N^{-2}$ for the $L^2$ error of the velocity and $N^{-3/2}$ for the other two errors.
%[text] - Use the Mini element, the first commented line in the parameters: the cubic bubbles enrich the $P\_1$ velocity just enough for stability with the continuous $P\_1$ pressure. The velocity errors converge like $N^{-1}$ and $N^{-1/2}$; the pressure converges faster than the guaranteed $N^{-1/2}$ on these structured meshes.
%[text] - Use the Crouzeix-Raviart element, the second commented line: the rates are $N^{-1}$ for the $L^2$ error of the velocity and $N^{-1/2}$ for the other two errors. The divergence of $u\_h$ is constant on every triangle, and the second block row with $\\theta = 0$ tests it with all piecewise constants: the discrete velocity is exactly divergence free on every triangle.
%[text] - Use the unstable pair $P\_1\\times P\_0$, the third commented line: on criss-cross meshes the velocity still converges, but spurious pressure modes make the matrix singular, MATLAB warns about it, and the pressure error grows under refinement. `exPoissonMixedP1P0Unstable` shows the same kind of instability for the mixed Poisson problem.
%[text] - Compute the lid-driven cavity flow of `problemStokesDrivenCavity`, which has no exact solution: remove the errors step and the convergence plots, keep the stream lines, and remove the line `p.geom = ...` to compute in the cavity $(-1,1)^2$ of the problem object instead of the unit square. `exStokesDG` discretizes the Stokes problem with a discontinuous Galerkin method. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
