%[text] # Arnold-Winther mixed finite elements for linear elasticity
%[text] Approximates the stress and the displacement of a nearly incompressible elastic body simultaneously with the Arnold-Winther elements, whose stresses are symmetric and have continuous normal components across the edges. In contrast to the displacement method of `exElasticityPkAdaptive` the errors do not deteriorate as the Poisson ratio approaches 1/2.
%[text] **Builds on:** `exPoissonRTAdaptive`, `exElasticityPkAdaptive`. **Next:** `exBiharmonicHHJ`, `exStokesTaylorHood`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] As in `exElasticityPkAdaptive`, the displacement $u$ of a linear elastic body and its stress $\\sigma = \\mathbb{C}\\varepsilon(u)$ satisfy the equilibrium equation $-\\mathrm{div}\\,\\sigma = f$. The Hellinger-Reissner formulation keeps the stress as an unknown and writes Hooke's law with the compliance tensor $\\mathbb{C}^{-1}$:
%[text]{"align":"center"} $ \\mathbb{C}^{-1}\\sigma = \\varepsilon(u), \\qquad -\\mathrm{div}\\,\\sigma = f \\textrm{ in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\Gamma\_D, \\qquad \\sigma n = g \\textrm{ on } \\Gamma\_N. $
%[text] In two dimensions the compliance tensor is
%[text]{"align":"center"} $ \\mathbb{C}^{-1}\\sigma = \\frac{1}{2\\mu}\\left(\\sigma - \\frac{\\lambda}{2(\\lambda+\\mu)}\\,\\mathrm{tr}(\\sigma)\\,I\\right) \\quad\\longrightarrow\\quad \\frac{1}{2\\mu}\\left(\\sigma - \\frac{1}{2}\\mathrm{tr}(\\sigma)\\,I\\right) \\quad \\textrm{for } \\lambda\\to\\infty. $
%[text] Unlike $\\mathbb{C}$, it stays bounded in the incompressible limit $\\nu\\to 1/2$, where $\\lambda = E\\nu/((1+\\nu)(1-2\\nu))\\to\\infty$; its matrix on the components $\[\\sigma\_{11}, \\sigma\_{22}, \\sigma\_{12}\]$ is `p.CInv` (see `ElasticityParameters`).
%[text] The problem object `problemElasticitySquareExactLocking` prescribes $u\_D = 0$ on the whole boundary of the unit square, $\\Gamma\_D = \\partial\\Omega$, and the body force $f$ of the divergence-free exact displacement
%[text]{"align":"center"} $ u = \\pi\\left(\\begin{array}{c} \\cos(\\pi x\_2)\\sin^2(\\pi x\_1)\\sin(\\pi x\_2) \\\\ -\\cos(\\pi x\_1)\\sin(\\pi x\_1)\\sin^2(\\pi x\_2) \\end{array}\\right). $
%[text] Since $\\mathrm{div}\\,u = 0$, the stress $\\sigma = 2\\mu\\,\\varepsilon(u)$ and the body force do not depend on $\\lambda$, which makes the problem a test for locking.
%[text] ### Weak formulation
%[text] Multiplication of $\\mathbb{C}^{-1}\\sigma = \\varepsilon(u)$ with a symmetric test function $\\tau$ and the integration by parts
%[text]{"align":"center"} $ \\int\_\\Omega \\varepsilon(u):\\tau\\,dx = -\\int\_\\Omega u\\cdot\\mathrm{div}\\,\\tau\\,dx + \\int\_{\\partial\\Omega} u\\cdot\\tau n\\,ds $
%[text] give the Hellinger-Reissner formulation: find $\\sigma\\in H(\\mathrm{div},\\Omega;\\mathbb{S})$ with $\\sigma n = g$ on $\\Gamma\_N$ and $u\\in L^2(\\Omega)^2$ with
%[text]{"align":"center"} $ \\int\_\\Omega \\mathbb{C}^{-1}\\sigma:\\tau\\,dx + \\int\_\\Omega u\\cdot\\mathrm{div}\\,\\tau\\,dx = \\int\_{\\Gamma\_D} u\_D\\cdot\\tau n\\,ds, $
%[text]{"align":"center"} $ \\int\_\\Omega \\mathrm{div}\\,\\sigma\\cdot v\\,dx = -\\int\_\\Omega f\\cdot v\\,dx $
%[text] for all $\\tau\\in H(\\mathrm{div},\\Omega;\\mathbb{S})$ with $\\tau n = 0$ on $\\Gamma\_N$ and all $v\\in L^2(\\Omega)^2$. Here $H(\\mathrm{div},\\Omega;\\mathbb{S})$ consists of the symmetric matrix fields with square-integrable rows and row-wise divergence. The roles of the boundary conditions are exchanged compared with the displacement method: the displacement condition is natural and enters the right-hand side, the traction condition is essential and restricts the stress space. The displacement needs no derivatives and lies in $L^2$ only.
%[text] ### Discretization
%[text] A stable discretization needs a stress space of symmetric, $H(\\mathrm{div})$-conforming matrix fields, i.e. with continuous normal components $\\sigma n\_E$ across the edges, whose divergence exactly fills the displacement space. Arnold and Winther constructed the first such pair with polynomials: for $k\\ge 1$
%[text]{"align":"center"} $ \\Sigma\_h = \\{\\tau\\in H(\\mathrm{div},\\Omega;\\mathbb{S}) : \\tau|\_T\\in P\_{k+2}(T;\\mathbb{S}),\\ \\mathrm{div}\\,\\tau|\_T\\in P\_k(T)^2 \\textrm{ for all } T\\in\\mathcal{T},\\ \\tau \\textrm{ continuous at the vertices}\\}, \\qquad V\_h = P\_k(\\mathcal{T})^2, $
%[text] with discontinuous displacements (`AWk(degree+2,mesh)` and `DGDkVector(degree,mesh)`). The degrees of freedom of $\\Sigma\_h$ are the values of $\\sigma$ at the vertices (3 per vertex), shared by the triangles of a vertex, the moments of degree $0,\\dots,k$ of $\\sigma n$ on the edges ($2k+2$ per edge) and interior moments (3 per triangle for $k = 1$); for $k = 1$ a triangle carries 24 degrees of freedom in total.
%[text] The Dirichlet data enter the right-hand side $b\_D$. On $\\Gamma\_N$ the traction is imposed on the stress: `feSigma.interpolateBoundary(@p.g)` sets the edge moments of $\\sigma\_h n$ on $\\Gamma\_N$ (`feSigma.fixedNodes`) and, at a corner between two non-parallel Neumann edges, the vertex value of $\\sigma\_h$. At the other vertices of $\\Gamma\_N$ the condition $\\sigma\_h(z) n = g(z)$ is a linear constraint $C\\sigma = b\_N$ (`feSigma.Constraints` and the second output of `interpolateBoundary`), imposed with a Lagrange multiplier. The coefficient vectors $\\sigma$, $u$ and the multiplier solve
%[text]{"align":"center"} $ \\left(\\begin{array}{ccc} A & B^T & C^T \\\\ B & 0 & 0 \\\\ C & 0 & 0 \\end{array}\\right) \\left(\\begin{array}{c} \\sigma \\\\ u \\\\ \\theta \\end{array}\\right) = \\left(\\begin{array}{c} b\_D \\\\ -b \\\\ b\_N \\end{array}\\right). $
%[text] The square of this example has $\\Gamma\_N = \\emptyset$, so the traction rows are empty and all stress degrees of freedom are free.
%[text:table]
%[text] | Matrix | Integral | fem4m form |
%[text] | --- | --- | --- |
%[text] | $A$ | $\\int\_\\Omega \\mathbb{C}^{-1}\\sigma:\\tau\\,dx$ | `Cinv(x)*sigma*tau dx` |
%[text] | $B$ | $\\int\_\\Omega \\mathrm{div}\\,\\sigma\\cdot v\\,dx$ | `div(sigma)*v dx` |
%[text] | $b\_D$ | $\\int\_{\\Gamma\_D} u\_D\\cdot\\tau n\\,ds$ | `uD(x)*(tau*n) dsd` |
%[text] | $b$ | $\\int\_\\Omega f\\cdot v\\,dx$ | `f(x)*v dx` |
%[text] |  | $\\Vert u-u\_h\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - u )^2 dx` |
%[text] |  | $\\Vert \\sigma-\\sigma\_h\\Vert\_{\\mathbb{C}^{-1}}$ | `( exact(x) - sigma )^2 Cin(x) dx` |
%[text:table]
%[text] The stress error is measured in the norm $\\|\\tau\\|\_{\\mathbb{C}^{-1}}^2 = \\int\_\\Omega \\mathbb{C}^{-1}\\tau:\\tau\\,dx$ of `exElasticityPkAdaptive`.
%[text] ### Expected convergence
%[text] The discrete inf-sup condition holds with a constant independent of the mesh size, and $\\mathbb{C}^{-1}$ is bounded uniformly in $\\lambda$. For a smooth solution
%[text]{"align":"center"} $ \\|\\sigma-\\sigma\_h\\|\_{L^2(\\Omega)} \\le C\\,h^{k+2}, \\qquad \\|u-u\_h\\|\_{L^2(\\Omega)} \\le C\\,h^{k+1}, $
%[text] with constants that do not depend on $\\lambda$. On uniformly refined meshes $N\\approx h^{-2}$, so the stress error decays like $N^{-(k+2)/2}$ and the displacement error like $N^{-(k+1)/2}$, also for $\\nu$ close to $1/2$.
%%
%[text] ## Parameters
degree = 1;     % degree k of the displacement, the stress has degree k+2
stressSpace = @(mesh) AWk(degree+2,mesh);              % Arnold-Winther stresses
displacementSpace = @(mesh) DGDkVector(degree,mesh);   % discontinuous P_k displacements
% stressSpace = @(mesh) AW0(3,mesh); displacementSpace = @(mesh) RBM(1,mesh);  % reduced element, piecewise rigid motions
E = 1e5;        % Young's modulus
nu = 0.4999;    % Poisson ratio, close to the incompressible limit 1/2
nrLevels = 4;   % number of meshes
%%
%[text] ## Problem data and initial mesh
%[text] The coarse mesh of the square has two triangles; one uniform refinement by newest vertex bisection gives the first mesh.
p = problemElasticitySquareExactLocking('E',E,'nu',nu);
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform','nvb');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
nrDoF4lvl = zeros(1,nrLevels);
L2Error4lvl = zeros(1,nrLevels);
sigmaError4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element spaces
    feSigma = stressSpace(mesh);
    feU = displacementSpace(mesh);
    nrDoF4lvl(lvl) = feSigma.nrDoF + feU.nrDoF;
    % assemble
    a = BilinearForm('Cinv(x)*sigma*tau dx',{'sigma','tau'});
    A = a.assemble(feSigma,'Cinv',p.CInv);
    bDiv = BilinearForm('div(sigma)*v dx',{'sigma','v'});
    B = bDiv.assemble(feSigma,feU);
    lD = LinearForm('uD(x)*(tau*n) dsd',{'sigma','tau'});
    bD = lD.assemble(feSigma,'uD',@p.u4Db);
    l = LinearForm('f(x)*v dx');
    b = l.assemble(feU,'f',@p.f);
    % solve: the traction as essential boundary condition, the vertex constraints with a multiplier
    [x,bN] = feSigma.interpolateBoundary(@p.g);
    C = feSigma.Constraints;
    nrC = size(C,1);
    S = [A, B', C'; B, sparse(feU.nrDoF,feU.nrDoF+nrC); C, sparse(nrC,feU.nrDoF+nrC)];
    x = [x; zeros(feU.nrDoF+nrC,1)];
    rhs = [bD; -b; bN] - S*x;
    free = [feSigma.freeNodes, feSigma.nrDoF+1:size(S,1)];
    x(free) = S(free,free)\rhs(free);
    sigma = x(1:feSigma.nrDoF);
    u = x(feSigma.nrDoF+(1:feU.nrDoF));
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(lvl) = l2.assemble(u,feU,12,'exact',@p.exactU);
    l2Sigma = ErrorForm('( exact(x) - sigma )^2 Cin(x) dx','sigma');
    sigmaError4lvl(lvl) = l2Sigma.assemble(sigma,feSigma,12,'exact',@p.exactSigma,'Cin',p.CInv);
    % refine
    mesh = mesh.refine('uniform','nvb');
end
%%
%[text] ## Results
%[text] ### Convergence table
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
sigmaRate = [NaN, -diff(log(sigmaError4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',sigmaError4lvl',sigmaRate', ...
    VariableNames=["N" "u error" "u rate" "sigma error" "sigma rate"])
%[text] The rates approach $(k+1)/2$ for the displacement and $(k+2)/2$ for the stress, although $\\nu$ is close to $1/2$.
%%
%[text] ### Discrete solution
%[text] The color shows the piecewise polynomial von Mises stress on the finest mesh.
clf
feU.plotDisplacement(u,1e-1,sigma,feSigma);
colorbar;
shading interp;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,sigmaError4lvl,'||\sigma-\sigma_h||_{C^{-1}}');
hold off
%[text] Both errors converge with the rates of a smooth solution. The conforming $P\_1$ discretization of the same problem locks: its relative error stays close to one on these meshes (see the extensions of `exElasticityPkAdaptive`).
%%
%[text] ## Possibilities for extensions
%[text] - Set `nu = 0.3`: the relative errors of the displacement agree with those for `nu = 0.4999` to three or more digits, those of the stress within about 10 percent, and the rates are the same; the method is robust in $\\lambda$.
%[text] - Set `degree = 2` for the stresses `AWk(4,mesh)` of degree 4 and `nrLevels = 5`: on the finer meshes the rates grow to $N^{-3/2}$ for the displacement and $N^{-2}$ for the stress.
%[text] - Use the reduced Arnold-Winther element, the commented line in the parameters: `AW0(3,mesh)` has the cubic stresses whose divergence is a rigid motion, without the interior degrees of freedom, and `RBM(1,mesh)` the piecewise rigid motions $a + b\\,(-x\_2, x\_1)$ as displacements. It converges with the rates $N^{-1}$ for the stress and $N^{-1/2}$ for the displacement. `AW1(3,mesh)` is another basis of the space of `AWk(3,mesh)` and gives the same results.
%[text] - Impose tractions: `problemElasticitySquareNeumannExact` prescribes $g$ on three edges of the square. With `nu = 0.3` the stress error converges with the rate $3/2$; without the constraint rows `C` the vertex values on $\\Gamma\_N$ are not tied to $g$, and the rate drops to about $3/4$.
%[text] - The application `applications/elasticity/afemElasticitAW.m` can also refine adaptively (argument `refine = 'adaptive'`) with a residual-based error estimator for the stress. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
