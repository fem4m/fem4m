%[text] # $P\_1$ finite elements for the Poisson problem
%[text] The first example: the Poisson problem on the unit square, discretized with continuous piecewise linear finite elements on a sequence of uniformly refined meshes. It shows the basic workflow of fem4m (mesh, finite element space, forms, assembly, solve, errors) and the convergence rates of the method.
%[text] **Builds on:** `exMeshEnumeration`, `exElementBasisFunctions`. **Next:** `exPoissonPkAdaptive`, `exPoissonDG`, `exPoissonMixedP1P0Unstable`, `exSolverSmoothing`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] Let $\\Omega = (0,1)^2$ be the unit square. Given the source $f$, find $u$ with
%[text]{"align":"center"} $ -\\Delta u = f \\quad \\textrm{in } \\Omega, \\qquad u = 0 \\quad \\textrm{on } \\partial\\Omega. $
%[text] The problem object `problemScalarEllipticSquareExact` provides the source $f = 2x\_1(1-x\_1) + 2x\_2(1-x\_2)$ of the exact solution $u = x\_1(1-x\_1)\\,x\_2(1-x\_2)$, so the errors can be computed exactly.
%[text] ### Weak formulation
%[text] Multiplication with a test function $v$ that vanishes on $\\partial\\Omega$ and integration by parts give the weak formulation: find $u\\in H^1\_0(\\Omega)$ with
%[text]{"align":"center"} $ \\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx = \\int\_\\Omega f v\\,dx \\quad \\textrm{for all } v\\in H^1\_0(\\Omega). $
%[text] The bilinear form on the left is the scalar product of $H^1\_0(\\Omega)$, so the Riesz representation theorem gives a unique solution.
%[text] ### Discretization
%[text] Let $\\mathcal{T}$ be a regular triangulation of $\\Omega$ and $V\_h = \\{v\_h\\in C(\\overline\\Omega) : v\_h|\_T\\in P\_1(T) \\textrm{ for all } T\\in\\mathcal{T}\\}$ the continuous piecewise linear functions (`Pk(1,mesh)`). Their nodal basis functions $\\varphi\_1,\\dots,\\varphi\_N$ belong to the $N$ nodes of $\\mathcal{T}$. The boundary condition fixes the coefficients of the boundary nodes (`fe.fixedNodes`) to zero; the coefficients $x\_j$ of the interior nodes (`fe.freeNodes`) of the discrete solution $u\_h = \\sum\_j x\_j\\varphi\_j$ solve the linear system $Ax = b$ with
%[text]{"align":"center"} $ A\_{jk} = \\int\_\\Omega \\nabla\\varphi\_k\\cdot\\nabla\\varphi\_j\\,dx, \\qquad b\_j = \\int\_\\Omega f\\varphi\_j\\,dx. $
%[text] fem4m assembles $A$, $b$ and the errors from form strings that read like the integrals. A name with `(x)`, such as `f(x)`, is a coefficient whose function is passed to `assemble`, e.g. `l.assemble(fe,'f',@p.f)`.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx$ | `BilinearForm('grad(u) * grad(v) dx')` |
%[text] | $\\int\_\\Omega f v\\,dx$ | `LinearForm('f(x) v dx')` |
%[text] | $\\Vert u-u\_h\\Vert\_{L^2(\\Omega)}$ | `ErrorForm('( exact(x) - u )^2 dx')` |
%[text] | $\\Vert \\nabla(u-u\_h)\\Vert\_{L^2(\\Omega)}$ | `ErrorForm('( exact(x) - grad(u) )^2 dx')` |
%[text:table]
%[text] ### Expected convergence
%[text] The solution is smooth, $u\\in H^2(\\Omega)$, and the square is convex, so with the maximal mesh size $h$
%[text]{"align":"center"} $ \\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)} \\le C\\,h\\,|u|\_{H^2(\\Omega)}, \\qquad \\|u-u\_h\\|\_{L^2(\\Omega)} \\le C\\,h^2\\,|u|\_{H^2(\\Omega)}. $
%[text] Every uniform refinement halves $h$ and multiplies the number of degrees of freedom by about four, $N\\approx h^{-2}$. In terms of $N$ the energy error therefore decays like $N^{-1/2}$ and the $L^2$ error like $N^{-1}$.
%%
%[text] ## Parameters
nrLevels = 6;  % number of meshes
%%
%[text] ## Problem data and initial mesh
%[text] The coarse mesh of the square has two triangles; one uniform refinement gives the first mesh with an interior node.
p = problemScalarEllipticSquareExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
%[text] On every mesh the loop assembles and solves the system, computes the errors and refines all triangles uniformly (each triangle into four).
nrDoF4lvl = zeros(1,nrLevels);
L2Error4lvl = zeros(1,nrLevels);
H1Error4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element space
    fe = Pk(1,mesh);
    nrDoF4lvl(lvl) = fe.nrDoF;
    % assemble
    a = BilinearForm('grad(u) * grad(v) dx');
    A = a.assemble(fe);
    l = LinearForm('f(x) v dx');
    b = l.assemble(fe,'f',@p.f);
    % solve
    x = zeros(fe.nrDoF,1);
    free = fe.freeNodes;
    x(free) = A(free,free)\b(free);
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(lvl) = l2.assemble(x,fe,'exact',@p.exactU);
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(lvl) = h1.assemble(x,fe,'exact',@p.exactGradU);
    % refine
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Convergence table
%[text] The rates are the slopes of the errors over $N$ in a log-log scale between two consecutive meshes.
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
H1Rate = [NaN, -diff(log(H1Error4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',H1Error4lvl',H1Rate', ...
    VariableNames=["N" "L2 error" "L2 rate" "H1 error" "H1 rate"])
%[text] The $L^2$ rates approach 1 and the energy rates 1/2, as predicted.
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
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(u-u_h)||_{L^2}');
hold off
%[text] The legend shows the least-squares slopes, the observed rates with respect to $N$.
%%
%[text] ## Possibilities for extensions
%[text] - Use `Pk(k,mesh)` with $k = 2$ or $3$: the energy error decays like $N^{-k/2}$ and the $L^2$ error like $N^{-(k+1)/2}$. For $k\\ge 4$ the polynomial solution lies in $V\_h$ and the errors drop to rounding level.
%[text] - Solve with the conjugate gradient method, `x(free) = pcg(A(free,free),b(free),1e-10,1000);`: the number of iterations grows like $h^{-1}$. `exSolverMultigrid` solves the same systems with a number of iterations independent of $h$.
%[text] - `spy(A)` shows the sparsity of the stiffness matrix: a node couples only with the nodes of its neighbouring triangles.
%[text] - On the L-shaped domain (`problemScalarEllipticLshapeExact`, with the Neumann term `g(x) v dsn` of `exPoissonPkAdaptive`) the corner singularity reduces the energy rate on uniform meshes to $N^{-1/3}$; `exPoissonPkAdaptive` restores $N^{-1/2}$ by adaptive refinement. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
