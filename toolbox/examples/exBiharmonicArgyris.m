%[text] # $C^1$ Argyris finite elements for the biharmonic problem
%[text] Solves the biharmonic equation of a clamped plate with the quintic Argyris element, whose piecewise polynomials are continuously differentiable across the edges and therefore conforming in $H^2$. The example checks the $C^1$ continuity of the discrete solution and the high convergence rates of the element; the TUBA elements of higher degree are an alternative in the parameters.
%[text] **Builds on:** `exPoissonP1`. **Next:** `exBiharmonicHHJ`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] The vertical deflection $u$ of a thin plate that is clamped along its boundary and loaded by the force density $f$ solves the biharmonic equation with clamped boundary conditions
%[text]{"align":"center"} $ \\Delta^2 u = f \\textrm{ in } \\Omega, \\qquad u = 0 \\textrm{ and } \\nabla u\\cdot n = 0 \\textrm{ on } \\partial\\Omega. $
%[text] The problem object `problemBiharmonicSquareClampedExact` provides, on the unit square $\\Omega = (0,1)^2$, the load $f = \\Delta^2 u$ of the polynomial exact solution
%[text]{"align":"center"} $ u = x\_1^2(1-x\_1)^2\\,x\_2^2(1-x\_2)^2 $
%[text] of degree 8, whose value and normal derivative vanish on the boundary.
%[text] ### Weak formulation
%[text] Since $\\Delta^2 u = \\mathrm{div}\\,\\mathrm{div}\\,D^2 u$ with the Hessian $D^2 u$, two integrations by parts against a test function $v$ with $v = \\nabla v\\cdot n = 0$ on $\\partial\\Omega$ give: find $u\\in H^2\_0(\\Omega)$ with
%[text]{"align":"center"} $ \\int\_\\Omega D^2 u : D^2 v\\,dx = \\int\_\\Omega f v\\,dx \\quad \\textrm{for all } v\\in H^2\_0(\\Omega). $
%[text] Both boundary conditions are essential and belong to the space $H^2\_0(\\Omega)$. On this space $\\|D^2 v\\|\_{L^2(\\Omega)}$ is a norm, the bilinear form is coercive, and the Lax-Milgram lemma gives a unique solution.
%[text] ### Why $C^1$ elements
%[text] A piecewise polynomial belongs to $H^2(\\Omega)$ if and only if it is continuous and has continuous first derivatives across all edges, i.e. it is $C^1$. The Lagrange elements of `exPoissonP1` are only continuous; their gradients jump across the edges, so they are not conforming for the biharmonic problem. On a general triangulation, $C^1$ piecewise polynomials with local degrees of freedom require at least degree 5.
%[text] ### Discretization
%[text] The Argyris element (`Argyris(5,mesh)`) gives the space $V\_h$ of the functions $v\\in C^1(\\overline\\Omega)$ with $v|\_T\\in P\_5(T)$ for all $T\\in\\mathcal{T}$ that are, in addition, twice continuously differentiable at the vertices, a subspace of $P\_5(\\mathcal{T})\\cap C^1(\\overline\\Omega)$. Its 21 local degrees of freedom, the dimension of $P\_5(T)$, are
%[text] - the value, the two first and the three second derivatives $\\partial\_{11}u$, $\\partial\_{22}u$, $\\partial\_{12}u$ at each vertex (6 per vertex),
%[text] - the normal derivative at each edge midpoint (1 per edge). \
%[text] The values and all derivatives at a vertex are shared by its triangles, which makes $u\_h$ and $\\nabla u\_h$ continuous at the vertices and, together with the normal derivatives at the midpoints, $C^1$ across the edges. Globally there are $N = 6\\,|\\mathcal{N}| + |\\mathcal{E}|$ degrees of freedom for the nodes $\\mathcal{N}$ and edges $\\mathcal{E}$. The clamped conditions fix (`fe.fixedNodes`), at the boundary vertices, the value, the gradient, the mixed derivative $\\partial\_{12}u$ and the second tangential derivative, and the normal derivatives at the midpoints of the boundary edges; this assumes boundary edges parallel to the axes, as on the square. The remaining degrees of freedom (`fe.freeNodes`) are the unknowns.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega D^2 u : D^2 v\\,dx$ | `D2(u):D2(v) dx` |
%[text] | $\\int\_\\Omega f v\\,dx$ | `f(x) v dx` |
%[text] | $\\Vert u-u\_h\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - u )^2 dx` |
%[text] | $\\Vert \\nabla(u-u\_h)\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - grad(u) )^2 dx` |
%[text] | $\\Vert D^2(u-u\_h)\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - D2(u) )^2 dx` |
%[text:table]
%[text] ### Checking the $C^1$ conformity
%[text] The jumps of $u\_h$ and of its normal derivative across the interior edges $\\mathcal{E}(\\Omega)$, scaled like the terms of a residual-based error estimator,
%[text]{"align":"center"} $ \\Big(\\sum\_{E\\in\\mathcal{E}(\\Omega)} h\_E^{-1}\\|\[u\_h\]\\|^2\_{L^2(E)} + h\_E\\|\[\\nabla u\_h\\cdot n\_E\]\\|^2\_{L^2(E)}\\Big)^{1/2}, $
%[text] assembled by `EstimatorForm('h^-1*([u])^2 dsi + h^1*([grad(u)*n])^2 dsi')`, vanish for every $C^1$ function. For the discrete solutions of the Argyris element they are zero up to rounding errors, which checks the construction of the basis.
%[text] ### Expected convergence
%[text] Céa's lemma and the interpolation error estimates give $\\|D^2(u-u\_h)\\|\_{L^2(\\Omega)} \\le C h^4$ for the smooth solution; duality arguments on the square add one order for each lower derivative. With $N\\sim h^{-2}$, uniform refinement gives
%[text]{"align":"center"} $ \\|D^2(u-u\_h)\\|\_{L^2(\\Omega)} = O(N^{-2}), \\qquad \\|\\nabla(u-u\_h)\\|\_{L^2(\\Omega)} = O(N^{-5/2}), \\qquad \\|u-u\_h\\|\_{L^2(\\Omega)} = O(N^{-3}). $
%[text] For $C^1$ elements of degree $k\\ge 5$ the rates are $N^{-(k-1)/2}$, $N^{-k/2}$ and $N^{-(k+1)/2}$.
%%
%[text] ## Parameters
space = @(mesh) Argyris(5,mesh);   % quintic Argyris element
% space = @(mesh) TUBAk(6,mesh);   % C^1 element of degree 6 of the TUBA family TUBAk(k), k >= 5
nrLevels = 4;                      % number of meshes
%%
%[text] ## Problem data and initial mesh
p = problemBiharmonicSquareClampedExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
nrDoF4lvl = zeros(1,nrLevels);
L2Error4lvl = zeros(1,nrLevels);
H1Error4lvl = zeros(1,nrLevels);
H2Error4lvl = zeros(1,nrLevels);
jump4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element space
    fe = space(mesh);
    nrDoF4lvl(lvl) = fe.nrDoF;
    % assemble
    a = BilinearForm('D2(u):D2(v) dx');
    A = a.assemble(fe);
    l = LinearForm('f(x) v dx');
    b = l.assemble(fe,'f',@p.f);
    % solve
    x = zeros(fe.nrDoF,1);
    free = fe.freeNodes;
    x(free) = A(free,free)\b(free);
    % check the C1 conformity by the jumps of u and of its normal derivative
    jump = EstimatorForm('h^-1*([u])^2 dsi + h^1*([grad(u)*n])^2 dsi');
    jump4lvl(lvl) = norm(jump.assemble(x,fe));
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(lvl) = l2.assemble(x,fe,12,'exact',@p.exactU);
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(lvl) = h1.assemble(x,fe,12,'exact',@p.exactGradU);
    h2 = ErrorForm('( exact(x) - D2(u) )^2 dx');
    H2Error4lvl(lvl) = h2.assemble(x,fe,12,'exact',@p.exactD2);
    % refine
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Convergence table
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
H1Rate = [NaN, -diff(log(H1Error4lvl))./diff(log(nrDoF4lvl))];
H2Rate = [NaN, -diff(log(H2Error4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',jump4lvl',L2Error4lvl',L2Rate',H1Error4lvl',H1Rate',H2Error4lvl',H2Rate', ...
    VariableNames=["N" "C1 jumps" "L2 error" "L2 rate" "H1 error" "H1 rate" "H2 error" "H2 rate"])
%[text] The jumps are of the size of rounding errors on all meshes: the discrete solutions are $C^1$. The rates lie above 3, 5/2 and 2 on these coarse meshes, where one refinement multiplies $N$ only by about 3 to 3.6 instead of the asymptotic factor 4, and decrease towards them.
%%
%[text] ### Discrete solution
%[text] The third argument of `plot` refines every triangle twice for the plot, which resolves the quintic polynomials better.
clf
plot(fe,x,2);
view([40,35]);
colorbar;
shading interp;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(u-u_h)||_{L^2}');
plotConvergence(nrDoF4lvl,H2Error4lvl,'||D^2(u-u_h)||_{L^2}');
hold off
%[text] Each derivative less gains half an order in $N$, one order in $h$.
%%
%[text] ## Possibilities for extensions
%[text] - Use the TUBA element of degree 6, the commented line in the parameters. `TUBAk(k,mesh)` extends the Argyris element to all degrees $k\\ge 5$ by further degrees of freedom on the edges and in the interior, `TUBAk(5,mesh)` spans the Argyris space. For $k = 6$ the rates grow to $N^{-7/2}$, $N^{-3}$ and $N^{-5/2}$.
%[text] - Use `problemBiharmonicSquareClampedOscExact` with the oscillating solution $\\sin^2(\\pi x\_1)\\sin^2(\\pi x\_2)$: the rates are the same, the errors larger.
%[text] - Replace the bilinear form by `'laplace(u)*laplace(v) dx'`, i.e. $\\int\_\\Omega\\Delta u\\,\\Delta v\\,dx$: for functions in $H^2\_0(\\Omega)$ both forms coincide, and the discrete solutions agree up to rounding errors.
%[text] - `exBiharmonicHHJ` avoids $C^1$ elements altogether: its mixed formulation needs only continuous Lagrange elements for $u$ and normal-normal continuous matrix fields for the Hessian. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
