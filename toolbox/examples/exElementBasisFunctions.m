%[text] # Basis functions on the reference triangle
%[text] Plots the nodal basis functions of the Lagrange element $P\_k$ on the reference triangle and explains how fem4m builds a finite element space from them: barycentric coordinates, Lagrange points, the affine element map and the nodal interpolation, whose errors decay with the optimal rates on uniformly refined meshes.
%[text] **Builds on:** `exMeshEnumeration`. **Next:** `exPoissonP1`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### Reference triangle and barycentric coordinates
%[text] The reference triangle $\\hat T$ has the vertices $\\hat a\_1 = (0,0)$, $\\hat a\_2 = (1,0)$ and $\\hat a\_3 = (0,1)$. Its barycentric coordinates
%[text]{"align":"center"} $ \\lambda\_1 = 1 - \\hat x\_1 - \\hat x\_2, \\qquad \\lambda\_2 = \\hat x\_1, \\qquad \\lambda\_3 = \\hat x\_2 $
%[text] are affine functions with $\\lambda\_j(\\hat a\_i) = \\delta\_{ij}$ and $\\lambda\_1 + \\lambda\_2 + \\lambda\_3 = 1$; they are nonnegative exactly on $\\hat T$, and every point is the weighted mean $\\hat x = \\sum\_j \\lambda\_j(\\hat x)\\,\\hat a\_j$ of the vertices. The edge opposite to $\\hat a\_j$ is the set where $\\lambda\_j = 0$.
%[text] ### Polynomials of degree k
%[text] The space $P\_k(\\hat T)$ of the polynomials of total degree at most $k$ is spanned by the monomials $\\hat x\_1^i\\hat x\_2^j$ with $i + j\\le k$. There are $m+1$ monomials of degree exactly $m$, so
%[text]{"align":"center"} $ \\dim P\_k(\\hat T) = \\sum\_{m=0}^{k} (m+1) = \\frac{(k+1)(k+2)}{2}. $
%[text] ### Lagrange points and the nodal basis
%[text] The Lagrange points of degree $k$ are the points of $\\hat T$ whose barycentric coordinates are multiples of $1/k$: the 3 vertices, $k-1$ points on every edge and $(k-1)(k-2)/2$ interior points, together $n = (k+1)(k+2)/2$ points $\\hat z\_1,\\dots,\\hat z\_n$. They are unisolvent: a polynomial $q\\in P\_k(\\hat T)$ that vanishes at all of them restricts on every edge to a polynomial of degree $k$ in one variable with $k+1$ zeros, so $q = \\lambda\_1\\lambda\_2\\lambda\_3\\,r$ with $r\\in P\_{k-3}(\\hat T)$, and $r$ vanishes at the interior points, the Lagrange points of degree $k-3$ of a smaller triangle; by induction $q = 0$. Hence there is a unique nodal basis $\\hat\\varphi\_1,\\dots,\\hat\\varphi\_n$ of $P\_k(\\hat T)$ with
%[text]{"align":"center"} $ \\hat\\varphi\_j(\\hat z\_i) = \\delta\_{ij} \\quad \\textrm{for } i,j = 1,\\dots,n. $
%[text] `Pk` writes the basis functions in the monomials, $\\hat\\varphi\_j = \\sum\_m C\_{jm}\\,\\hat x^{\\alpha\_m}$ with the exponents $\\alpha\_m$, and computes the coefficient matrix $C$ as the inverse of the Vandermonde matrix $V\_{mi} = \\hat z\_i^{\\alpha\_m}$, which is exactly the condition above. The Lagrange points `fe.c4refDof` are ordered as follows: the vertices $\\hat a\_1, \\hat a\_2, \\hat a\_3$; then, for $j = 1,\\dots,k-1$, the points $(j/k,0)$, $(1-j/k,j/k)$ and $(0,1-j/k)$ on the edges from $\\hat a\_1$ to $\\hat a\_2$, from $\\hat a\_2$ to $\\hat a\_3$ and from $\\hat a\_3$ to $\\hat a\_1$; finally the interior points. For $k = 1$ the basis functions are the barycentric coordinates; for $k = 2$ they are $\\lambda\_j(2\\lambda\_j - 1)$ at the vertex $\\hat a\_j$ and $4\\lambda\_i\\lambda\_j$ at the midpoint of the edge from $\\hat a\_i$ to $\\hat a\_j$.
%[text] ### Element map and global basis
%[text] A triangle $T$ of a mesh with the nodes $P\_1, P\_2, P\_3$, in the order of `n4e`, is the image of $\\hat T$ under the affine map
%[text]{"align":"center"} $ F\_T(\\hat x) = P\_1 + (P\_2 - P\_1)\\,\\hat x\_1 + (P\_3 - P\_1)\\,\\hat x\_2 $
%[text] (`mesh.point4e`) with the constant Jacobian $B\_T = (P\_2 - P\_1, P\_3 - P\_1)$ and $\\det B\_T = 2|T|$. The basis functions on $T$ are $\\varphi\_j = \\hat\\varphi\_j\\circ F\_T^{-1}$ with the Lagrange points $F\_T(\\hat z\_j)$ (`fe.c4dof`), and the chain rule gives their gradients $\\nabla\\varphi\_j = B\_T^{-T}\\,(\\hat\\nabla\\hat\\varphi\_j)\\circ F\_T^{-1}$ with the matrix $B\_T^{-T}$ of `mesh.transformation`. Local functions of neighbouring elements that belong to the same Lagrange point form one global basis function; `fe.dof4e(T,j)` is the global number of the local function $j$ of $T$. The global functions are continuous, since on a common edge the traces from both sides are polynomials of degree $k$ in one variable that agree at the $k+1$ Lagrange points of the edge. The space $V\_h = P\_k(\\mathcal{T})\\cap C(\\overline\\Omega)$ (`Pk(k,mesh)`) therefore has the dimension
%[text]{"align":"center"} $ N = |\\mathcal{N}| + (k-1)\\,|\\mathcal{E}| + \\frac{(k-1)(k-2)}{2}\\,|\\mathcal{T}|. $
%[text] ### Nodal interpolation
%[text] The nodal interpolation $I\_h v = \\sum\_{j=1}^N v(z\_j)\\,\\varphi\_j$ of a continuous function $v$ (`fe.interpolate(v)`) takes the values at the Lagrange points $z\_j$ as coefficients. It reproduces the polynomials of degree $k$ on every element, and the Bramble-Hilbert lemma together with the transformation to $\\hat T$ gives for $v\\in H^{k+1}(\\Omega)$ and the maximal mesh size $h$
%[text]{"align":"center"} $ \\|v - I\_h v\\|\_{L^2(\\Omega)} \\le C\\,h^{k+1}\\,|v|\_{H^{k+1}(\\Omega)}, \\qquad \\|\\nabla(v - I\_h v)\\|\_{L^2(\\Omega)} \\le C\\,h^{k}\\,|v|\_{H^{k+1}(\\Omega)}, $
%[text] where $C$ depends only on $k$ and on the shape regularity of the mesh. Uniform refinement halves $h$ and multiplies $N$ by about four, $N\\approx h^{-2}$, so the two errors decay like $N^{-(k+1)/2}$ and $N^{-k/2}$. These are the rates of the finite element method in `exPoissonP1`: by Céa's lemma the Galerkin solution is at least as accurate as $I\_h u$ in the energy norm.
%%
%[text] ## Parameters
degree = 2;    % polynomial degree k
nrLevels = 6;  % number of meshes for the interpolation errors
%%
%[text] ## Reference triangle
%[text] The mesh consists of the reference triangle with the nodes in the order $\\hat a\_1, \\hat a\_2, \\hat a\_3$, so the element map is the identity and the global basis functions are the reference ones; `Mesh(c4n,n4e)` makes all boundary edges Dirichlet edges. The table lists the reference Lagrange points `fe.c4refDof` in the order of the basis functions, the plot their images `fe.c4dof` on the element.
mesh = Mesh([0,0;1,0;0,1],[1 2 3]);
fe = Pk(degree,mesh);
table((1:fe.nrDoF)',fe.c4refDof,VariableNames=["j" "Lagrange point"])
plot(mesh);
hold on
plot(fe.c4dof(:,1),fe.c4dof(:,2),'o','MarkerFaceColor','auto');
text(fe.c4dof(:,1)+0.02,fe.c4dof(:,2)+0.03,string(1:fe.nrDoF),'FontSize',14);
hold off
axis equal
%%
%[text] ## Basis functions
%[text] The coefficient vector with the entry 1 at position $j$ and 0 elsewhere is the basis function $\\varphi\_j$. `fe.plot(c,4)` evaluates the function on the reference triangle refined four times.
clf
tiledlayout('flow');
for j = 1:fe.nrDoF
    c = zeros(fe.nrDoF,1);
    c(j) = 1;
    nexttile
    fe.plot(c,4);
    shading interp;
    view(3);
    title("$\varphi_{" + j + "}$",'Interpreter','latex');
end
%[text] Every basis function is 1 at its own Lagrange point and vanishes at all other Lagrange points. For $k\\ge 2$ some of them change sign, e.g. the vertex functions $\\lambda\_j(2\\lambda\_j - 1)$ of $P\_2$, which are negative where $0 \< \\lambda\_j \< 1/2$.
%%
%[text] ## Nodal basis
%[text] The matrix of the values $\\hat\\varphi\_j(\\hat z\_i)$ is the identity up to rounding errors, and the basis functions sum up to one, the interpolant of the constant function 1.
values = zeros(fe.nrDoF);
for i = 1:fe.nrDoF
    values(i,:) = fe.getBasis([],fe.c4refDof(i,:),[]);
end
deviationFromIdentity = norm(values - eye(fe.nrDoF),inf)
sumOfBasisFunctions = sum(fe.getBasis([],[0.2,0.3],[]))
%%
%[text] ## Nodal interpolation on uniform meshes
%[text] The loop interpolates the smooth function $v(x) = \\sin(2\\pi x\_1)\\cos(\\pi x\_2)$ on the reference triangle and its uniform refinements and computes the errors in the $L^2$ norm and in the energy norm.
v = @(x) sin(2*pi*x(:,1)).*cos(pi*x(:,2));
gradV = @(x) [2*pi*cos(2*pi*x(:,1)).*cos(pi*x(:,2)), -pi*sin(2*pi*x(:,1)).*sin(pi*x(:,2))];
nrDoF4lvl = zeros(1,nrLevels);
L2Error4lvl = zeros(1,nrLevels);
H1Error4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element space
    fe = Pk(degree,mesh);
    nrDoF4lvl(lvl) = fe.nrDoF;
    % interpolate
    x = fe.interpolate(v);
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(lvl) = l2.assemble(x,fe,'exact',v);
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(lvl) = h1.assemble(x,fe,'exact',gradV);
    % refine
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Interpolation errors
%[text] The rates are the slopes of the errors over $N$ in a log-log scale between two consecutive meshes.
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
H1Rate = [NaN, -diff(log(H1Error4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',H1Error4lvl',H1Rate', ...
    VariableNames=["N" "L2 error" "L2 rate" "H1 error" "H1 rate"])
%[text] The rates approach $(k+1)/2$ for the $L^2$ error and $k/2$ for the energy error, as predicted.
%%
%[text] ### Interpolant
clf
plot(fe,x);
view([40,35]);
colorbar;
shading interp;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||v-I_hv||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(v-I_hv)||_{L^2}');
hold off
%[text] The legend shows the least-squares slopes over all meshes. They include the coarsest meshes with only a few elements and are therefore somewhat larger than the asymptotic rates of the table.
%%
%[text] ## Possibilities for extensions
%[text] - Set `degree = 1`, `3` or `4`: there are $(k+1)(k+2)/2$ basis functions, and the interpolation errors decay like $N^{-(k+1)/2}$ and $N^{-k/2}$.
%[text] - Map the basis to a physical triangle, e.g. `mesh = Mesh([0,0;2,0;0.5,1],[1 2 3]);`: the plots show $\\hat\\varphi\_j\\circ F\_T^{-1}$ at the Lagrange points `fe.c4dof` $= F\_T(\\hat z\_j)$, while the table and the check of the nodal basis stay on the reference triangle.
%[text] - Plot, with the loop of the section Basis functions only, the hierarchical basis of `HPk(degree,mesh)`, which consists of the $P\_1$ functions, integrated Legendre polynomials on the edges and bubble functions, or the nonconforming Crouzeix-Raviart basis of `CRk(1,mesh)`, whose degrees of freedom are the mean values on the edges. These elements have no Lagrange points `c4refDof` and no nodal interpolation.
%[text] - `Pk3d` is the Lagrange element on tetrahedra with $(k+1)(k+2)(k+3)/6$ basis functions per element; `exPoisson3dPkAdaptive` uses it. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
