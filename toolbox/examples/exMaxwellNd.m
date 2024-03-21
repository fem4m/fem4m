%[text] # Nédélec finite elements for a 2D Maxwell problem
%[text] Solves the two-dimensional model problem $\\mathrm{curl}\\,\\mathrm{rot}\\,u + \\alpha u = f$ of electromagnetics with the $H(\\mathrm{rot})$-conforming Nédélec elements of the first and second family. Unlike the Lagrange elements of `exPoissonP1`, these vector-valued elements are continuous only in the tangential components across the edges, and the tangential boundary condition is imposed on the edge degrees of freedom.
%[text] **Builds on:** `exPoissonP1`. **Next:** `exPoissonRTAdaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] In two dimensions fem4m uses two curl operators: the scalar rot of a vector field and the vector curl of a scalar function,
%[text]{"align":"center"} $ \\mathrm{rot}\\,v = \\partial\_1 v\_2 - \\partial\_2 v\_1, \\qquad \\mathrm{curl}\\,w = (\\partial\_2 w, -\\partial\_1 w), $
%[text] in the form strings `rot(v)` and `curl(w)`. Maxwell's equations for an electric field $(u\_1, u\_2, 0)$ that does not depend on $x\_3$ lead, after an implicit time discretization of the eddy-current problem, to the problem with a positive coefficient $\\alpha$: given $f$ and the tangential boundary data $u\_D\\cdot t$, find $u:\\Omega\\to\\mathbb{R}^2$ with
%[text]{"align":"center"} $ \\mathrm{curl}\\,\\mathrm{rot}\\,u + \\alpha u = f \\textrm{ in } \\Omega, \\qquad u\\cdot t = u\_D\\cdot t \\textrm{ on } \\partial\\Omega, $
%[text] where $t$ is the unit tangent of the boundary. The problem object `problemMaxwellExactSquare(alpha)` provides on the unit square $\\Omega = (0,1)^2$ the source $f = (\\pi^2+\\alpha)\\,u$ of the exact solution
%[text]{"align":"center"} $ u = \\left(\\sin(\\pi x\_2), \\sin(\\pi x\_1)\\right), \\qquad \\mathrm{rot}\\,u = \\pi\\cos(\\pi x\_1) - \\pi\\cos(\\pi x\_2), $
%[text] whose tangential component vanishes on $\\partial\\Omega$, so $u\_D\\cdot t = 0$.
%[text] ### Weak formulation
%[text] For a scalar $w$ and a vector field $v$ integration by parts gives, with the counterclockwise tangent $t$,
%[text]{"align":"center"} $ \\int\_\\Omega \\mathrm{curl}\\,w\\cdot v\\,dx = \\int\_\\Omega w\\,\\mathrm{rot}\\,v\\,dx - \\int\_{\\partial\\Omega} w\\,(v\\cdot t)\\,ds. $
%[text] The natural space is
%[text]{"align":"center"} $ H(\\mathrm{rot},\\Omega) = \\{v\\in L^2(\\Omega)^2 : \\mathrm{rot}\\,v\\in L^2(\\Omega)\\}, $
%[text] the two-dimensional version of $H(\\mathrm{curl},\\Omega)$, and $H\_0(\\mathrm{rot},\\Omega)$ its subspace with $v\\cdot t = 0$ on $\\partial\\Omega$. The weak formulation reads: find $u\\in H(\\mathrm{rot},\\Omega)$ with $u\\cdot t = u\_D\\cdot t$ on $\\partial\\Omega$ and
%[text]{"align":"center"} $ \\int\_\\Omega \\mathrm{rot}\\,u\\,\\mathrm{rot}\\,v\\,dx + \\alpha\\int\_\\Omega u\\cdot v\\,dx = \\int\_\\Omega f\\cdot v\\,dx \\quad \\textrm{for all } v\\in H\_0(\\mathrm{rot},\\Omega). $
%[text] For $\\alpha \> 0$ the bilinear form is coercive on $H(\\mathrm{rot},\\Omega)$. The rot part alone is not: every gradient $\\nabla\\varphi$ satisfies $\\mathrm{rot}\\,\\nabla\\varphi = 0$, so the term $\\alpha\\,u$ is essential for the well-posedness.
%[text] ### Discretization
%[text] A piecewise polynomial vector field belongs to $H(\\mathrm{rot},\\Omega)$ if and only if its tangential component $v\\cdot t\_E$ is continuous across every edge $E$; the normal component may jump. The Nédélec elements have exactly this continuity:
%[text] - the first family `Nd1k(k,mesh)` with $P\_{k-1}(T)^2\\subset V\_h|\_T\\subset P\_k(T)^2$, $k(k+2)$ functions per triangle, $k$ of them per edge, and $\\mathrm{rot}\\,V\_h|\_T = P\_{k-1}(T)$; for $k = 1$ the lowest-order element with one degree of freedom per edge,
%[text] - the second family `Nd2k(k,mesh)` with $V\_h|\_T = P\_k(T)^2$, $(k+1)(k+2)$ functions per triangle, $k+1$ of them per edge, and also $\\mathrm{rot}\\,V\_h|\_T = P\_{k-1}(T)$. \
%[text] Both use a hierarchical basis of edge functions, whose tangential traces live on one edge, and interior functions with vanishing tangential trace. The degrees of freedom of the boundary edges (`fe.fixedNodes`) are fixed by the $L^2$ projection of the tangential data $u\_D\\cdot t$ onto the edge functions, `fe.projectBoundary(mesh.ed4Db,g)` with the function `g(pts,tangents)` of the points and unit tangents; the other degrees of freedom (`fe.freeNodes`) are the unknowns.
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega \\mathrm{rot}\\,u\\,\\mathrm{rot}\\,v\\,dx + \\alpha\\int\_\\Omega u\\cdot v\\,dx$ | `rot(u)*rot(v) dx + alpha(x)*u*v dx` |
%[text] | $\\int\_\\Omega f\\cdot v\\,dx$ | `f(x)*v dx` |
%[text] | $\\Vert u-u\_h\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - u )^2 dx` |
%[text] | $\\Vert \\mathrm{rot}(u-u\_h)\\Vert\_{L^2(\\Omega)}$ | `( exact(x) - rot(u) )^2 dx` |
%[text:table]
%[text] ### Expected convergence
%[text] For the smooth solution and uniform refinement with $N\\sim h^{-2}$:
%[text:table]
%[text] | Element | $\\Vert\\mathrm{rot}(u-u\_h)\\Vert\_{L^2(\\Omega)}$ | $\\Vert u-u\_h\\Vert\_{L^2(\\Omega)}$ |
%[text] | --- | --- | --- |
%[text] | `Nd1k(k,mesh)` | $O(h^k) = O(N^{-k/2})$ | $O(h^k) = O(N^{-k/2})$ |
%[text] | `Nd2k(k,mesh)` | $O(h^k) = O(N^{-k/2})$ | $O(h^{k+1}) = O(N^{-(k+1)/2})$ |
%[text:table]
%[text] The rot of both families is a piecewise polynomial of degree $k-1$, so both converge with the same order in the rot; the full space $P\_k(T)^2$ of the second family gains one order in the $L^2$ error of $u$.
%%
%[text] ## Parameters
degree = 2;                      % polynomial degree k
space = @(mesh) Nd1k(degree,mesh);   % Nedelec elements of the first family
% space = @(mesh) Nd2k(degree,mesh);  % Nedelec elements of the second family, full P_k
alpha = 1;                       % coefficient of the zero-order term, alpha > 0
nrLevels = 5;                    % number of meshes
%%
%[text] ## Problem data and initial mesh
p = problemMaxwellExactSquare(alpha);
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
nrDoF4lvl = zeros(1,nrLevels);
L2Error4lvl = zeros(1,nrLevels);
rotError4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element space
    fe = space(mesh);
    nrDoF4lvl(lvl) = fe.nrDoF;
    % assemble
    a = BilinearForm('rot(u)*rot(v) dx + alpha(x)*u*v dx');
    A = a.assemble(fe,'alpha',alpha);
    l = LinearForm('f(x)*v dx');
    b = l.assemble(fe,'f',@p.f);
    % solve: tangential boundary values by projection onto the boundary edges
    x = fe.projectBoundary(mesh.ed4Db,@(pts,tangents)sum(p.u4Db(pts).*tangents,2));
    b = b - A*x;
    free = fe.freeNodes;
    x(free) = A(free,free)\b(free);
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(lvl) = l2.assemble(x,fe,12,'exact',@p.exactU);
    l2Rot = ErrorForm('( exact(x) - rot(u) )^2 dx');
    rotError4lvl(lvl) = l2Rot.assemble(x,fe,12,'exact',@p.exactRotU);
    % refine
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Convergence table
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
rotRate = [NaN, -diff(log(rotError4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',rotError4lvl',rotRate', ...
    VariableNames=["N" "L2 error" "L2 rate" "rot error" "rot rate"])
%[text] For the first family with $k = 2$ both rates approach $k/2 = 1$.
%%
%[text] ### Discrete solution
%[text] The magnitude $|u\_h|$ of the discrete field on the finest mesh.
clf
plot(fe,x);
view([40,35]);
colorbar;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,rotError4lvl,'||rot(u-u_h)||_{L^2}');
hold off
%[text] The slopes in the legend agree with the expected rates of the chosen family.
%%
%[text] ## Possibilities for extensions
%[text] - Use the second family, the commented line in the parameters: the rot error is the same as for the first family, while the $L^2$ error of $u$ converges like $N^{-3/2}$ for $k = 2$.
%[text] - Set `degree = 1` for the lowest-order Nédélec element of the first family: both errors converge like $N^{-1/2}$.
%[text] - Set `alpha = 0`: the gradients of the continuous $P\_k$ functions (for `Nd2k` of the $P\_{k+1}$ functions) that vanish on $\\partial\\Omega$ lie in the discrete space and have zero rot, so the matrix is singular and MATLAB warns about it.
%[text] - Use `problemMaxwellExactLshape(alpha)` on the L-shaped domain, whose field behaves like $r^{-1/3}$ at the reentrant corner: the $L^2$ error converges only like $N^{-1/3}$, while $\\mathrm{rot}\\,u$ stays bounded (it behaves like $r^{2/3}$) and its error converges almost with the optimal rate.
%[text] - In two dimensions the first and second Nédélec families are the Raviart-Thomas elements of `exPoissonRTAdaptive` and the Brezzi-Douglas-Marini elements, rotated by 90 degrees: tangential continuity becomes normal continuity, $H(\\mathrm{rot})$ becomes $H(\\mathrm{div})$. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
