%[text] # Adaptive $P\_k$ finite elements for an anisotropic diffusion eigenvalue problem
%[text] Computes a multiple eigenvalue of a diffusion operator with a smooth, full and anisotropic diffusion matrix on the unit square, whose eigenpairs are known analytically, with conforming $P\_k$ finite elements on adaptively refined meshes. New are a matrix-valued coefficient and its divergence in the forms, an error estimator for a cluster of eigenvalues and the errors of the eigenfunctions with respect to the exact eigenspace.
%[text] **Builds on:** `exEigenvaluePkAdaptive`. **Next:** `exElasticityPkAdaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] On the unit square $\\Omega = (0,1)^2$ find the eigenpairs $(\\lambda,u)$ with $u\\neq 0$ and
%[text]{"align":"center"} $ -\\mathrm{div}(\\alpha\\nabla u) = \\lambda u \\textrm{ in } \\Omega, \\qquad u = 0 \\textrm{ on } \\partial\\Omega, $
%[text] for the symmetric positive definite diffusion matrix $\\alpha$ of the problem object `problemScalarEllipticEigenvalueSquareTwist`.
%[text] ### Diffusion coefficient
%[text] The smooth, full and anisotropic diffusion matrix $\\alpha$ stems from an area preserving twist $\\Phi$ of the unit square, which rotates the circles around the center $c = (1/2,1/2)$ by an angle $g(r)$ that depends on their radius $r = |x-c|$:
%[text]{"align":"center"} $ \\Phi(x) = c + Q(g(r))\\,(x-c), \\qquad Q(\\theta) = \\left(\\begin{array}{cc} \\cos\\theta & -\\sin\\theta \\\\ \\sin\\theta & \\cos\\theta \\end{array}\\right), \\qquad g(r) = \\left\\{\\begin{array}{ll} \\theta\_0\\,(1-r^2/\\rho^2)^4 & \\textrm{for } r \\leq \\rho, \\\\ 0 & \\textrm{otherwise,} \\end{array}\\right. $
%[text] with the maximal twist angle $\\theta\_0$ (`theta0`, default $\\pi/2$) and the radius $\\rho = 1/2$. The twist is the identity outside the disk $B\_\\rho(c)$, hence $\\Phi(\\Omega) = \\Omega$ and $\\Phi = \\mathrm{id}$ on $\\partial\\Omega$. With $(\\xi,\\eta) = x-c$, $J = \\left(\\begin{array}{cc} 0 & -1 \\\\ 1 & 0 \\end{array}\\right)$ and
%[text]{"align":"center"} $ \\kappa(r) = \\frac{g'(r)}{r} = -\\frac{8\\theta\_0}{\\rho^2}\\Big(1-\\frac{r^2}{\\rho^2}\\Big)^3 \\quad \\textrm{for } r \\leq \\rho, \\qquad \\kappa(r) = 0 \\quad \\textrm{otherwise}, $
%[text] the Jacobian is $D\\Phi = Q(g)\\,(I + \\kappa\\, J (x-c)(x-c)^T)$ with $\\det D\\Phi = 1$, since $(J(x-c)(x-c)^T)^2 = 0$. The diffusion matrix
%[text]{"align":"center"} $ \\alpha = D\\Phi^{-1} D\\Phi^{-T} = \\left(\\begin{array}{cc} 1 + 2\\kappa\\xi\\eta + \\kappa^2 r^2 \\eta^2 & -\\kappa(\\xi^2-\\eta^2) - \\kappa^2 r^2 \\xi\\eta \\\\ -\\kappa(\\xi^2-\\eta^2) - \\kappa^2 r^2 \\xi\\eta & 1 - 2\\kappa\\xi\\eta + \\kappa^2 r^2 \\xi^2 \\end{array}\\right) $
%[text] is a polynomial in $\\xi,\\eta$ inside the disk, equals the identity outside, belongs to $C^2(\\overline{\\Omega})$ and is symmetric positive definite with $\\det\\alpha = 1$. In the basis of the radial and the tangential direction $\\alpha = \\left(\\begin{array}{cc} 1 & -s \\\\ -s & 1+s^2 \\end{array}\\right)$ with $s = r\\,g'(r)$, so the eigenvalues of $\\alpha$ are $\\big((\\sqrt{s^2+4} \\pm s)/2\\big)^2$ and the condition number of $\\alpha$ is about 12 for $\\theta\_0 = \\pi/2$.
%[text] ### Exact eigenpairs
%[text] For $u = \\tilde u\\circ\\Phi$ the transformation rule with $\\det D\\Phi = 1$ gives
%[text]{"align":"center"} $ \\int\_\\Omega \\alpha\\nabla u\\cdot\\nabla v\\,dx = \\int\_\\Omega \\nabla\\tilde u\\cdot\\nabla\\tilde v\\,dy, \\qquad \\int\_\\Omega u v\\,dx = \\int\_\\Omega \\tilde u\\tilde v\\,dy, $
%[text] so the eigenvalues $\\lambda = \\pi^2(j^2+k^2)$ are those of the Laplacian on the unit square and the $L^2$-normalized eigenfunctions are $u = 2\\sin(j\\pi\\Phi\_1)\\sin(k\\pi\\Phi\_2)$ with the gradients $\\nabla u = D\\Phi^T\\,\\nabla\\tilde u(\\Phi)$. The multiplicity of $\\lambda = \\pi^2 n$ is the number of pairs $j,k\\geq 1$ with $j^2+k^2 = n$, e.g. two for $5\\pi^2$. `p.eigenvalue(m)` returns the $m$-th eigenvalue, counted with multiplicity, `p.eigenspace(m)` the pairs $\[j\\ k\]$ of its eigenspace, `p.cluster(m)` the numbers of all modes with the same eigenvalue, and `p.exactU(x,j,k)` and `p.exactGradU(x,j,k)` the eigenfunctions and their gradients.
%[text] ### Weak formulation and discretization
%[text] Find $\\lambda\\in\\mathbb{R}$ and $u\\in H^1\_0(\\Omega)$ with $\\|u\\|\_{L^2(\\Omega)} = 1$ and
%[text]{"align":"center"} $ \\int\_\\Omega \\alpha\\nabla u\\cdot\\nabla v\\,dx = \\lambda\\int\_\\Omega u v\\,dx \\quad \\textrm{for all } v\\in H^1\_0(\\Omega). $
%[text] As in `exEigenvaluePkAdaptive`, the space $V\_h = P\_k(\\mathcal{T})\\cap H^1\_0(\\Omega)$ (`Pk(degree,mesh)`) leads to the generalized matrix eigenvalue problem $Ax = \\lambda\_h Bx$ with the stiffness matrix $A$ of the diffusion matrix and the mass matrix $B$. A matrix-valued coefficient, here `alpha(x)` with the function `@p.alpha`, multiplies the gradient in front of which it stands.
%[text:table]
%[text] | Matrix | Integral | fem4m form |
%[text] | --- | --- | --- |
%[text] | $A$ | $\\int\_\\Omega \\alpha\\nabla u\\cdot\\nabla v\\,dx$ | `alpha(x)grad(u)*grad(v) dx` |
%[text] | $B$ | $\\int\_\\Omega u v\\,dx$ | `u * v dx` |
%[text:table]
%[text] ### Eigenvalue clusters and error estimator
%[text] On the discrete level the double eigenvalue $\\lambda = 5\\pi^2$ of the default mode $m = 2$ splits into two close discrete eigenvalues $\\lambda\_{h,2}\\le\\lambda\_{h,3}$, and a single discrete eigenfunction need not converge to a fixed eigenfunction. The estimator therefore treats the whole cluster $\\mathcal{C}$ = `p.cluster(m)` of the modes of the eigenvalue at once: it sums the squared residual indicators of all discrete eigenpairs $(\\lambda\_{h,i},u\_{h,i})$, $i\\in\\mathcal{C}$, with the right-hand sides $\\lambda\_{h,i}u\_{h,i}$,
%[text]{"align":"center"} $ \\eta^2 = \\sum\_{i\\in\\mathcal{C}} \\Big( \\sum\_{T\\in\\mathcal{T}} h\_T^2 \\|\\lambda\_{h,i} u\_{h,i} + \\mathrm{div}(\\alpha\\nabla u\_{h,i})\\|^2\_{L^2(T)} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} h\_E \\|\[\\alpha\\nabla u\_{h,i}\\cdot n\_E\]\\|^2\_{L^2(E)} \\Big). $
%[text] The volume term needs the divergence of the discrete flux, $\\mathrm{div}(\\alpha\\nabla u\_h) = (\\mathrm{div}\\,\\alpha)\\cdot\\nabla u\_h + \\alpha : D^2 u\_h$ with the row-wise divergence of $\\alpha$, and thus the derivatives of $\\alpha$: the coefficient `alphac(x)` inside `div` receives the cell array `{@p.alpha,@p.dalpha}` of the matrix and its derivative.
%[text:table]
%[text] | Term | fem4m form |
%[text] | --- | --- |
%[text] | $h\_T^2\\Vert \\lambda\_h u\_h+\\mathrm{div}(\\alpha\\nabla u\_h)\\Vert^2\_{L^2(T)}$ | `h^2*(f(x)+div(alphac(x)grad(u)))^2 dx` |
%[text] | $h\_E\\Vert\[\\alpha\\nabla u\_h\\cdot n\_E\]\\Vert^2\_{L^2(E)}$ | `h^1*([alpha(x)*grad(u)*n])^2 dsi` |
%[text:table]
%[text] As for a simple eigenvalue, the eigenvalue errors behave like $\\eta^2$. The adaptive algorithm is that of `exEigenvaluePkAdaptive`, with bulk marking of the cluster indicators $\\eta\_T$ (`MarkBulk` with $\\theta = 1/2$) and red-green-blue refinement.
%[text] ### Eigenfunction errors
%[text] The eigenfunction errors are the distances of the discrete eigenfunction $u\_h$ to the exact eigenspace $E(\\lambda)$, measured against the $L^2$-projection $u^\* = P\_E u\_h = \\sum\_l (u\_h,u\_l)\\,u\_l$ with the $L^2$-orthonormal eigenfunctions $u\_l$ of $E(\\lambda)$. This fixes the sign for simple eigenvalues and the unknown limit in $E(\\lambda)$ for multiple eigenvalues. The coefficients $(u\_h,u\_l) = \\sum\_n x\_n\\int\_\\Omega u\_l\\varphi\_n\\,dx$ follow from the load vector of `f(x) v dx` with $f = u\_l$.
%[text] ### Expected convergence
%[text] The diffusion matrix and the eigenfunctions are sufficiently smooth, $\\alpha\\in C^2(\\overline\\Omega)$ and $u\\in W^{4,\\infty}(\\Omega)$ since $g\\in W^{4,\\infty}$, and the square is convex, so for $k\\le 3$ the rates are optimal already on uniform meshes: $\\|\\nabla(u^\*-u\_h)\\|\_{L^2(\\Omega)} = O(N^{-k/2})$, $\\|u^\*-u\_h\\|\_{L^2(\\Omega)} = O(N^{-(k+1)/2})$ and $\\lambda\_h - \\lambda = O(N^{-k})$. The adaptive meshes adapt to the variation of $\\alpha$ in the twisted disk.
%%
%[text] ## Parameters
degree = 1;       % polynomial degree k
m = 2;            % compute the m-th eigenvalue lambda_m, here the double eigenvalue 5*pi^2
theta0 = pi/2;    % maximal twist angle
minNrDoF = 1e4;   % stop as soon as the number of degrees of freedom exceeds minNrDoF
mark = MarkBulk;  % bulk (Doerfler) marking with theta = 1/2
% mark = MarkUniform;  % uniform refinement, for comparison
%%
%[text] ## Problem data and initial mesh
p = problemScalarEllipticEigenvalueSquareTwist(theta0);
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform',2);
plot(mesh);
pairs4l = p.eigenspace(m);  % index pairs [j k] of the exact eigenspace
modes = p.cluster(m);       % discrete modes approximating the same eigenvalue
%%
%[text] ## Adaptive loop
%[text] Every iteration runs SOLVE, ESTIMATE, MARK and REFINE; the loop stops when the number of degrees of freedom exceeds `minNrDoF`.
nrDoF4lvl = [];
lambda4lvl = [];
lambdaError4lvl = [];
L2Error4lvl = [];
H1Error4lvl = [];
eta4lvl = [];
while true
    % finite element space
    fe = Pk(degree,mesh);
    nrDoF4lvl(end+1) = fe.nrDoF;
    % assemble
    a = BilinearForm('alpha(x)grad(u)*grad(v) dx');
    A = a.assemble(fe,'alpha',@p.alpha);
    b = BilinearForm('u * v dx');
    B = b.assemble(fe);
    % solve
    free = fe.freeNodes;
    [V,D] = eigs(A(free,free),B(free,free),max(modes),'smallestabs');
    [lambda4mode,ind] = sort(diag(D));
    V = V(:,ind);
    V = V./sqrt(sum(V.*(B(free,free)*V),1));  % ||u_h||_L2 = 1
    X = zeros(fe.nrDoF,max(modes));
    X(free,:) = V;
    x = X(:,m);
    lambda4lvl(end+1) = lambda4mode(m);
    % estimate
    eta = EstimatorForm(['h^2*(f(x)+div(alphac(x)grad(u)))^2 dx', ...
        '+ h^1*([alpha(x)*grad(u)*n])^2 dsi']);
    eta4e = zeros(mesh.nrElems,1);
    for iMode = modes
        etaMode4e = eta.assemble(X(:,iMode),fe, ...
            'f',@(pts,ptsRef,parts)(lambda4mode(iMode)*fe.getValue(pts,ptsRef,parts,X(:,iMode))), ...
            'alphac',{@p.alpha,@p.dalpha},'alpha',@p.alpha);
        eta4e = eta4e + etaMode4e.^2;
    end
    eta4e = sqrt(eta4e);
    eta4lvl(end+1) = norm(eta4e);
    % errors
    a4l = zeros(size(pairs4l,1),1);
    for iPair = 1:size(pairs4l,1)
        l = LinearForm('f(x) v dx');
        b4l = l.assemble(fe,2*degree+2,'f',@(y)p.exactU(y,pairs4l(iPair,1),pairs4l(iPair,2)));
        a4l(iPair) = b4l'*x;  % (u_h,u_l)
    end
    uStar = @(y)sumModes(@p.exactU,y,a4l,pairs4l);
    gradUStar = @(y)sumModes(@p.exactGradU,y,a4l,pairs4l);
    lambdaError4lvl(end+1) = lambda4lvl(end) - p.eigenvalue(m);
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(end+1) = l2.assemble(x,fe,2*degree+2,'exact',uStar);
    h1 = ErrorForm('( exact(x) - grad(u) )^2 dx');
    H1Error4lvl(end+1) = h1.assemble(x,fe,2*degree+2,'exact',gradUStar);
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
%[text] The eigenfunction errors are measured against $u^\*$; the last column compares the eigenvalue error with $\\eta^2$.
table(nrDoF4lvl',lambda4lvl',lambdaError4lvl',L2Error4lvl',H1Error4lvl',(eta4lvl.^2)',(lambdaError4lvl./eta4lvl.^2)', ...
    VariableNames=["N" "lambda_h" "lambda_h - lambda" "L2 error" "H1 error" "eta^2" "error/eta^2"])
%[text] The discrete eigenvalues approach $5\\pi^2$ from above, and the ratio of the eigenvalue error and $\\eta^2$ stays bounded.
%%
%[text] ### Adaptive mesh
clf
plot(mesh);
%[text] The refinement is strongest in the twisted disk, where $\\alpha$ varies and the eigenfunctions are twisted.
%%
%[text] ### Discrete eigenfunction
clf
plot(fe,x);
view([40,35]);
colorbar;
shading interp;
%[text] The eigenfunction is a twisted combination of $2\\sin(\\pi y\_1)\\sin(2\\pi y\_2)$ and $2\\sin(2\\pi y\_1)\\sin(\\pi y\_2)$ in the coordinates $y = \\Phi(x)$.
%%
%[text] ### Condition number of the diffusion matrix
%[text] Since $\\det\\alpha = 1$, the condition number of $\\alpha$ is the square of its largest eigenvalue, plotted at the nodes of the final mesh.
clf
alpha4n = p.alpha(mesh.c4n);
lambdaMax4n = (alpha4n(:,1,1)+alpha4n(:,2,2))/2 ...
    + sqrt(((alpha4n(:,1,1)-alpha4n(:,2,2))/2).^2 + alpha4n(:,1,2).^2);
trisurf(mesh.n4e,mesh.c4n(:,1),mesh.c4n(:,2),lambdaMax4n.^2);
view(2);
colorbar;
shading interp;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,lambdaError4lvl,'\lambda_h-\lambda');
hold on
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u^*-u_h||_{L^2}');
plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(u^*-u_h)||_{L^2}');
plotConvergence(nrDoF4lvl,eta4lvl.^2,'\eta^2');
hold off
%[text] The eigenvalue error and $\\eta^2$ decrease with the rate $k$, the $L^2$ and the energy error of the eigenfunction with the rates $(k+1)/2$ and $k/2$.
%%
%[text] ## Possibilities for extensions
%[text] - Set `mark = MarkUniform;` to refine uniformly and compare the errors for the same number of degrees of freedom.
%[text] - Set `m = 1` for the simple eigenvalue $2\\pi^2$: the cluster has one mode, and $u^\*$ only fixes the sign of $u\_h$.
%[text] - Set `theta0 = 0`: then $\\alpha = I$, and the problem is the Laplace eigenvalue problem on the square with the eigenfunctions $2\\sin(j\\pi x\_1)\\sin(k\\pi x\_2)$. A larger angle, e.g. `theta0 = pi`, increases the anisotropy of $\\alpha$.
%[text] - Set `degree = 2`: the eigenvalue error decreases like $N^{-2}$. \
%%
%[text] ## Local functions
function val = sumModes(fun,y,a4l,pairs4l)
% linear combination sum_l a_l fun(y,j_l,k_l) of eigenfunctions
val = 0;
for l = 1:size(pairs4l,1)
    val = val + a4l(l)*fun(y,pairs4l(l,1),pairs4l(l,2));
end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
