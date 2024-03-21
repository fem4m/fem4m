%[text] # Mixed DG method with numerical fluxes for the Poisson problem
%[text] Discretizes the Poisson problem in mixed form with discontinuous approximations of the flux and of the solution that are coupled across the edges by numerical fluxes. Three flux parameters span the family of mixed DG methods of the unified framework of Arnold, Brezzi, Cockburn and Marini; the example compares six choices and shows how they determine the convergence of the flux.
%[text] **Builds on:** `exPoissonDG`, `exPoissonLDG`. **Next:** `exPoissonHDG`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] On the unit square $\\Omega = (0,1)^2$ the Poisson problem $-\\Delta u = f$ with $u = 0$ on $\\partial\\Omega$ reads in mixed form, with the flux $\\sigma = -\\nabla u$,
%[text]{"align":"center"} $ \\sigma + \\nabla u = 0, \\qquad \\mathrm{div}\\,\\sigma = f \\quad \\textrm{in } \\Omega, \\qquad u = 0 \\quad \\textrm{on } \\partial\\Omega. $
%[text] The problem object `problemScalarEllipticSquareExact` provides the source of the exact solution $u = x\_1(1-x\_1)\\,x\_2(1-x\_2)$.
%[text] ### Numerical fluxes
%[text] The notation for jumps $\[\\cdot\]$, averages $\\{\\cdot\\}$, normals $n\_E$ and the edges $\\mathcal{E}$, $\\mathcal{E}(\\Omega)$ is that of `exPoissonDG`. The flux $\\sigma\_h$ and the solution $u\_h$ are discontinuous piecewise polynomials of degree $k$, $\\sigma\_h\\in P\_k(\\mathcal{T})^2$ (`DGDkVector(degree,mesh)`) and $u\_h\\in P\_k(\\mathcal{T})$ (`DGDk(degree,mesh)`). As in `exPoissonLDG`, integration by parts on every triangle with numerical fluxes $\\hat u\_h$ and $\\hat\\sigma\_h$ in place of the traces gives, summed over the triangles,
%[text]{"align":"center"} $ (\\sigma\_h,\\tau)\_{\\mathcal{T}} - (u\_h,\\mathrm{div}\\,\\tau)\_{\\mathcal{T}} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} (\\hat u\_h, \[\\tau\]\\cdot n\_E)\_E = 0, \\qquad -(\\sigma\_h,\\nabla\_h v)\_{\\mathcal{T}} + \\sum\_{E\\in\\mathcal{E}} (\\hat\\sigma\_h\\cdot n\_E, \[v\])\_E = (f,v) $
%[text] for all test functions $\\tau$ and $v$; on the boundary $\\hat u\_h = 0$. The unified framework collects the consistent fluxes of the form
%[text]{"align":"center"} $ \\hat u\_h = \\{u\_h\\} - (C\_{12}\\cdot n\_E)\[u\_h\] + C\_{22}\[\\sigma\_h\]\\cdot n\_E, \\qquad \\hat\\sigma\_h = \\{\\sigma\_h\\} + C\_{11}\[u\_h\]\\,n\_E + C\_{12}\[\\sigma\_h\]\\cdot n\_E $
%[text] on the interior edges and $\\hat\\sigma\_h = \\sigma\_h + C\_{11}u\_h\\,n$ on the boundary edges, with the parameters $C\_{11} \> 0$, $C\_{22}\\ge 0$ and the vector $C\_{12}$; for $C\_{11} = 0$ the linear system below is singular. $C\_{11}$ penalizes the jumps of $u\_h$, $C\_{22}$ the jumps of the normal component of $\\sigma\_h$. $C\_{12}$ shifts both fluxes towards one side of the edge: for $C\_{12}\\cdot n\_E = 1/2$, $\\hat u\_h$ is the trace of $u\_h$ from $T\_-$ and $\\hat\\sigma\_h\\cdot n\_E$ is the trace of $\\sigma\_h\\cdot n\_E$ from $T\_+$ plus the penalty $C\_{11}\[u\_h\]$. This example uses $C\_{12} = 0$. The parameters may depend on the edge length $h\_E$; in the code they are functions of the vector `h` of the edge lengths.
%[text] ### Linear system
%[text] Inserting the fluxes and integrating $(\\sigma\_h,\\nabla\_h v)\_{\\mathcal{T}}$ back by parts gives the system
%[text]{"align":"center"} $ \\left(\\begin{array}{cc} A & -B^T \\\\ B & S \\end{array}\\right) \\left(\\begin{array}{c} \\sigma \\\\ u \\end{array}\\right) = \\left(\\begin{array}{c} 0 \\\\ b \\end{array}\\right) $
%[text] for the coefficient vectors $\\sigma$ and $u$, with the matrices of the forms
%[text:table]
%[text] | Matrix | Integral | fem4m form |
%[text] | --- | --- | --- |
%[text] | $A$ | $(\\sigma,\\tau)\_{\\mathcal{T}} + \\sum\_{E\\in\\mathcal{E}(\\Omega)} (C\_{22}\[\\sigma\]\\cdot n\_E, \[\\tau\]\\cdot n\_E)\_E$ | `sigma*tau dx + C(x)*[sigma]*n [tau]*n dsi` |
%[text] | $B$ | $(\\mathrm{div}\\,\\sigma, v)\_{\\mathcal{T}} - \\sum\_{E\\in\\mathcal{E}(\\Omega)} (\[\\sigma\]\\cdot n\_E, \\{v\\})\_E$ | `div(sigma) v dx - [sigma]*n{v} dsi` |
%[text] | $B$ (part) | $\\sum\_{E\\in\\mathcal{E}(\\Omega)} ((C\_{12}\\cdot n\_E)\[v\], \[\\sigma\]\\cdot n\_E)\_E$ | `C(x)[u]*n [tau]*n dsi`, transposed |
%[text] | $S$ | $\\sum\_{E\\in\\mathcal{E}} (C\_{11}\[u\], \[v\])\_E$ | `C(x)*[u][v] ds` |
%[text] | $b$ | $(f,v)$ | `f(x) v dx` |
%[text:table]
%[text] With a vector-valued coefficient, `C(x)[u]*n` is $(C\_{12}\\cdot n\_E)\[u\]$. For $C\_{22} = 0$ the matrix $A$ is the block diagonal mass matrix and $\\sigma$ can be eliminated triangle by triangle: these are the LDG methods of `exPoissonLDG` (there with $\\sigma\_h$ of degree $k-1$). For $C\_{22} \> 0$ the matrix $A$ couples neighbouring triangles.
%[text] ### The six choices
%[text] The parameters below are the alternatives of the section Parameters. For $k\\ge 1$ the computations on the uniform meshes of this example give the following orders of convergence in $h$; for $C\_{22} = 0$ they are the orders of the a priori analysis. In terms of the number $N$ of unknowns the rates are half of them.
%[text:table]
%[text] | $C\_{11}$ | $C\_{22}$ | $\\Vert u-u\_h\\Vert\_{L^2(\\Omega)}$ | $\\Vert\\sigma-\\sigma\_h\\Vert\_{L^2(\\Omega)}$ |
%[text] | --- | --- | --- | --- |
%[text] | $1$ or $1/h$ | $0$ (LDG) | $k+1$ | $k$ |
%[text] | $1$ or $1/h$ | $h$ | $k+1$ | $k$ |
%[text] | $1/h$ | $1$ | $k+1$ | $k$ for $k\\ge 2$, about $k+1/2$ for $k = 1$ |
%[text] | $1$ | $1$ | $k+1$ | $k+1$ |
%[text:table]
%[text] On coarse meshes the flux may converge faster than with the order $k$. Only the choice with both parameters of order one, the default, approximates the flux with the full order $k+1$. It resembles the HDG method of `exPoissonHDG`: eliminating the trace unknown $\\lambda\_h = \\{u\_h\\} + \\frac{1}{2\\tau}\[\\sigma\_h\]\\cdot n\_E$ from the HDG method with the stabilization $\\tau$ gives exactly the fluxes above with $C\_{11} = \\tau/2$ on the interior and $C\_{11} = \\tau$ on the boundary edges, $C\_{22} = 1/(2\\tau)$ and $C\_{12} = 0$.
%[text] For $k = 0$ the gradient of $u\_h$ vanishes on every triangle and all information sits in the jumps. The choices with $C\_{11} = 1$ converge with the order $1$ in both variables. The jumps of a piecewise constant approximation of $u$ are of order $h$, and with $C\_{11} = 1/h$ their penalty $\\sum\_{E\\in\\mathcal{E}} C\_{11}\\|\[u\_h\]\\|^2\_{L^2(E)}$ does not become small: these methods do not converge for $k = 0$.
%[text] ### Expected convergence
%[text] The default is $k = 0$ with $C\_{11} = C\_{22} = 1$: both errors decay with the order $1$ in $h$, i.e. like $N^{-1/2}$.
%%
%[text] ## Parameters
degree = 0;          % polynomial degree k of sigma_h and u_h
C11 = @(h) 1;        % penalty of the jumps of u_h, a function of the edge lengths h
% C11 = @(h) 1./h;   % penalty of interior penalty type: no convergence for k = 0
C22 = @(h) 1;        % penalty of the normal jumps of sigma_h
% C22 = @(h) 0;      % LDG: sigma_h is eliminated locally, the flux converges with order k
% C22 = @(h) h;      % small penalty: the flux converges with order k
C12 = [0 0];         % vector that shifts the fluxes towards one side of the edges
nrLevels = 6;        % number of meshes
%%
%[text] ## Problem data and initial mesh
p = problemScalarEllipticSquareExact();
mesh = Mesh(p.geom);
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
nrDoF4lvl = zeros(1,nrLevels);
L2Error4lvl = zeros(1,nrLevels);
sigmaError4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element spaces
    feSigma = DGDkVector(degree,mesh);
    feU = DGDk(degree,mesh);
    nrDoF4lvl(lvl) = feSigma.nrDoF + feU.nrDoF;
    % assemble
    h = mesh.length4ed;
    a = BilinearForm('sigma*tau dx + C(x)*[sigma]*n [tau]*n dsi',{'sigma','tau'});
    A = a.assemble(feSigma,'C',C22(h));
    bDiv = BilinearForm('div(sigma) v dx - [sigma]*n{v} dsi',{'sigma','v'});
    B = bDiv.assemble(feSigma,feU);
    c12 = BilinearForm('C(x)[u]*n [tau]*n dsi',{'u','tau'});
    B = B + c12.assemble(feU,feSigma,'C',C12)';
    s = BilinearForm('C(x)*[u][v] ds');
    S = s.assemble(feU,'C',C11(h));
    l = LinearForm('f(x) v dx');
    b = l.assemble(feU,'f',@p.f);
    % solve
    x = [A, -B'; B, S]\[zeros(feSigma.nrDoF,1); b];
    sigma = x(1:feSigma.nrDoF);
    u = x(feSigma.nrDoF+1:end);
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    L2Error4lvl(lvl) = l2.assemble(u,feU,2*degree+2,'exact',@p.exactU);
    l2Sigma = ErrorForm('( exact(x) + sigma )^2 dx','sigma');
    sigmaError4lvl(lvl) = l2Sigma.assemble(sigma,feSigma,2*degree+2,'exact',@p.exactGradU);
    % refine
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Convergence table
L2Rate = [NaN, -diff(log(L2Error4lvl))./diff(log(nrDoF4lvl))];
sigmaRate = [NaN, -diff(log(sigmaError4lvl))./diff(log(nrDoF4lvl))];
table(nrDoF4lvl',L2Error4lvl',L2Rate',sigmaError4lvl',sigmaRate', ...
    VariableNames=["N" "u error" "u rate" "sigma error" "sigma rate"])
%[text] Both rates approach $1/2$ with respect to $N$, i.e. the order $1 = k+1$ in $h$.
%%
%[text] ### Discrete solution
clf
plot(feU,u);
view([40,35]);
colorbar;
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}');
hold on
plotConvergence(nrDoF4lvl,sigmaError4lvl,'||\sigma-\sigma_h||_{L^2}');
hold off
%%
%[text] ## Possibilities for extensions
%[text] - Set `degree = 1`: with the default parameters both errors converge with the order 2 in $h$, i.e. like $N^{-1}$.
%[text] - Set `degree = 1` and `C22 = @(h) 0`: the LDG method; the error of the flux drops to the order 1, i.e. $N^{-1/2}$, while the error of $u\_h$ keeps the order 2.
%[text] - Set `C11 = @(h) 1./h` with `degree = 0`: the errors do not converge, the error of $u\_h$ even grows on the finer meshes.
%[text] - Choose the parameters of the HDG method with $\\tau = 1$: `C22 = @(h) 1/2` and, in the assembly of `S`, the vector `c11` with `c11 = 1/2*ones(mesh.nrEdges,1); c11(mesh.ed4Db) = 1;` in place of `C11(h)`. The solution coincides with that of `exPoissonHDG` for `degree = degreeU` and `tau = 1`, which computes it from a smaller system for the trace unknowns.
%[text] - Set `degree = 1` and `C12 = [1 1]/2`: the fluxes become one-sided on the horizontal and vertical edges, the rates do not change. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
