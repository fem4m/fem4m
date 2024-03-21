%[text] # An unstable mixed $P\_1^2\\times P\_0$ discretization of the Poisson problem
%[text] Discretizes the Neumann problem for the Poisson equation in mixed form with continuous piecewise linear fluxes and piecewise constant solutions. This pair of spaces violates the discrete inf-sup condition: the linear system is singular and in general has no solution, and the discrete flux does not converge. The stable Raviart-Thomas pair $RT\_0\\times P\_0$ on the same meshes converges.
%[text] **Builds on:** `exPoissonP1`. **Next:** `exPoissonRTAdaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] On the unit square $\\Omega = (0,1)^2$ the Neumann problem reads: given $f$ with $\\int\_\\Omega f\\,dx = 0$, find $u$ with
%[text]{"align":"center"} $ -\\Delta u = f \\textrm{ in } \\Omega, \\qquad \\nabla u\\cdot n = 0 \\textrm{ on } \\partial\\Omega, \\qquad \\int\_\\Omega u\\,dx = 0. $
%[text] Without the last condition $u$ would be unique only up to a constant. The problem object `problemScalarEllipticSquareNeumannExact` provides the source $f = -\\Delta u$ of the exact solution
%[text]{"align":"center"} $ u(x) = x\_1^2(1-x\_1)^2\\,x\_2^2(1-x\_2)^2 - 1/900, $
%[text] whose normal derivative vanishes on $\\partial\\Omega$ and whose mean value is zero.
%[text] ### Mixed formulation
%[text] With the flux $\\sigma = \\nabla u$ the problem becomes $\\sigma - \\nabla u = 0$, $\\mathrm{div}\\,\\sigma = -f$ in $\\Omega$ and $\\sigma\\cdot n = 0$ on $\\partial\\Omega$. As in `exPoissonRTAdaptive` the Neumann condition is essential: the flux lies in $V = H\_0(\\mathrm{div},\\Omega) = \\{\\tau\\in H(\\mathrm{div},\\Omega) : \\tau\\cdot n = 0 \\textrm{ on } \\partial\\Omega\\}$. The mean value of $u$ is fixed by a Lagrange multiplier $\\lambda\\in\\mathbb{R}$. Find $\\sigma\\in V$, $u\\in L^2(\\Omega)$ and $\\lambda\\in\\mathbb{R}$ with
%[text]{"align":"center"} $ \\int\_\\Omega \\sigma\\cdot\\tau\\,dx + \\int\_\\Omega u\\,\\mathrm{div}\\,\\tau\\,dx = 0 \\quad \\textrm{for all } \\tau\\in V, $
%[text]{"align":"center"} $ \\int\_\\Omega \\mathrm{div}\\,\\sigma\\,v\\,dx + \\lambda\\int\_\\Omega v\\,dx = -\\int\_\\Omega f v\\,dx \\quad \\textrm{for all } v\\in L^2(\\Omega), \\qquad \\int\_\\Omega u\\,dx = 0. $
%[text] Testing with $v = 1$ shows $\\lambda = 0$, because $\\int\_\\Omega\\mathrm{div}\\,\\sigma\\,dx = \\int\_{\\partial\\Omega}\\sigma\\cdot n\\,ds = 0$ and $\\int\_\\Omega f\\,dx = 0$.
%[text] ### The inf-sup condition
%[text] By the theory of Brezzi this saddle point problem, and also its discretization in subspaces $V\_h\\subset V$ and $Q\_h\\subset L^2(\\Omega)$, has a unique and stable solution if two conditions hold. First, $\\int\_\\Omega\\sigma\\cdot\\tau\\,dx$ must be coercive in the norm of $H(\\mathrm{div},\\Omega)$ on the discrete kernel $\\{\\tau\_h\\in V\_h : \\int\_\\Omega \\mathrm{div}\\,\\tau\_h\\,v\_h\\,dx = 0 \\textrm{ for all } v\_h\\in Q\_h\\}$; this holds whenever $\\mathrm{div}\\,V\_h\\subset Q\_h$, because then the kernel consists of divergence-free fields. Second, the discrete inf-sup condition
%[text]{"align":"center"} $ \\beta\_h = \\inf\_{v\_h\\in Q\_h\\cap L^2\_0(\\Omega)\\setminus\\{0\\}}\\ \\sup\_{\\tau\_h\\in V\_h\\setminus\\{0\\}} \\frac{\\int\_\\Omega \\mathrm{div}\\,\\tau\_h\\,v\_h\\,dx}{\\|\\tau\_h\\|\_{H(\\mathrm{div},\\Omega)}\\|v\_h\\|\_{L^2(\\Omega)}} \\ge \\beta\_0 \> 0 $
%[text] must hold with a constant $\\beta\_0$ independent of the mesh, where $L^2\_0(\\Omega)$ are the functions with mean value zero. It requires the discrete divergence to map $V\_h$ onto $Q\_h\\cap L^2\_0(\\Omega)$, so $V\_h$ must be rich enough compared to $Q\_h$.
%[text] ### Discretization
%[text] This example takes the continuous piecewise linear vector fields for $\\sigma\_h$ (`PkVector(1,mesh)`), which are in $H(\\mathrm{div},\\Omega)$, and the piecewise constants for $u\_h$ (`DGDk(0,mesh)`). The geometry with the slip boundary `'S'` stores the coordinate direction of the normal of every boundary edge in `mesh.n4Sb(:,3)`; the essential condition $\\tau\_h\\cdot n = 0$ fixes this component of $\\tau\_h$ at the two nodes of the edge, both components at the corners. The forms are those of `exPoissonRTAdaptive`, plus the mean value:
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega \\sigma\\cdot\\tau\\,dx$ | `sigma*tau dx` |
%[text] | $\\int\_\\Omega \\mathrm{div}\\,\\sigma\\,v\\,dx$ | `div(sigma)*v dx` |
%[text] | $\\int\_\\Omega f v\\,dx$ | `f(x) v dx` |
%[text] | $\\int\_\\Omega v\\,dx$ | `v dx` |
%[text:table]
%[text] The divergence of a piecewise linear field is piecewise constant, so $\\mathrm{div}\\,V\_h\\subset Q\_h$ and the first condition holds. The inf-sup condition fails. The coarse mesh cuts the square by the diagonal from the lower left to the upper right corner, and uniform refinement keeps this pattern: the mesh of level $\\ell$ consists of $n\\times n$ squares with $n = 2^\\ell$, each cut into two triangles. Of the $2(n+1)^2$ nodal values of $\\tau\_h$ the boundary condition fixes $4n+4$, so
%[text]{"align":"center"} $ \\dim V\_h = 2n^2 - 2 \< 2n^2 - 1 = \\dim(Q\_h\\cap L^2\_0(\\Omega)). $
%[text] The divergence cannot map onto $Q\_h\\cap L^2\_0(\\Omega)$: there are spurious modes $v\_h\\ne 0$ with mean value zero and $\\int\_\\Omega \\mathrm{div}\\,\\tau\_h\\,v\_h\\,dx = 0$ for all $\\tau\_h\\in V\_h$, and $\\beta\_h = 0$. A spurious mode can be added to $u\_h$ without changing the equations, and the second equation has a solution only if $f$ is orthogonal to all spurious modes: the system matrix is singular, and for a general $f$ the system has no solution at all. The code computes the minimum-norm least-squares solution with `lsqminnorm`.
%[text] The Raviart-Thomas pair $RT\_0\\times P\_0$ (`RTk(1,mesh)` and `DGDk(0,mesh)`) has one degree of freedom per edge, the normal flux, which is zero on the boundary edges. The divergence maps $RT\_0$ onto the piecewise constants with mean value zero, with a right inverse bounded independently of the mesh (Fortin operator), so the inf-sup condition holds.
%[text] ### Expected convergence
%[text] The stable pair $RT\_0\\times P\_0$ converges like $\\|\\sigma-\\sigma\_h\\|\_{L^2(\\Omega)} + \\|u-u\_h\\|\_{L^2(\\Omega)} = O(h) = O(N^{-1/2})$. For $P\_1^2\\times P\_0$ the theory gives no convergence; the results show that the flux error does not decrease.
%%
%[text] ## Parameters
nrLevels = 5;  % number of meshes
%%
%[text] ## Problem data and initial mesh
%[text] The slip boundary `'S'` replaces the Neumann boundary of `p.geom`; `'criss'` cuts the square by the diagonal from the lower left to the upper right corner. Both pairs use the same meshes.
p = problemScalarEllipticSquareNeumannExact();
mesh = Mesh(Geometry(struct('name','Square','bdry','S','reftype','criss')));
mesh = mesh.refine('uniform');
plot(mesh);
%%
%[text] ## Convergence study on uniform meshes
%[text] On every mesh the loop solves both discretizations, computes the rank of the matrix $B$ of the discrete divergence of $P\_1^2\\times P\_0$ and the errors, and refines all triangles uniformly.
nrDoF4lvl = zeros(1,nrLevels);
nrDoFRT4lvl = zeros(1,nrLevels);
dimV4lvl = zeros(1,nrLevels);
dimQ4lvl = zeros(1,nrLevels);
rankB4lvl = zeros(1,nrLevels);
L2Error4lvl = zeros(1,nrLevels);
sigmaError4lvl = zeros(1,nrLevels);
L2ErrorRT4lvl = zeros(1,nrLevels);
sigmaErrorRT4lvl = zeros(1,nrLevels);
for lvl = 1:nrLevels
    % finite element spaces
    feSigma = PkVector(1,mesh);
    feSigmaRT = RTk(1,mesh);
    feU = DGDk(0,mesh);
    ed4Sb = mesh.ed4Sb;
    isNormalX = mesh.n4Sb(:,3) == 1;  % normal +-e1: tau_1 = 0
    isNormalY = mesh.n4Sb(:,3) == 2;  % normal +-e2: tau_2 = 0
    fixed = unique([feSigma.dof4ed(ed4Sb(isNormalX),1:feSigma.degree+1), ...
        feSigma.dof4ed(ed4Sb(isNormalY),feSigma.degree+2:end)]);
    free = setdiff(1:feSigma.nrDoF,fixed);
    freeRT = setdiff(1:feSigmaRT.nrDoF,feSigmaRT.dof4ed(ed4Sb,:));
    nrDoF4lvl(lvl) = numel(free) + feU.nrDoF;
    nrDoFRT4lvl(lvl) = numel(freeRT) + feU.nrDoF;
    dimV4lvl(lvl) = numel(free);
    dimQ4lvl(lvl) = feU.nrDoF;
    % assemble
    a = BilinearForm('sigma*tau dx',{'sigma','tau'});
    bDiv = BilinearForm('div(sigma)*v dx',{'sigma','v'});
    A = a.assemble(feSigma);
    B = bDiv.assemble(feSigma,feU);
    ART = a.assemble(feSigmaRT);
    BRT = bDiv.assemble(feSigmaRT,feU);
    l = LinearForm('f(x) v dx');
    b = l.assemble(feU,'f',@p.f);
    m = LinearForm('v dx');
    c = m.assemble(feU);
    % rank of the discrete divergence
    rankB4lvl(lvl) = rank(full(B(:,free)));
    % solve: least squares for the singular P1 x P0 system
    nrFree = numel(free);
    S = [A(free,free), B(:,free)', sparse(nrFree,1);
        B(:,free), sparse(feU.nrDoF,feU.nrDoF), c;
        sparse(1,nrFree), c', 0];
    y = lsqminnorm(S,[zeros(nrFree,1); -b; 0]);
    sigma = zeros(feSigma.nrDoF,1);
    sigma(free) = y(1:nrFree);
    u = y(nrFree+(1:feU.nrDoF));
    nrFreeRT = numel(freeRT);
    SRT = [ART(freeRT,freeRT), BRT(:,freeRT)', sparse(nrFreeRT,1);
        BRT(:,freeRT), sparse(feU.nrDoF,feU.nrDoF), c;
        sparse(1,nrFreeRT), c', 0];
    yRT = SRT\[zeros(nrFreeRT,1); -b; 0];
    sigmaRT = zeros(feSigmaRT.nrDoF,1);
    sigmaRT(freeRT) = yRT(1:nrFreeRT);
    uRT = yRT(nrFreeRT+(1:feU.nrDoF));
    % errors
    l2 = ErrorForm('( exact(x) - u )^2 dx');
    l2Sigma = ErrorForm('( exact(x) - sigma )^2 dx','sigma');
    L2Error4lvl(lvl) = l2.assemble(u,feU,8,'exact',@p.exactU);
    sigmaError4lvl(lvl) = l2Sigma.assemble(sigma,feSigma,8,'exact',@p.exactSigma);
    L2ErrorRT4lvl(lvl) = l2.assemble(uRT,feU,8,'exact',@p.exactU);
    sigmaErrorRT4lvl(lvl) = l2Sigma.assemble(sigmaRT,feSigmaRT,8,'exact',@p.exactSigma);
    % refine
    mesh = mesh.refine('uniform');
end
%%
%[text] ## Results
%[text] ### Spurious modes
%[text] The number of spurious modes of $P\_1^2\\times P\_0$ is $\\dim(Q\_h\\cap L^2\_0(\\Omega)) - \\mathrm{rank}\\,B$, the dimension of the space of the mean-free $v\_h$ with $B^Tv\_h = 0$.
table(dimV4lvl',dimQ4lvl',rankB4lvl',(dimQ4lvl-1-rankB4lvl)', ...
    VariableNames=["dim V_h" "dim Q_h" "rank B" "spurious modes"])
%[text] The rank of $B$ falls short of $\\dim Q\_h - 1$ by more than the dimension count requires: the number of spurious modes doubles with every refinement, like $n$. The kernel of $B$ in $V\_h$, the fields with $\\mathrm{div}\\,\\tau\_h = 0$, has the dimension $\\dim V\_h - \\mathrm{rank}\\,B$, which is one less than the number of spurious modes. The second equation thus leaves $\\sigma\_h$ only these few degrees of freedom: the divergence constraint almost determines the discrete flux, which cannot approximate $\\sigma$ (locking).
%%
%[text] ### The checkerboard mode
%[text] One spurious mode is the checkerboard $v\_h = 1$ on the upper and $v\_h = -1$ on the lower triangle of every square; the upper triangle has two vertices on the upper edge of its square. Around an interior node the six triangles alternate between upper and lower ones, and the alternating sum of the gradients of the hat function of the node vanishes; the same holds for the tangential derivative at a boundary node. Hence $\\int\_\\Omega \\mathrm{div}\\,\\tau\_h\\,v\_h\\,dx = 0$ for every $\\tau\_h\\in V\_h$, and $B^Tv\_h$ vanishes on the finest mesh:
y4e = reshape(feU.mesh.c4n(feU.mesh.n4e,2),[],3);
isUpper = sum(y4e == max(y4e,[],2),2) == 2;
checkerboard = 2*isUpper - 1;
norm(B(:,free)'*checkerboard)
clf
plot(feU,checkerboard);
view(2);
colorbar;
%[text] Adding a multiple of the checkerboard to $u\_h$ does not change the discrete equations.
%%
%[text] ### Discrete solutions
clf
tiledlayout(1,2);
nexttile
plot(feU,u);
view([40,35]);
title('P_1^2\times P_0');
nexttile
plot(feU,uRT);
view([40,35]);
title('RT_0\times P_0');
%[text] The least-squares solution of $P\_1^2\\times P\_0$ is much larger than the exact solution and has a different shape, while $RT\_0\\times P\_0$ approximates $u$ by piecewise constants.
%%
%[text] ### Convergence table
table(nrDoF4lvl',sigmaError4lvl',L2Error4lvl',nrDoFRT4lvl',sigmaErrorRT4lvl',L2ErrorRT4lvl', ...
    VariableNames=["N P1xP0" "sigma error P1xP0" "u error P1xP0" "N RT0xP0" "sigma error RT0xP0" "u error RT0xP0"])
%%
%[text] ### Convergence
clf
plotConvergence(nrDoF4lvl,sigmaError4lvl,'||\sigma-\sigma_h||_{L^2}, P_1^2\times P_0');
hold on
plotConvergence(nrDoF4lvl,L2Error4lvl,'||u-u_h||_{L^2}, P_1^2\times P_0');
plotConvergence(nrDoFRT4lvl,sigmaErrorRT4lvl,'||\sigma-\sigma_h||_{L^2}, RT_0\times P_0');
plotConvergence(nrDoFRT4lvl,L2ErrorRT4lvl,'||u-u_h||_{L^2}, RT_0\times P_0');
hold off
legend('Location','southwest');
%[text] Both errors of $RT\_0\\times P\_0$ decay like $N^{-1/2}$. The flux error of $P\_1^2\\times P\_0$ does not decrease, and its error in $u$ even grows: a discretization that violates the inf-sup condition does not converge, although its spaces approximate $\\sigma$ and $u$ well.
%%
%[text] ## Possibilities for extensions
%[text] - Use the coarse mesh with both diagonals of the square, `'reftype','crisscross'`, and remove the checkerboard section: the pair is unstable on these meshes as well, the table still shows spurious modes, but the checkerboard is no longer one of them.
%[text] - Compute the discrete inf-sup constant $\\beta\_h$, the square root of the second smallest generalized eigenvalue of $B\_fM\_V^{-1}B\_f^T$ and the mass matrix $M\_Q$ of $Q\_h$ (the smallest one, zero, belongs to the constants), with the matrix $M\_V$ of the scalar product of $H(\\mathrm{div},\\Omega)$: `MV = A + BilinearForm('div(sigma)*div(tau) dx',{'sigma','tau'}).assemble(feSigma); MQ = BilinearForm('u*v dx').assemble(feU); ev = sort(eig(full(B(:,free)*(MV(free,free)\B(:,free)')),full(MQ)));` and `sqrt(abs(ev(2)))`. For $P\_1^2\\times P\_0$ the eigenvalue `ev(2)` vanishes up to rounding, so $\\beta\_h = 0$; for $RT\_0\\times P\_0$ (with `ART`, `BRT`, `feSigmaRT` and `freeRT`) $\\beta\_h$ is bounded away from zero, independently of the mesh.
%[text] - Replace $RT\_0\\times P\_0$ by $RT\_1\\times P\_1$, `RTk(2,mesh)` with `DGDk(1,mesh)`: both errors decay like $N^{-1}$.
%[text] - The same pair is unstable for the Stokes problem, where $P\_1^2\\times P\_0$ locks; `exStokesTaylorHood` uses the stable Taylor-Hood pair.
%[text] - `exPoissonRTAdaptive` solves the Poisson problem with Raviart-Thomas elements of any degree on adaptive meshes. \

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---
