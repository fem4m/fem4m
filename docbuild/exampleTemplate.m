%[text] # \<Adaptive\> \<method\> for the \<problem\>
%[text] \<One or two sentences: what is solved, with which method, and what this example shows that the previous ones do not. The overview of the examples reuses this paragraph.\>
%[text] **Builds on:** `exPoissonP1`. **Next:** `exPoissonPkAdaptive`.
%[text]{"align":"right"} Copyright (c) 2026 by the fem4m authors
%[text:tableOfContents]{"heading":"Contents"}
%%
%[text] ## Introduction
%[text] ### The problem
%[text] \<Domain, data and boundary conditions of the problem object used below, e.g. `problemScalarEllipticLshapeExact`, and its exact solution.\>
%[text]{"align":"center"} $ -\\Delta u = f \\textrm{ in } \\Omega, \\qquad u = u\_D \\textrm{ on } \\Gamma\_D, \\qquad \\nabla u\\cdot n = g \\textrm{ on } \\Gamma\_N $
%[text] ### Weak formulation
%[text] \<Spaces and the variational problem; mixed or DG formulations here.\>
%[text] ### Discretization
%[text] \<Discrete spaces with their fem4m classes, the discrete problem or the linear system, and the integrals next to their form strings.\>
%[text:table]
%[text] | Integral | fem4m form |
%[text] | --- | --- |
%[text] | $\\int\_\\Omega \\nabla u\\cdot\\nabla v\\,dx$ | `grad(u) * grad(v) dx` |
%[text] | $\\int\_\\Omega f v\\,dx + \\int\_{\\Gamma\_N} g v\\,ds$ | `f(x) v dx + g(x) v dsn` |
%[text:table]
%[text] ### Error estimator
%[text] \<Adaptive examples: the estimator, its terms as form strings, the marking. Solver examples: a section Solver instead.\>
%[text] ### Expected convergence
%[text] \<Rates in terms of the number of degrees of freedom N; the Results section compares them with the computed ones.\>
%%
%[text] ## Parameters
%[text] \<Every value a reader may change, one per line with a comment. Alternatives only here, as commented lines with a reason.\>
% degree = 1;          % polynomial degree k
% minNrDoF = 1e4;      % stop as soon as the number of degrees of freedom exceeds minNrDoF
% mark = MarkBulk;     % bulk (Doerfler) marking
% % mark = MarkUniform;  % uniform refinement, for comparison
%%
%[text] ## Problem data and initial mesh
% p = problemScalarEllipticLshapeExact();
% mesh = Mesh(p.geom);
% plot(mesh);
%%
%[text] ## Adaptive loop
%[text] \<Or: Convergence study on uniform meshes. One loop, no output inside, the steps in this order with these comments.\>
% nrDoF4lvl = [];
% H1Error4lvl = [];
% eta4lvl = [];
% while true
%     % finite element space
%     % assemble
%     % solve
%     % estimate
%     % errors
%     % stop
%     % mark
%     % refine
% end
%%
%[text] ## Results
%[text] ### Convergence table
% table(nrDoF4lvl',H1Error4lvl',eta4lvl',VariableNames=["N" "H1 error" "eta"])
%%
%[text] ### Adaptive mesh
% clf
% plot(mesh);
%%
%[text] ### Discrete solution
% clf
% plot(fe,x);
%%
%[text] ### Convergence
% clf
% plotConvergence(nrDoF4lvl,H1Error4lvl,'||\nabla(u-u_h)||_{L^2}');
% hold on
% plotConvergence(nrDoF4lvl,eta4lvl,'\eta');
% hold off
%[text] \<Compare the observed rates (legend) with the expected ones; no numbers that depend on the run.\>
%%
%[text] ## Possibilities for extensions
%[text] - \<Concrete variation with the code to change, e.g. `mark = MarkUniform;`, and what to observe.\>
%[text] - \<Related example, e.g. `exPoissonHPAdaptive`.\> \
%%
%[text] ## Conventions (not part of an example)
%[text] - Title "\<Adaptive\> \<method\> for the \<problem\>", sentence case, at most one formula ($P\_k$). The first heading is the title of the page and of the overview.
%[text] - Display equations: `%[text]{"align":"center"} $ ... $`, LaTeX commands with two backslashes, markdown characters inside formulas with one (`\_`, `\*`, `\[`); no `$$`, no headings around formulas, no line breaks inside formulas.
%[text] - Notation: triangulation $\\mathcal{T}$, edges $\\mathcal{E}$, interior edges $\\mathcal{E}(\\Omega)$, $h\_T$, $h\_E$, $n\_E$, jump $\[\\cdot\]$, average $\\{\\cdot\\}$, $P\_k(\\mathcal{T})$, $\\Gamma\_D$, $\\Gamma\_N$, $u\_D$, $g$, number of degrees of freedom $N$. In a table write $\\Vert\\cdot\\Vert$ instead of $\\|\\cdot\\|$.
%[text] - Variables: `p` the problem, `mesh`, `fe` or `feU`, `feP`, `feSigma`, `feLambda`, `A`, `b`, `x`, `free`, `lvl`, `nrDoF4lvl`, `L2Error4lvl`, `H1Error4lvl`, `eta4lvl`; forms named by role (`a`, `l`, `l2`, `h1`, `eta`), never `error`.
%[text] - Loop steps with the comments finite element space, assemble, solve, estimate, errors, stop, mark, refine; no output inside the loop, a table after it.
%[text] - Figures: one per section, `clf` before every figure but the first, `hold on` ... `hold off`, no `figure(n)`; `plotConvergence` with TeX labels such as `'||u-u_h||_{L^2}'`.
%[text] - No unreachable branches; choose by a parameter. Local functions at the end. \

% Copyright (c) 2026 by the fem4m authors

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline"}
%---