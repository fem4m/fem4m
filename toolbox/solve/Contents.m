% SOLVE  Multigrid, Schwarz preconditioners, active set and inverse iteration
%
%   Classes
%   Multigrid    -  Geometric multigrid solver
%   Schwarz      -  Base class of additive and multiplicative Schwarz methods
%   SchwarzBDDC  -  BDDC preconditioner
%   SchwarzBPS   -  BPS preconditioner
%   SchwarzBPX   -  BPX multilevel preconditioner
%   SchwarzHB    -  Hierarchical basis preconditioner
%   SchwarzNN    -  Neumann-Neumann preconditioner
%   SchwarzVO    -  Vertex patch preconditioner
%
%   Functions
%   pdaSet       -  Primal-dual active set algorithm
%   pinvit       -  Inexact preconditioned inverse iteration