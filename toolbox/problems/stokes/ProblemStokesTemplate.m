classdef ProblemStokesTemplate < StokesParameters
    %ProblemStokesTemplate - Abstract Stokes problem
    %   ProblemStokesTemplate is an abstract template for Stokes flow
    %   problems -NU*Delta u + grad p = f, div u = 0; subclasses supply the
    %   geometry and the data. It inherits the viscosity NU and the
    %   compliance matrices A, AWeakSym and A3d from StokesParameters.
    %
    %   The data are evaluated at the points X (N-by-d, one row per point):
    %   f and u4Db return the body force and the Dirichlet velocity as
    %   N-by-d arrays, and g(X,NORMALS) gets the outer unit normals NORMALS
    %   (N-by-d) and returns the Neumann data as an N-by-d array
    %   (ProblemStokesTemplateExact: N-by-1-by-d). The subclasses define
    %   which quantity g is (eps(u)*n or (grad u)'*n, without NU and the
    %   pressure).
    %
    %   ProblemStokesTemplate functions:
    %       f     - Body force (abstract)
    %       g     - Boundary flux eps(u)*n or (grad u)'*n (abstract)
    %       u4Db  - Dirichlet velocity data (abstract)
    %
    %   ProblemStokesTemplate properties:
    %       geom  - Domain geometry (abstract)
    %
    %   See also StokesParameters, ProblemStokesTemplateExact,
    %   problemStokesDrivenCavity, problemStokesBackwardFacingStep

    % Copyright (c) 2026 by the fem4m authors

    properties(Abstract)
        geom
    end

    methods(Abstract)
        val = f(x);
        val = g(x,normals);
        val = u4Db(x);
    end

end