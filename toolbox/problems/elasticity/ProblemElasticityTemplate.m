classdef ProblemElasticityTemplate < ElasticityParameters
    %ProblemElasticityTemplate - Abstract linear elasticity problem
    %   ProblemElasticityTemplate is an abstract template for 2d linear
    %   elasticity problems -div(sigma(u)) = f with sigma(u) = C eps(u);
    %   subclasses supply geometry, scale, and data. It inherits the
    %   material parameters and the constructor of ElasticityParameters,
    %   so a subclass without its own constructor requires the material
    %   parameters NAME1,VALUE1,NAME2,VALUE2 as in ElasticityParameters.
    %
    %   The data functions are evaluated at the points X (N-by-2, one row
    %   per point) and return N-by-2 arrays, one vector per point: f(X)
    %   is the body force, u4Db(X) the displacement on the Dirichlet
    %   boundary, and g(X,NORMALS) the traction sigma*n on the Neumann
    %   boundary for the outer unit normals NORMALS (N-by-2).
    %
    %   ProblemElasticityTemplate functions:
    %       f     - Body force (abstract)
    %       g     - Traction sigma*n on the Neumann boundary (abstract)
    %       u4Db  - Displacement on the Dirichlet boundary (abstract)
    %
    %   ProblemElasticityTemplate properties:
    %       geom   - Domain geometry (abstract)
    %       scale  - Amplification of the displacement plots
    %
    %   See also ElasticityParameters, ProblemElasticityTemplateExact

    % Copyright (c) 2026 by the fem4m authors

    properties(Abstract)
        geom
        scale
    end

    methods(Abstract)
        val = f(x);
        val = g(x,normals);
        val = u4Db(x);
    end

end