classdef ProblemScalarEllipticTemplate
    %ProblemScalarEllipticTemplate - Abstract scalar elliptic problem
    %   ProblemScalarEllipticTemplate is an abstract template for scalar
    %   elliptic problems -div(alpha*grad(u)) + beta*grad(u) + gamma*u = f
    %   with the Dirichlet data u = u4Db and the Neumann data
    %   alpha*grad(u)*n = g; subclasses supply the geometry and the data
    %   and override the coefficients that differ from the defaults.
    %
    %   The data and coefficients are evaluated at the points X (N-by-d,
    %   one row per point). f, u4Db, g and gamma return an N-by-1 vector,
    %   beta an N-by-d array, alpha an N-by-d-by-d array and dalpha an
    %   N-by-d-by-d-by-d array with VAL(:,i,j,k) = d alpha_ij/dx_k. The
    %   defaults of alpha and dalpha return a single row (1-by-d-by-d and,
    %   2d only, 1-by-2-by-2-by-2) that holds at all points.
    %   g(OBJ,X,NORMALS) also takes the outer unit normals NORMALS (N-by-d)
    %   at the points X.
    %
    %   ProblemScalarEllipticTemplate functions:
    %       f       - Source term (abstract)
    %       u4Db    - Dirichlet data (abstract)
    %       g       - Neumann data alpha*grad(u)*n (abstract)
    %       alpha   - Diffusion coefficient (default: identity matrix)
    %       dalpha  - Derivative of the diffusion coefficient (default: 0)
    %       beta    - Convection coefficient (default: zero vector field)
    %       gamma   - Reaction coefficient (default: 0)
    %
    %   ProblemScalarEllipticTemplate properties:
    %       geom    - Domain geometry (abstract)
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticSquare, Geometry

    % Copyright (c) 2026 by the fem4m authors

    properties(Abstract)
        geom
    end


    methods(Abstract)
        val = f(obj,x);
        val = u4Db(obj,x);
        val = g(obj,x,normals);
    end

    methods

        function val = alpha(~,x)
            d = size(x,2);
            val = reshape(eye(d,d),[1 d d]);
        end

        function val = dalpha(~,~)
            val = zeros(1,2,2,2);
        end

        function val = beta(~,x)
            val = zeros(size(x,1),size(x,2));
        end

        function val = gamma(~,x)
            val = zeros(size(x,1),1);
        end

    end
end