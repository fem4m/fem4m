classdef ProblemBiharmonicTemplate
    %ProblemBiharmonicTemplate - Abstract biharmonic plate problem
    %   ProblemBiharmonicTemplate is an abstract template for biharmonic
    %   (plate) problems such as Delta^2 u = f with clamped, simply
    %   supported or free boundary conditions; subclasses supply the
    %   geometry and data.
    %
    %   The data are evaluated at the points X (N-by-d, one row per
    %   point) and return an N-by-1 vector (N-by-3 in the vector-valued
    %   problemBendingSquare); g, g2 and g3 also take the outer unit
    %   normals NORMALS (N-by-d) at these points.
    %
    %   ProblemBiharmonicTemplate functions:
    %       f             - Source term (abstract)
    %       g             - Normal derivative grad(u)*n on the boundary
    %                       (abstract)
    %       u4Db          - Dirichlet data (abstract)
    %       laplaceU4Db   - Laplacian of the solution on the boundary
    %                       (default: 0)
    %       g2            - Second normal derivative n'*D2u*n on the
    %                       boundary (default: 0)
    %       g3            - Normal derivative grad(Delta u)*n on the
    %                       boundary (default: 0)
    %
    %   ProblemBiharmonicTemplate properties:
    %       geom          - Domain geometry (abstract)
    %
    %   See also ProblemBiharmonicTemplateExact, problemBendingSquare,
    %   problemBiharmonicLshapeClamped

    % Copyright (c) 2026 by the fem4m authors

    properties(Abstract)
        geom
    end


    methods(Abstract)
        val = f(x);
        val = g(x,normals);
        val = u4Db(x);
    end

    methods

        function val = laplaceU4Db(~,x)
            val = zeros(size(x,1),1);
        end

        function val = g2(~,x,~)
            val = zeros(size(x,1),1);
        end

        function val = g3(~,x,~)
            val = zeros(size(x,1),1);
        end

    end

end