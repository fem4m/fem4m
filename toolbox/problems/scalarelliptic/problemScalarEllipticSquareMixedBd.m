classdef problemScalarEllipticSquareMixedBd < ProblemScalarEllipticTemplate
    %problemScalarEllipticSquareMixedBd - Elliptic problem with mixed boundary
    %   problemScalarEllipticSquareMixedBd is a Poisson problem on the
    %   unit square with f = 0, the homogeneous Dirichlet condition u = 0
    %   on the right edge x = 1 and the Neumann condition grad(u)*n = 1 on
    %   the other three edges.
    %
    %   problemScalarEllipticSquareMixedBd functions:
    %       f     - Source term
    %       u4Db  - Dirichlet data
    %       g     - Neumann data
    %
    %   problemScalarEllipticSquareMixedBd properties:
    %       geom  - Unit square geometry
    %
    %   See also ProblemScalarEllipticTemplate, problemScalarEllipticSquare

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry',['N','D','N','N']))
    end


    methods

        function val = f(~,x)
            val = zeros(size(x,1),1);
        end

        function val = u4Db(~,x)
            val = zeros(size(x,1),1);
        end

        function val = g(~,x,~)
            val = ones(size(x,1),1);
        end


    end

end