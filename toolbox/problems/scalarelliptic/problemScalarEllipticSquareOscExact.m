classdef problemScalarEllipticSquareOscExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticSquareOscExact - Oscillating solution on a square
    %   problemScalarEllipticSquareOscExact is a Poisson problem on the
    %   unit square with the oscillatory exact solution
    %   u = sin(K*pi*x)*sin(K*pi*y) of frequency K.
    %
    %   OBJ = problemScalarEllipticSquareOscExact() uses K = 1 and
    %   Dirichlet boundary conditions on the whole boundary.
    %
    %   OBJ = problemScalarEllipticSquareOscExact(K) sets the frequency K.
    %
    %   OBJ = problemScalarEllipticSquareOscExact(K,BDRY) also sets the
    %   boundary conditions BDRY, one letter 'D' (Dirichlet) or 'N'
    %   (Neumann) per edge in the order bottom, right, top, left (see
    %   Geometry).
    %
    %   problemScalarEllipticSquareOscExact functions:
    %       f           - Source term of the oscillatory solution
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticSquareOscExact properties:
    %       geom        - Unit square geometry
    %       k           - Oscillation frequency
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticLshapeOscExact, Geometry

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        k
    end

    methods

        function obj = problemScalarEllipticSquareOscExact(k,bdry)
            if nargin < 1
                k = 1;
            end
            obj.k = k;
            if nargin > 1 && ~isempty(bdry)
                obj.geom = Geometry(struct('name','Square','bdry',bdry));
            end
        end

        function val = f(obj,x)
            val = 2*pi^2*obj.k^2*sin(obj.k*pi*x(:,1)).*sin(obj.k*pi*x(:,2));
        end

        function val = exactU(obj,x)
            val = sin(obj.k*pi*x(:,1)).*sin(obj.k*pi*x(:,2));
        end

        function val = exactGradU(obj,x)
            val = [obj.k*pi*cos(obj.k*pi*x(:,1)).*sin(obj.k*pi*x(:,2)),...
                obj.k*pi*sin(obj.k*pi*x(:,1)).*cos(obj.k*pi*x(:,2))];
        end

    end

end