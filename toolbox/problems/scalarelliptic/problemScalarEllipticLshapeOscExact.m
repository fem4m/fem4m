classdef problemScalarEllipticLshapeOscExact < ProblemScalarEllipticTemplateExact
    %problemScalarEllipticLshapeOscExact - Oscillating solution on an L-shape
    %   problemScalarEllipticLshapeOscExact is a Poisson problem on the
    %   L-shape (-1,1)^2 \ [0,1]x[-1,0] with the oscillatory exact
    %   solution u = sin(K*pi*x)*sin(K*pi*y) of frequency K.
    %
    %   OBJ = problemScalarEllipticLshapeOscExact() uses K = 1 and
    %   Dirichlet boundary conditions on the whole boundary.
    %
    %   OBJ = problemScalarEllipticLshapeOscExact(K) sets the frequency K.
    %
    %   OBJ = problemScalarEllipticLshapeOscExact(K,BDRY) also sets the
    %   boundary conditions BDRY, one letter 'D' (Dirichlet) or 'N'
    %   (Neumann) per boundary part: part 1 is the two edges at the
    %   re-entrant corner, part 2 the remaining edges (see Geometry).
    %
    %   problemScalarEllipticLshapeOscExact functions:
    %       f           - Source term of the oscillatory solution
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %
    %   problemScalarEllipticLshapeOscExact properties:
    %       geom        - L-shaped domain geometry
    %       k           - Oscillation frequency
    %
    %   See also ProblemScalarEllipticTemplateExact,
    %   problemScalarEllipticLshape, problemScalarEllipticSquareOscExact,
    %   Geometry

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape','bdry','D'))
        k
    end


    methods

        function obj = problemScalarEllipticLshapeOscExact(k,bdry)
            if nargin < 1
                k = 1;
            end
            obj.k = k;
            if nargin > 1 && ~isempty(bdry)
                obj.geom = Geometry(struct('name','Lshape','bdry',bdry));
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