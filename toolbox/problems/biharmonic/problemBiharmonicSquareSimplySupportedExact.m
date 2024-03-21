classdef problemBiharmonicSquareSimplySupportedExact < ProblemBiharmonicTemplateExact
    %problemBiharmonicSquareSimplySupportedExact - Simply supported plate
    %   problemBiharmonicSquareSimplySupportedExact defines a simply
    %   supported biharmonic problem on the unit square with an oscillatory
    %   solution sin(M*x*pi)*sin(N*y*pi) of modes M and N.
    %
    %   OBJ = problemBiharmonicSquareSimplySupportedExact() uses M = N = 1.
    %
    %   OBJ = problemBiharmonicSquareSimplySupportedExact(M,N) sets the
    %   modes M and N independently.
    %
    %   problemBiharmonicSquareSimplySupportedExact functions:
    %       f           - Source term of the oscillatory solution
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %       exactD2     - Exact second derivatives
    %
    %   problemBiharmonicSquareSimplySupportedExact properties:
    %       geom        - Unit square geometry
    %       m           - Mode number in x
    %       n           - Mode number in y
    %
    %   See also ProblemBiharmonicTemplateExact,
    %   problemBiharmonicSquareClampedOscExact,
    %   problemBiharmonicSimplySupportedSingularRHSCubeExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        m
        n
    end


    methods

        function obj = problemBiharmonicSquareSimplySupportedExact(m,n)
            if nargin > 0
                obj.m = m;
                obj.n = n;
            else
                obj.m = 1;
                obj.n = 1;
            end
        end


        function val = f(obj,x)
            val = pi^4*(obj.m^2+obj.n^2)^2.*sin(obj.m*pi*x(:,1)).*sin(obj.n*pi*x(:,2));
        end

        function val = exactU(obj,x)
            val = sin(obj.m*pi*x(:,1)).*sin(obj.n*pi*x(:,2));
        end

        function val = exactGradU(obj,x)
            val = [obj.m*pi.*cos(obj.m*pi*x(:,1)).*sin(obj.n*pi*x(:,2)),...
                obj.n*pi.*sin(obj.m*pi*x(:,1)).*cos(obj.n*pi*x(:,2))];
        end

        function val = exactD2(obj,x)
            val = [-(obj.m*pi)^2.*sin(obj.m*pi*x(:,1)).*sin(obj.n*pi*x(:,2)),...
                -(obj.n*pi)^2.*sin(obj.m*pi*x(:,1)).*sin(obj.n*pi*x(:,2)),...
                obj.m*pi*obj.n*pi.*cos(obj.m*pi*x(:,1)).*cos(obj.n*pi*x(:,2))];
        end
    end
end