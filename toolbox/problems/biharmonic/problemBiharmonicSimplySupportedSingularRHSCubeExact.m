classdef problemBiharmonicSimplySupportedSingularRHSCubeExact < ProblemBiharmonicTemplateExact
    %problemBiharmonicSimplySupportedSingularRHSCubeExact - Simply supported cube, point load
    %   problemBiharmonicSimplySupportedSingularRHSCubeExact is a
    %   biharmonic problem with simply supported boundary conditions (u
    %   and its second normal derivative g2, or its Laplacian
    %   laplaceU4Db, prescribed) on the cube (-1/2,1/2)^3 with a singular
    %   right-hand side concentrated at the origin and the fundamental
    %   solution u = -|x|/(8*pi) as exact solution. The boundary data
    %   u4Db, g, g2 and laplaceU4Db are taken from the exact solution; g3
    %   is not available (exactDlaplaceU is not implemented).
    %
    %   problemBiharmonicSimplySupportedSingularRHSCubeExact functions:
    %       f             - Point load concentrated at the origin
    %       exactU        - Fundamental solution -|x|/(8*pi) (singular)
    %       exactGradU    - Gradient of the exact solution
    %       exactD2       - Second derivatives of the solution
    %       exactLaplaceU - Laplacian of the exact solution
    %
    %   problemBiharmonicSimplySupportedSingularRHSCubeExact properties:
    %       geom          - Cube geometry centred at the origin
    %
    %   See also ProblemBiharmonicTemplateExact,
    %   problemBiharmonicCubeClampedExact,
    %   problemBiharmonicSquareSimplySupportedExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Cube','midpoint',[0,0,0],'length',1));
    end

    methods

        function val = f(~,x,varargin)
            val = zeros(size(x,1),1);
            val(x(:,1)==0 & x(:,2)==0 & x(:,3)==0) = 1;
        end

        function val = exactU(~,x)
            val = -1./(8*pi).*vecnorm(x,2,2);
        end

        function val = exactGradU(~,x)
            val = -1./(8*pi).*x./vecnorm(x,2,2);
        end

        function val = exactD2(~,x)
            norm = vecnorm(x,2,2);
            y = x(:,2); z = x(:,3); x = x(:,1);
            val = -1./(8*pi).*[ (y.^2 + z.^2),(x.^2 + z.^2),(x.^2 + y.^2), -x.*y, -x.*z, -y.*z ]./norm.^3;
        end

        function val = exactLaplaceU(~,x)
            val = -1./(4*pi)./vecnorm(x,2,2);
        end

    end
end