classdef problemBiharmonicCubeClampedExact < ProblemBiharmonicTemplateExact
    %problemBiharmonicCubeClampedExact - Clamped on unit cube, homogeneous data
    %   problemBiharmonicCubeClampedExact defines a clamped biharmonic
    %   problem on the unit cube with the polynomial exact solution
    %   u = x^2(1-x)^2*y^2(1-y)^2*z^2(1-z)^2, whose value and normal
    %   derivative vanish on the boundary.
    %
    %   problemBiharmonicCubeClampedExact functions:
    %       f           - Source term from the exact solution
    %       g           - Zero normal-derivative data
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %       exactD2     - Exact second derivatives
    %
    %   problemBiharmonicCubeClampedExact properties:
    %       geom        - Unit cube geometry
    %
    %   See also ProblemBiharmonicTemplateExact,
    %   problemBiharmonicCubeClampedExact2

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Cube','midpoint',[0.5,0.5,0.5],'length',1))
    end


    methods

        function val = f(~,pts)
            x = pts(:,1);
            y = pts(:,2);
            z = pts(:,3);
            val = 8.*(3.*x.^4.*y.^4 - 6.*x.^4.*y.^3 + 36.*x.^4.*y.^2.*z.^2 - 36.*x.^4.*y.^2.*z + 9.*x.^4.*y.^2 ...
                - 36.*x.^4.*y.*z.^2 + 36.*x.^4.*y.*z - 6.*x.^4.*y + 3.*x.^4.*z.^4 - 6.*x.^4.*z.^3 + 9.*x.^4.*z.^2 ...
                - 6.*x.^4.*z + x.^4 - 6.*x.^3.*y.^4 + 12.*x.^3.*y.^3 - 72.*x.^3.*y.^2.*z.^2 + 72.*x.^3.*y.^2.*z ...
                - 18.*x.^3.*y.^2 + 72.*x.^3.*y.*z.^2 - 72.*x.^3.*y.*z + 12.*x.^3.*y - 6.*x.^3.*z.^4 + 12.*x.^3.*z.^3 ...
                - 18.*x.^3.*z.^2 + 12.*x.^3.*z - 2.*x.^3 + 36.*x.^2.*y.^4.*z.^2 - 36.*x.^2.*y.^4.*z + 9.*x.^2.*y.^4 ...
                - 72.*x.^2.*y.^3.*z.^2 + 72.*x.^2.*y.^3.*z - 18.*x.^2.*y.^3 + 36.*x.^2.*y.^2.*z.^4 - 72.*x.^2.*y.^2.*z.^3 ...
                + 108.*x.^2.*y.^2.*z.^2 - 72.*x.^2.*y.^2.*z + 15.*x.^2.*y.^2 - 36.*x.^2.*y.*z.^4 + 72.*x.^2.*y.*z.^3 ...
                - 72.*x.^2.*y.*z.^2 + 36.*x.^2.*y.*z - 6.*x.^2.*y + 9.*x.^2.*z.^4 - 18.*x.^2.*z.^3 + 15.*x.^2.*z.^2 ...
                - 6.*x.^2.*z + x.^2 - 36.*x.*y.^4.*z.^2 + 36.*x.*y.^4.*z - 6.*x.*y.^4 + 72.*x.*y.^3.*z.^2 - 72.*x.*y.^3.*z ...
                + 12.*x.*y.^3 - 36.*x.*y.^2.*z.^4 + 72.*x.*y.^2.*z.^3 - 72.*x.*y.^2.*z.^2 + 36.*x.*y.^2.*z - 6.*x.*y.^2 ...
                + 36.*x.*y.*z.^4 - 72.*x.*y.*z.^3 + 36.*x.*y.*z.^2 - 6.*x.*z.^4 + 12.*x.*z.^3 - 6.*x.*z.^2 + 3.*y.^4.*z.^4 ...
                - 6.*y.^4.*z.^3 + 9.*y.^4.*z.^2 - 6.*y.^4.*z + y.^4 - 6.*y.^3.*z.^4 + 12.*y.^3.*z.^3 - 18.*y.^3.*z.^2 ...
                + 12.*y.^3.*z - 2.*y.^3 + 9.*y.^2.*z.^4 - 18.*y.^2.*z.^3 + 15.*y.^2.*z.^2 - 6.*y.^2.*z + y.^2 - 6.*y.*z.^4 ...
                + 12.*y.*z.^3 - 6.*y.*z.^2 + z.^4 - 2.*z.^3 + z.^2);
        end

        function val = g(~,pts,~)
            val = zeros(size(pts,1),1);
        end

        function val = exactU(~,pts)
            x = pts(:,1);
            y = pts(:,2);
            z = pts(:,3);
            val = x.^2.*(1-x).^2.*y.^2.*(1-y).^2.*z.^2.*(1-z).^2;
        end

        function val = exactGradU(~,pts)
            nrPts = size(pts,1);
            x = pts(:,1);
            y = pts(:,2);
            z = pts(:,3);
            val = zeros(nrPts,3);
            val(:,1) = x.*y.^2.*z.^2.*(x-1.0).^2.*(y-1.0).^2.*(z-1.0).^2.*2.0+x.^2.*y.^2.*z.^2.*(x.*2.0-2.0).*(y-1.0).^2.*(z-1.0).^2;
            val(:,2) = x.^2.*y.*z.^2.*(x-1.0).^2.*(y-1.0).^2.*(z-1.0).^2.*2.0+x.^2.*y.^2.*z.^2.*(y.*2.0-2.0).*(x-1.0).^2.*(z-1.0).^2;
            val(:,3) = x.^2.*y.^2.*z.*(x-1.0).^2.*(y-1.0).^2.*(z-1.0).^2.*2.0+x.^2.*y.^2.*z.^2.*(z.*2.0-2.0).*(x-1.0).^2.*(y-1.0).^2;
        end

        function val = exactD2(~,pts)
            nrPts = size(pts,1);
            x = pts(:,1);
            y = pts(:,2);
            z = pts(:,3);
            val = zeros(nrPts,6);
            val(:,1) = x.^2.*y.^2.*z.^2.*(y-1.0).^2.*(z-1.0).^2.*2.0+y.^2.*z.^2.*(x-1.0).^2.*(y-1.0).^2.*(z-1.0).^2.*2.0+x.*y.^2.*z.^2.*(x.*2.0-2.0).*(y-1.0).^2.*(z-1.0).^2.*4.0;
            val(:,2) = x.^2.*y.^2.*z.^2.*(x-1.0).^2.*(z-1.0).^2.*2.0+x.^2.*z.^2.*(x-1.0).^2.*(y-1.0).^2.*(z-1.0).^2.*2.0+x.^2.*y.*z.^2.*(y.*2.0-2.0).*(x-1.0).^2.*(z-1.0).^2.*4.0;
            val(:,3) = x.^2.*y.^2.*z.^2.*(x-1.0).^2.*(y-1.0).^2.*2.0+x.^2.*y.^2.*(x-1.0).^2.*(y-1.0).^2.*(z-1.0).^2.*2.0+x.^2.*y.^2.*z.*(z.*2.0-2.0).*(x-1.0).^2.*(y-1.0).^2.*4.0;
            val(:,4) = x.*y.^2.*z.^2.*(y.*2.0-2.0).*(x-1.0).^2.*(z-1.0).^2.*2.0+x.^2.*y.*z.^2.*(x.*2.0-2.0).*(y-1.0).^2.*(z-1.0).^2.*2.0+x.^2.*y.^2.*z.^2.*(x.*2.0-2.0).*(y.*2.0-2.0).*(z-1.0).^2+x.*y.*z.^2.*(x-1.0).^2.*(y-1.0).^2.*(z-1.0).^2.*4.0;
            val(:,5) = x.*y.^2.*z.^2.*(z.*2.0-2.0).*(x-1.0).^2.*(y-1.0).^2.*2.0+x.^2.*y.^2.*z.*(x.*2.0-2.0).*(y-1.0).^2.*(z-1.0).^2.*2.0+x.^2.*y.^2.*z.^2.*(x.*2.0-2.0).*(z.*2.0-2.0).*(y-1.0).^2+x.*y.^2.*z.*(x-1.0).^2.*(y-1.0).^2.*(z-1.0).^2.*4.0;
            val(:,6) = x.^2.*y.*z.^2.*(z.*2.0-2.0).*(x-1.0).^2.*(y-1.0).^2.*2.0+x.^2.*y.^2.*z.*(y.*2.0-2.0).*(x-1.0).^2.*(z-1.0).^2.*2.0+x.^2.*y.^2.*z.^2.*(y.*2.0-2.0).*(z.*2.0-2.0).*(x-1.0).^2+x.^2.*y.*z.*(x-1.0).^2.*(y-1.0).^2.*(z-1.0).^2.*4.0;
        end

    end

end