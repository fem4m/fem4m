classdef problemBiharmonicSquareClampedOscExact < ProblemBiharmonicTemplateExact
    %problemBiharmonicSquareClampedOscExact - Clamped plate, oscillatory
    %   problemBiharmonicSquareClampedOscExact defines a clamped biharmonic
    %   problem on the unit square with an oscillatory solution
    %   sin(M*x*pi)^2*sin(M*y*pi)^2 of mode M; for integer M its value
    %   and normal derivative vanish on the boundary.
    %
    %   OBJ = problemBiharmonicSquareClampedOscExact() uses mode M = 1.
    %
    %   OBJ = problemBiharmonicSquareClampedOscExact(M) sets the mode M
    %   (an empty M uses the default 1).
    %
    %   problemBiharmonicSquareClampedOscExact functions:
    %       f           - Source term of the oscillatory solution
    %       exactU      - Exact solution
    %       exactGradU  - Exact gradient of the solution
    %       exactD2     - Exact second derivatives
    %
    %   problemBiharmonicSquareClampedOscExact properties:
    %       geom        - Unit square geometry
    %       m           - Mode number of the solution
    %
    %   See also ProblemBiharmonicTemplateExact,
    %   problemBiharmonicSquareClampedExact,
    %   problemBiharmonicSquareSimplySupportedExact

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        m
    end


    methods

        function obj = problemBiharmonicSquareClampedOscExact(m)
            if nargin > 0 && ~isempty(m)
                obj.m = m;
            else
                obj.m = 1;
            end
        end

        function val = f(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.m;
            val = k.^4.*pi.^4.*cos(k.*x.*pi).^2.*cos(k.*y.*pi).^2.*8.0-k.^4.*pi.^4.*sin(k.*x.*pi).^2.*cos(k.*y.*pi).^2.*1.6e+1...
                -k.^4.*pi.^4.*sin(k.*y.*pi).^2.*cos(k.*x.*pi).^2.*1.6e+1+k.^4.*pi.^4.*sin(k.*x.*pi).^2.*sin(k.*y.*pi).^2.*2.4e+1;
        end

        function val = exactU(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.m;
            val = sin(k.*x.*pi).^2.*sin(k.*y.*pi).^2;
        end

        function val = exactGradU(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.m;
            val = [k.*pi.*sin(k.*x.*pi).*sin(k.*y.*pi).^2.*cos(k.*x.*pi).*2.0,k.*pi.*sin(k.*x.*pi).^2.*sin(k.*y.*pi).*cos(k.*y.*pi).*2.0];
        end

        function val = exactD2(obj,x)
            y = x(:,2);
            x = x(:,1);
            k = obj.m;
            val = [k.^2.*pi.^2.*sin(k.*y.*pi).^2.*cos(k.*x.*pi).^2.*2.0-k.^2.*pi.^2.*sin(k.*x.*pi).^2.*sin(k.*y.*pi).^2.*2.0,...
                k.^2.*pi.^2.*sin(k.*x.*pi).^2.*cos(k.*y.*pi).^2.*2.0-k.^2.*pi.^2.*sin(k.*x.*pi).^2.*sin(k.*y.*pi).^2.*2.0,...
                k.^2.*pi.^2.*sin(k.*x.*pi).*sin(k.*y.*pi).*cos(k.*x.*pi).*cos(k.*y.*pi).*4.0];
        end

    end

end