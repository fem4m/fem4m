classdef problemBiharmonicEigenvalueBucklingLshapeClamped < ProblemEigenvalue
    %problemBiharmonicEigenvalueBucklingLshapeClamped - Buckling on an L-shape
    %   OBJ = problemBiharmonicEigenvalueBucklingLshapeClamped() supplies
    %   reference eigenvalues lambda of the clamped-plate buckling problem
    %   Delta^2 u = -lambda*Delta u with u = du/dn = 0 on the boundary of the
    %   L-shaped domain (-1,1)^2 \ [0,1]x[-1,0]: the modes 1 and 2; other
    %   entries are NaN.
    %
    %   problemBiharmonicEigenvalueBucklingLshapeClamped functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemBiharmonicEigenvalueBucklingLshapeClamped properties:
    %       geom        - L-shaped domain geometry
    %
    %   See also ProblemEigenvalue,
    %   problemBiharmonicEigenvalueVibrationLshapeClamped,
    %   problemBiharmonicEigenvalueBucklingSlitClamped,
    %   problemBiharmonicEigenvalueBucklingSquareClamped

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 32.13269465;
            eigenvalue(2) = 37.018335;
            val = eigenvalue(m);
        end

    end

end