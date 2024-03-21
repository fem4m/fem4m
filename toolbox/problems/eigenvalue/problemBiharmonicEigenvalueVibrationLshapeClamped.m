classdef problemBiharmonicEigenvalueVibrationLshapeClamped < ProblemEigenvalue
    %problemBiharmonicEigenvalueVibrationLshapeClamped - Vibration on an L-shape
    %   OBJ = problemBiharmonicEigenvalueVibrationLshapeClamped() supplies
    %   reference eigenvalues lambda of the clamped-plate vibration problem
    %   Delta^2 u = lambda*u with u = du/dn = 0 on the boundary of the
    %   L-shaped domain (-1,1)^2 \ [0,1]x[-1,0]: the modes 1 and 2; other
    %   entries are NaN.
    %
    %   problemBiharmonicEigenvalueVibrationLshapeClamped functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemBiharmonicEigenvalueVibrationLshapeClamped properties:
    %       geom        - L-shaped domain geometry
    %
    %   See also ProblemEigenvalue,
    %   problemBiharmonicEigenvalueBucklingLshapeClamped,
    %   problemBiharmonicEigenvalueVibrationSlitClamped,
    %   problemBiharmonicEigenvalueVibrationSquareClamped

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 418.97529285; %418.9752928-474/508/513->average 500
            eigenvalue(2) = 690.90651;
            val = eigenvalue(m);
        end

    end

end