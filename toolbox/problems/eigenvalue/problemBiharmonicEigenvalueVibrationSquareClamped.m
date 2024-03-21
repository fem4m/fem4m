classdef problemBiharmonicEigenvalueVibrationSquareClamped < ProblemEigenvalue
    %problemBiharmonicEigenvalueVibrationSquareClamped - Vibration on a square
    %   OBJ = problemBiharmonicEigenvalueVibrationSquareClamped() supplies
    %   reference eigenvalues lambda of the clamped-plate vibration problem
    %   Delta^2 u = lambda*u with u = du/dn = 0 on the boundary of the unit
    %   square (0,1)^2: the modes 1, 2, 4 and 11; other entries are NaN.
    %
    %   problemBiharmonicEigenvalueVibrationSquareClamped functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemBiharmonicEigenvalueVibrationSquareClamped properties:
    %       geom        - Square domain geometry
    %
    %   See also ProblemEigenvalue,
    %   problemBiharmonicEigenvalueBucklingSquareClamped,
    %   problemBiharmonicEigenvalueVibrationLshapeClamped,
    %   problemBiharmonicEigenvalueVibrationSlitClamped

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 1294.93397959; %1294.93397959(2048)
            eigenvalue(2) = 5386.65658;
            eigenvalue(4) = 11710.811238;
            eigenvalue(11) = 48414.379725;
            val = eigenvalue(m);
        end

    end

end