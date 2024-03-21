classdef problemBiharmonicEigenvalueVibrationSlitClamped < ProblemEigenvalue
    %problemBiharmonicEigenvalueVibrationSlitClamped - Vibration on a slit
    %   OBJ = problemBiharmonicEigenvalueVibrationSlitClamped() supplies
    %   reference eigenvalues lambda of the clamped-plate vibration problem
    %   Delta^2 u = lambda*u with u = du/dn = 0 on the boundary of the slit
    %   domain (-1,1)^2 \ [0,1]x{0}: the modes 1 and 2; other entries are
    %   NaN.
    %
    %   problemBiharmonicEigenvalueVibrationSlitClamped functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemBiharmonicEigenvalueVibrationSlitClamped properties:
    %       geom        - Slit domain geometry
    %
    %   See also ProblemEigenvalue,
    %   problemBiharmonicEigenvalueBucklingSlitClamped,
    %   problemBiharmonicEigenvalueVibrationLshapeClamped,
    %   problemBiharmonicEigenvalueVibrationSquareClamped

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Slit'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 386.809225;
            eigenvalue(2) = 471.682903;
            val = eigenvalue(m);
        end

    end

end