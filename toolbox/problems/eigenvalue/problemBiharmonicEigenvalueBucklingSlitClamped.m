classdef problemBiharmonicEigenvalueBucklingSlitClamped < ProblemEigenvalue
    %problemBiharmonicEigenvalueBucklingSlitClamped - Buckling on a slit
    %   OBJ = problemBiharmonicEigenvalueBucklingSlitClamped() supplies
    %   reference eigenvalues lambda of the clamped-plate buckling problem
    %   Delta^2 u = -lambda*Delta u with u = du/dn = 0 on the boundary of the
    %   slit domain (-1,1)^2 \ [0,1]x{0}: the modes 1 and 2; other entries
    %   are NaN.
    %
    %   problemBiharmonicEigenvalueBucklingSlitClamped functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemBiharmonicEigenvalueBucklingSlitClamped properties:
    %       geom        - Slit domain geometry
    %
    %   See also ProblemEigenvalue,
    %   problemBiharmonicEigenvalueVibrationSlitClamped,
    %   problemBiharmonicEigenvalueBucklingLshapeClamped,
    %   problemBiharmonicEigenvalueBucklingSquareClamped

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Slit'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 29.9168629;
            eigenvalue(2) = 32.589993;
            val = eigenvalue(m);
        end

    end

end