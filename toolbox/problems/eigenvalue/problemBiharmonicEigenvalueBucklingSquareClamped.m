classdef problemBiharmonicEigenvalueBucklingSquareClamped < ProblemEigenvalue
    %problemBiharmonicEigenvalueBucklingSquareClamped - Buckling on a square
    %   OBJ = problemBiharmonicEigenvalueBucklingSquareClamped() supplies
    %   reference eigenvalues lambda of the clamped-plate buckling problem
    %   Delta^2 u = -lambda*Delta u with u = du/dn = 0 on the boundary of the
    %   unit square (0,1)^2: the modes 1 to 12; other entries are NaN.
    %
    %   problemBiharmonicEigenvalueBucklingSquareClamped functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemBiharmonicEigenvalueBucklingSquareClamped properties:
    %       geom        - Square domain geometry
    %
    %   See also ProblemEigenvalue,
    %   problemBiharmonicEigenvalueVibrationSquareClamped,
    %   problemBiharmonicEigenvalueBucklingLshapeClamped,
    %   problemBiharmonicEigenvalueBucklingSlitClamped

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 52.344691168416544538705330750365; %52.344691168416544538705330750365
            eigenvalue(2) =   92.124393972; %92.12439397187146
            eigenvalue(3) =   92.124393972; %92.12439397191328
            eigenvalue(4) =  128.2095843136; %128.209584313673
            eigenvalue(5) =  154.125463071; %154.1254630710968
            eigenvalue(6) =  167.029175273; %167.0291752737527
            eigenvalue(7) =  189.571868112; %189.5718681125406
            eigenvalue(8) =  189.571868112; %189.5718681126013
            eigenvalue(9) =  246.32226976; %246.3222697669958
            eigenvalue(10) = 246.32226976; %246.3222697670244
            eigenvalue(11) = 246.32719675; %246.327196757048
            eigenvalue(12)=  269.12098008; %269.1209880083672
            val = eigenvalue(m);
        end

    end

end