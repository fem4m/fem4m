classdef problemScalarEllipticEigenvalueLshape < ProblemEigenvalue
    %problemScalarEllipticEigenvalueLshape - L-shape Dirichlet eigenvalues
    %   OBJ = problemScalarEllipticEigenvalueLshape() supplies reference
    %   Dirichlet eigenvalues of the Laplacian on the L-shaped domain
    %   (-1,1)^2 \ [0,1]x[-1,0]: the modes 1 to 10 (mode 1 to double
    %   precision from a 100-digit value, modes 3, 8 and 9 exact, the others
    %   numerical approximations); other entries are NaN.
    %
    %   problemScalarEllipticEigenvalueLshape functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemScalarEllipticEigenvalueLshape properties:
    %       geom        - L-shaped domain geometry
    %
    %   See also ProblemEigenvalue,
    %   problemScalarEllipticEigenvalueLshapeNeumann,
    %   problemScalarEllipticEigenvalueSlit

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max(m),1);
            % First Dirichlet eigenvalue (100-digit reference value, stored as double)
            % cf. Robert Stephen Jones, Computing ultra-precise eigenvalues
            % of the Laplacian within polygons, Adv Comput Math (2017)
            % 43:1325–1354, DOI 10.1007/s10444-017-9527-y
            eigenvalue(1) = 9.639723844021941052711459262364823156267289525821906456109579700564035647863370390722873165008796788;
            % hp-adaptive equilibrated
            eigenvalue(2) = 15.197251926454;%[(4)34-(4)35]
            eigenvalue(3) = 8*pi^2/4;
            eigenvalue(4) = 29.521481114145;%[(4)81-(4)86]
            eigenvalue(5) = 31.912635957138;%[(6)77-(7)94]
            eigenvalue(6) = 41.474509890215;%[(4)89-(4)90]
            eigenvalue(7) = 44.948487781351;%[?    -(1)17];
            eigenvalue(8) = 5*pi^2;
            eigenvalue(9) = 5*pi^2;
            eigenvalue(10)= 56.709609887385;%[?    -(5)02];
            val = eigenvalue(m);
        end

    end

end