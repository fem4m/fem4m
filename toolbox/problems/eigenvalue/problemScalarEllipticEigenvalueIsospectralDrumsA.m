classdef problemScalarEllipticEigenvalueIsospectralDrumsA < ProblemEigenvalue
    %problemScalarEllipticEigenvalueIsospectralDrumsA - Isospectral drum, A
    %   OBJ = problemScalarEllipticEigenvalueIsospectralDrumsA() supplies
    %   reference Dirichlet eigenvalues of the Laplacian on the first of a
    %   pair of isospectral polygons with different shapes (Gordon, Webb and
    %   Wolpert): the modes 1 to 3 and 50; other entries are NaN. Both
    %   domains have the same eigenvalues.
    %
    %   problemScalarEllipticEigenvalueIsospectralDrumsA functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemScalarEllipticEigenvalueIsospectralDrumsA properties:
    %       geom        - Isospectral domain A geometry
    %
    %   See also ProblemEigenvalue,
    %   problemScalarEllipticEigenvalueIsospectralDrumsB

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','IsospectralDrumsA'))
    end

    methods

        function val = eigenvalue(~,m)
            eigenvalue = nan(max([m(:);50]),1);
            eigenvalue(1) = 2.5379440;
            eigenvalue(2) = 3.6555097;
            eigenvalue(3) = 5.1755594;
            eigenvalue(50) = 54.187935625408;
            val = eigenvalue(m);
        end

    end

end