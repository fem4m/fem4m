classdef ProblemEigenvalue
    %ProblemEigenvalue - Abstract base class for eigenvalue problems
    %   ProblemEigenvalue is an abstract template for eigenvalue problems;
    %   subclasses supply the geometry and reference eigenvalues.
    %
    %   The reference eigenvalues are returned by eigenvalue(OBJ,M) for the
    %   modes M (positive integers), one value per entry of M. The modes are
    %   numbered by increasing eigenvalue, and a multiple eigenvalue is
    %   repeated according to its multiplicity. The value is NaN for a mode
    %   without a reference value; each subclass names its known modes.
    %
    %   ProblemEigenvalue functions:
    %       eigenvalue  - Reference eigenvalues of the modes M (abstract)
    %
    %   ProblemEigenvalue properties:
    %       geom        - Domain geometry
    %
    %   See also problemBiharmonicEigenvalueBucklingLshapeClamped,
    %   problemScalarEllipticEigenvalueSquare,
    %   problemElasticityEigenvalueSquare, problemSurfaceEigenvalueSphere

    % Copyright (c) 2026 by the fem4m authors

    properties(Abstract)
        geom
    end

    methods(Abstract)
        val = eigenvalue(obj,m);
    end

end