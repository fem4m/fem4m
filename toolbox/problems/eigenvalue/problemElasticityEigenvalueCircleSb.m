classdef problemElasticityEigenvalueCircleSb < ElasticityParameters & ProblemEigenvalue
    %problemElasticityEigenvalueCircleSb - Elastic disk, symmetry boundary
    %   OBJ = problemElasticityEigenvalueCircleSb(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies reference eigenvalues for an elastic unit disk with
    %   symmetry boundary conditions on the whole boundary: the zero
    %   eigenvalue of the rigid rotation (mode 1) for all material
    %   parameters and the double eigenvalue 5.18806 (modes 2 and 3) for
    %   mu = 1 and lambda = 1; other entries are NaN. The material parameters
    %   NAME1,VALUE1,NAME2,VALUE2 are given as in ElasticityParameters and
    %   are required.
    %
    %   This is an academic benchmark, not a physical model, for testing
    %   eigenvalue approximations.
    %
    %   problemElasticityEigenvalueCircleSb functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemElasticityEigenvalueCircleSb properties:
    %       geom        - Disk geometry with symmetry boundary
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemEigenvalue, problemElasticityEigenvalueCircle

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Circle','bdry','S'))
        scale = 0.1
    end

    methods

        function val = eigenvalue(obj,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 0;
            if obj.mu == 1 && obj.lambda == 1
                eigenvalue(2) = 5.18806;
                eigenvalue(3) = 5.18806;
            end
            val = eigenvalue(m);
        end

    end

end