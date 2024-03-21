classdef problemElasticityEigenvalueLshapeSb < ElasticityParameters & ProblemEigenvalue
    %problemElasticityEigenvalueLshapeSb - Elastic L-shape, symmetry boundary
    %   OBJ = problemElasticityEigenvalueLshapeSb(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies reference eigenvalues for an elastic L-shape
    %   (-1,1)^2 \ [0,1]x[-1,0] with symmetry boundary conditions on all
    %   edges: the modes 1 and 2 for mu = 1 and lambda = 1; other entries
    %   are NaN. The material parameters NAME1,VALUE1,NAME2,VALUE2 are
    %   given as in ElasticityParameters and are required.
    %
    %   This is an academic benchmark, not a physical model, for testing
    %   eigenvalue approximations.
    %
    %   problemElasticityEigenvalueLshapeSb functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemElasticityEigenvalueLshapeSb properties:
    %       geom        - L-shape geometry with symmetry boundary
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemEigenvalue, problemElasticityEigenvalueLshape

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape','bdry','S'))
        scale = 0.1
    end

    methods

        function val = eigenvalue(obj,m)
            eigenvalue = nan(max(m),1);
            if obj.mu == 1 && obj.lambda == 1
                eigenvalue(1) = 6.4197092418;
                eigenvalue(2) = 10.6020941;
                %eigenvalue(3) = ;
            end
            val = eigenvalue(m);
        end

    end

end