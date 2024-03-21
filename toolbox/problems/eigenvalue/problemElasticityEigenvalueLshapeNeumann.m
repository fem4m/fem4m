classdef problemElasticityEigenvalueLshapeNeumann < ElasticityParameters & ProblemEigenvalue
    %problemElasticityEigenvalueLshapeNeumann - Elastic L-shape, Neumann
    %   OBJ = problemElasticityEigenvalueLshapeNeumann(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies reference eigenvalues for an elastic L-shape
    %   (-1,1)^2 \ [0,1]x[-1,0] with Neumann boundary: the three zero
    %   eigenvalues of the rigid body motions (modes 1 to 3) for all
    %   material parameters and the modes 4 and 5 for mu = 1 and lambda = 1;
    %   other entries are NaN. The material parameters
    %   NAME1,VALUE1,NAME2,VALUE2 are given as in ElasticityParameters and
    %   are required.
    %
    %   problemElasticityEigenvalueLshapeNeumann functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemElasticityEigenvalueLshapeNeumann properties:
    %       geom        - L-shape geometry
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemEigenvalue, problemElasticityEigenvalueLshape

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape','bdry','N'))
        scale = 0.1
    end

    methods

        function val = eigenvalue(obj,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 0;
            eigenvalue(2) = 0;
            eigenvalue(3) = 0;
            if obj.mu == 1 && obj.lambda == 1
                eigenvalue(4) = 0.833454280384;
                eigenvalue(5) = 3.960009701053;
            end
            val = eigenvalue(m);
        end

    end

end