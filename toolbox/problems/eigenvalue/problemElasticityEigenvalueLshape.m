classdef problemElasticityEigenvalueLshape < ElasticityParameters & ProblemEigenvalue
    %problemElasticityEigenvalueLshape - Elastic L-shape, Dirichlet boundary
    %   OBJ = problemElasticityEigenvalueLshape(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies reference eigenvalues for an elastic L-shape
    %   (-1,1)^2 \ [0,1]x[-1,0] with Dirichlet boundary: the modes 1 to 5
    %   for mu = 1 and lambda = 1; other entries are NaN. The material
    %   parameters NAME1,VALUE1,NAME2,VALUE2 are given as in
    %   ElasticityParameters and are required.
    %
    %   problemElasticityEigenvalueLshape functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemElasticityEigenvalueLshape properties:
    %       geom        - L-shape geometry
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemEigenvalue, problemElasticityEigenvalueLshapeNeumann,
    %   problemElasticityEigenvalueLshapeSb

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Lshape'))
        scale = 0.1
    end

    methods

        function val = eigenvalue(obj,m)
            eigenvalue = nan(max(m),1);
            if obj.mu == 1 && obj.lambda == 1
                eigenvalue(1:5) = [13.591920213185,17.2708242904,26.5296796275,32.2965878967,33.7303348734];
            end
            val = eigenvalue(m);
        end

    end

end