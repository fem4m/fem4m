classdef problemElasticityEigenvalueCooks < ElasticityParameters & ProblemEigenvalue
    %problemElasticityEigenvalueCooks - Eigenvalues of the Cook membrane
    %   OBJ = problemElasticityEigenvalueCooks(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies reference eigenvalues of the Cook membrane, the
    %   quadrilateral with the corners (0,-44), (48,0), (48,16) and (0,0),
    %   clamped (Dirichlet) on the left edge x = 0 and traction-free
    %   (Neumann) on the other edges: the modes 1 to 6 for mu = 1 and
    %   lambda = 1; other entries are NaN. The material parameters
    %   NAME1,VALUE1,NAME2,VALUE2 are given as in ElasticityParameters and
    %   are required.
    %
    %   problemElasticityEigenvalueCooks functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemElasticityEigenvalueCooks properties:
    %       geom        - Cook membrane geometry
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemEigenvalue, problemElasticityEigenvalueCircle

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Cooks'))
        scale = 0.1
    end

    methods

        function val = eigenvalue(obj,m)
            eigenvalue = nan(max(m),1);
            if obj.mu == 1 && obj.lambda == 1
                eigenvalue(1:6) = [0.0004029679680725,0.00274619089312,0.0037995216855,0.01113078873485,0.017878285454,0.02229140017];
            end
            val = eigenvalue(m);
        end

    end

end