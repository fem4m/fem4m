classdef problemElasticityEigenvalueCircleNeumann < ElasticityParameters & ProblemEigenvalue
    %problemElasticityEigenvalueCircleNeumann - Elastic disk, Neumann boundary
    %   OBJ = problemElasticityEigenvalueCircleNeumann(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies reference eigenvalues for an elastic unit disk with Neumann
    %   boundary conditions: only the three zero eigenvalues of the rigid
    %   body motions (modes 1 to 3), for all material parameters; other
    %   entries are NaN. The material parameters NAME1,VALUE1,NAME2,VALUE2
    %   are given as in ElasticityParameters and are required.
    %
    %   problemElasticityEigenvalueCircleNeumann functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemElasticityEigenvalueCircleNeumann properties:
    %       geom        - Disk geometry
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemEigenvalue, problemElasticityEigenvalueCircle

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Circle','bdry','N'))
        scale = 0.1
    end

    methods

        function val = eigenvalue(obj,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 0;
            eigenvalue(2) = 0;
            eigenvalue(3) = 0;
            if obj.mu == 1 && obj.lambda == 1
                %eigenvalue(4) = ;
                %eigenvalue(5) = ;
            end
            val = eigenvalue(m);
        end

    end

end