classdef problemElasticityEigenvalueStripeX < ElasticityParameters & ProblemEigenvalue
    %problemElasticityEigenvalueStripeX - Elastic square, symmetry at x = 0
    %   OBJ = problemElasticityEigenvalueStripeX(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies reference eigenvalues for the elastic unit square with a
    %   symmetry boundary on the left edge x = 0 and Neumann boundary on the
    %   other edges. Only the first eigenvalue, the zero eigenvalue of the
    %   translation in y, is available for all material parameters; other
    %   entries are NaN. The material parameters NAME1,VALUE1,NAME2,VALUE2
    %   are given as in ElasticityParameters and are required.
    %
    %   problemElasticityEigenvalueStripeX functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemElasticityEigenvalueStripeX properties:
    %       geom        - Unit square, symmetry boundary at x = 0
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemEigenvalue, problemElasticityEigenvalueStripeY

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry',['N';'N';'N';'S']))
        scale = 0.1
    end

    methods

        function val = eigenvalue(obj,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1) = 0;
            if obj.mu == 1 && obj.lambda == 1
                %eigenvalue(2) = ;
                %eigenvalue(3) = ;
            end
            val = eigenvalue(m);
        end

    end

end