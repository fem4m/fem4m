classdef problemElasticityEigenvalueSquare < ElasticityParameters & ProblemEigenvalue
    %problemElasticityEigenvalueSquare - Elastic square, Dirichlet boundary
    %   OBJ = problemElasticityEigenvalueSquare(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies numerically computed reference eigenvalues for an elastic
    %   unit square with Dirichlet boundary for mu = 1: the modes 1 to 5 for
    %   lambda = 1, 10 and 100, and the first mode for lambda = 1000 and the
    %   Stokes limit lambda = Inf; other entries are NaN. The material
    %   parameters NAME1,VALUE1,NAME2,VALUE2 are given as in
    %   ElasticityParameters and are required.
    %
    %   problemElasticityEigenvalueSquare functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemElasticityEigenvalueSquare properties:
    %       geom        - Square geometry
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemEigenvalue, problemElasticityEigenvalueSquareNeumann,
    %   problemElasticityEigenvalueSquareSb, problemElasticityEigenvalueStripeX

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square'))
        scale = 0.1
    end

    methods

        function val = eigenvalue(obj,m)
            eigenvalue = nan(max(m),1);
            if obj.mu == 1
                if obj.lambda == 1
                    eigenvalue(1:5) = [37.266072199199,37.266072199199,51.294997977322,78.72828218144,99.57599569];
                elseif obj.lambda == 1e1
                    eigenvalue(1:5) = [52.07785993298,80.9236774144,80.923677414,125.130006644,135.594384323];
                elseif obj.lambda == 1e2
                    eigenvalue(1:5) = [52.31315105098,91.477822733,91.477822733,127.906603532,153.04028626];
                elseif obj.lambda == 1e3
                    eigenvalue(1) = 52.341478702;
                elseif isinf(obj.lambda) %Stokes: 52.344691168416544538705330750365
                    eigenvalue(1) = 52.344691168416544538705330750365;
                end
            end
            val = eigenvalue(m);
        end

    end

end