classdef problemElasticityEigenvalueSquareNeumann < ElasticityParameters & ProblemEigenvalue
    %problemElasticityEigenvalueSquareNeumann - Elastic square, Neumann
    %   OBJ = problemElasticityEigenvalueSquareNeumann(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies reference eigenvalues of an elastic unit square with
    %   Neumann boundary: the three zero eigenvalues of the rigid body
    %   motions (modes 1 to 3) for all material parameters and, for mu = 1
    %   and lambda = 1, 10 or lambda >= 100, the analytic eigenvalues
    %   2*k^2*pi^2, k = 1,...,9, at their mode numbers, which depend on
    %   lambda; other entries are NaN. The material parameters
    %   NAME1,VALUE1,NAME2,VALUE2 are given as in ElasticityParameters and
    %   are required.
    %
    %   problemElasticityEigenvalueSquareNeumann functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemElasticityEigenvalueSquareNeumann properties:
    %       geom        - Square geometry
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemEigenvalue, problemElasticityEigenvalueSquare

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry','N'))
        scale = 0.1
    end

    methods

        function val = eigenvalue(obj,m)
            eigenvalue = nan(max(m),1);
            eigenvalue(1:3) = 0;
            if obj.mu == 1
                if obj.lambda == 1
                    eigenvalue([7,16,29,49,69,95,123,161,197]) = 2*(1:9).^2*pi^2;
                elseif obj.lambda == 1e1
                    eigenvalue([5,13,25,42,57,81,103,131,163]) = 2*(1:9).^2*pi^2;
                elseif obj.lambda >= 1e2
                    eigenvalue([5,13,25,39,56,75,97,125,155]) = 2*(1:9).^2*pi^2;
                end
            end
            val = eigenvalue(m);
        end

    end

end