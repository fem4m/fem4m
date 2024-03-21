classdef problemElasticityEigenvalueSquareSb < ElasticityParameters & ProblemEigenvalue
    %problemElasticityEigenvalueSquareSb - Elastic square, symmetry boundary
    %   OBJ = problemElasticityEigenvalueSquareSb(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies the analytic eigenvalues mu*pi^2*(k^2+l^2) with k,l >= 1
    %   and (lambda+2*mu)*pi^2*(k^2+l^2) with k,l >= 0 not both zero, sorted
    %   in ascending order with multiplicity, of an elastic unit square with
    %   symmetry boundary conditions on all edges, for all material
    %   parameters. The material parameters NAME1,VALUE1,NAME2,VALUE2 are
    %   given as in ElasticityParameters and are required.
    %
    %   problemElasticityEigenvalueSquareSb functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemElasticityEigenvalueSquareSb properties:
    %       geom        - Square geometry with symmetry boundary
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemEigenvalue, problemElasticityEigenvalueSquare

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Square','bdry','S'))
        scale = 0.1
    end

    methods

        function val = eigenvalue(obj,m)
            mMax = max(m(:));
            M = 1:(mMax+1);
            L = (1:(mMax+1))';
            smodes = obj.mu*pi^2*(M.^2 + L.^2);
            M = 0:mMax;
            L = (0:mMax)';
            pmodes = (obj.lambda+2*obj.mu)*pi^2*(M.^2 + L.^2);
            pmodes(1) = [];
            eigenvalue = sort([smodes(:);pmodes(:)]);
            val = eigenvalue(m);
        end

    end

end