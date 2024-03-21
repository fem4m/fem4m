classdef problemElasticityEigenvalueCircle < ElasticityParameters & ProblemEigenvalue
    %problemElasticityEigenvalueCircle - Elastic disk, Dirichlet boundary
    %   OBJ = problemElasticityEigenvalueCircle(NAME1,VALUE1,NAME2,VALUE2)
    %   supplies reference eigenvalues for an elastic unit disk with
    %   Dirichlet boundary conditions: the modes 1 to 10 for mu = 1 and
    %   lambda = 1; other entries are NaN. The material parameters
    %   NAME1,VALUE1,NAME2,VALUE2 are given as in ElasticityParameters and
    %   are required.
    %
    %   The eigenvalues w^2 of -div(C eps(u)) = w^2 u on the unit disk with
    %   u = 0 on the boundary are, for the Fourier modes n = 0,1,..., the
    %   roots of kp*ks*J_n'(kp)*J_n'(ks) - n^2*J_n(kp)*J_n(ks) = 0 with
    %   kp = w/sqrt(lambda+2*mu) and ks = w/sqrt(mu), double for n > 0.
    %
    %   problemElasticityEigenvalueCircle functions:
    %       eigenvalue  - Reference eigenvalue of mode M
    %
    %   problemElasticityEigenvalueCircle properties:
    %       geom        - Disk geometry
    %       scale       - Amplification of the displacement plots
    %
    %   See also ProblemEigenvalue, problemElasticityEigenvalueCircleNeumann,
    %   problemElasticityEigenvalueCircleSb

    % Copyright (c) 2026 by the fem4m authors

    properties
        geom = Geometry(struct('name','Circle','bdry','D'))
        scale = 0.1
    end

    methods

        function val = eigenvalue(obj,m)
            eigenvalue = nan(max(m),1);
            if obj.mu == 1 && obj.lambda == 1
                eigenvalue(1:2) = 11.3221446428006;   % n = 1
                eigenvalue(3) = 14.6819706421239;     % n = 0, ks = j_{1,1}
                eigenvalue(4:5) = 27.2731008528181;   % n = 2
                eigenvalue(6:7) = 28.9350715078811;   % n = 1
                eigenvalue(8) = 44.0459119263717;     % n = 0, kp = j_{1,1}
                eigenvalue(9:10) = 45.7626072410073;  % n = 3
            end
            val = eigenvalue(m);
        end

    end

end