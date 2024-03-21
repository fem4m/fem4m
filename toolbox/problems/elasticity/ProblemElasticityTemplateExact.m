classdef ProblemElasticityTemplateExact < ProblemElasticityTemplate
    %ProblemElasticityTemplateExact - Elasticity with a known exact solution
    %   ProblemElasticityTemplateExact extends ProblemElasticityTemplate
    %   with abstract methods for an exact solution; subclasses provide it.
    %
    %   The exact solution is evaluated at the points X (N-by-2, one row
    %   per point). exactU returns the displacement (N-by-2), exactGradU
    %   the gradient (N-by-2-by-2 with (i,j) = d_j u_i), exactSigma the
    %   stress [s11 s22 s12] and exactEpsU the strain [e11 e22 e12]
    %   (N-by-3 each), exactD2U the second derivatives (N-by-2-by-3 with
    %   (i,:) = [d_xx d_yy d_xy] u_i). The Dirichlet data gradu4Dbt,
    %   gradu4Dbtt and d2u4Dbntt also take the unit tangents TANGENTS
    %   (N-by-2); d2u4Dbntt uses the normal n = (t2,-t1). g(X,NORMALS)
    %   returns sigma*n for the outer unit normals NORMALS (N-by-2).
    %
    %   ProblemElasticityTemplateExact functions:
    %       u4Db          - Displacement, taken from exactU
    %       gradu4Dbtt    - Tangential strain t'*grad(u)*t (N-by-1)
    %       gradu4Dbt     - Tangential derivative grad(u)*t (N-by-2)
    %       d2u4Dbntt     - Second tangential derivative n'*(d_tt u) of the
    %                       displacement (N-by-1), needs exactD2U
    %       g             - Traction sigma*n, from exactSigma
    %       exactD2U      - Second derivatives of the solution
    %                       (default: zero)
    %       exactEpsU     - Symmetric strain of the solution, from
    %                       exactGradU
    %       exactSkewU    - Skew part s = (d_2 u_1-d_1 u_2)/2 of the
    %                       gradient, returned as [0 0 s -s] (N-by-4)
    %       exactU        - Exact displacement (abstract)
    %       exactSigma    - Exact stress (abstract)
    %       exactGradU    - Exact displacement gradient (abstract)
    %
    %   See also ProblemElasticityTemplate, problemElasticitySquareExact

    % Copyright (c) 2026 by the fem4m authors

    methods(Abstract)
        val = exactU(x);
        val = exactSigma(x);
        val = exactGradU(x);
    end

    methods
        function val = u4Db(obj,x)
            val = obj.exactU(x);
        end

        function val = gradu4Dbtt(obj,x,tangents)
            val = obj.exactGradU(x);
            val = sum([val(:,1,1),val(:,2,2),val(:,1,2)+val(:,2,1)].*[tangents.^2,prod(tangents,2)],2);
        end

        function val = gradu4Dbt(obj,x,tangents)
            val = obj.exactGradU(x);
            val = [val(:,1,1).*tangents(:,1)+val(:,1,2).*tangents(:,2),...
                val(:,2,1).*tangents(:,1)+val(:,2,2).*tangents(:,2)];
        end

        function val = d2u4Dbntt(obj,x,tangents)
            n1 = tangents(:,2); n2 = -tangents(:,1);
            u = obj.exactD2U(x);
            val = n2.^3.*u(:,2,1)+n1.^3.*u(:,1,2)...
                +n1.*n2.^2.*u(:,1,1)-n1.^2.*n2.*u(:,1,3).*2.0...
                -n1.*n2.^2.*u(:,2,3).*2.0+n1.^2.*n2.*u(:,2,2);
        end

        function val = g(obj,x,normals)
            sigma = obj.exactSigma(x);
            val = [sigma(:,1).*normals(:,1)+sigma(:,3).*normals(:,2),...
                sigma(:,3).*normals(:,1)+sigma(:,2).*normals(:,2)];
        end

        function val = exactD2U(~,x)
            val = zeros(size(x,1),2,3);
        end

        function val = exactEpsU(obj,x)
            val = obj.exactGradU(x);
            val = [val(:,1,1),val(:,2,2),1/2*(val(:,1,2)+val(:,2,1))];
        end

        function val = exactSkewU(obj,x)
            grad = obj.exactGradU(x);
            val = zeros(size(x,1),4);
            val(:,3) = 1/2*(grad(:,1,2)-grad(:,2,1));
            val(:,4) = - val(:,3);
        end

    end

end