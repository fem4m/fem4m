classdef ProblemStokesTemplateExact < ProblemStokesTemplate
    %ProblemStokesTemplateExact - Stokes problem with a known solution
    %   ProblemStokesTemplateExact extends ProblemStokesTemplate with an
    %   exact velocity and pressure; subclasses provide them, and the
    %   Dirichlet data are taken from the exact velocity.
    %
    %   The exact solution is evaluated at the points X (N-by-d, one row
    %   per point): exactU returns an N-by-d array, exactP and exactDivU
    %   return N-by-1 vectors, and exactGradU returns an N-by-d-by-d array
    %   with VAL(:,i,j) = d_j u_i.
    %
    %   ProblemStokesTemplateExact functions:
    %       exactP     - Exact pressure (abstract)
    %       exactU     - Exact velocity (abstract)
    %       exactGradU - Exact velocity gradient (abstract)
    %       u4Db       - Dirichlet velocity, taken from exactU
    %       g          - (grad u)'*n from exactGradU, without NU and
    %                    pressure, as an N-by-1-by-d array
    %       exactDivU  - Divergence of the exact velocity
    %
    %   See also ProblemStokesTemplate, problemStokesExactSlit,
    %   problemStokesExactSquareSinus, problemStokesSquareMonomExact

    % Copyright (c) 2026 by the fem4m authors

    methods(Abstract)
        val = exactP(x);
        val = exactU(x);
        val = exactGradU(x);
    end

    methods
        function val = u4Db(obj,x)
            val = obj.exactU(x);
        end

        function val = g(obj,x,normals)
            gVal = obj.exactGradU(x);
            val = sum(gVal.*normals,2);
        end

        function val = exactDivU(obj,x)
            grad = obj.exactGradU(x);
            val = zeros(size(grad,1),1);
            for k = 1 : size(grad,2)
                val = val + grad(:,k,k);
            end
        end

    end

end