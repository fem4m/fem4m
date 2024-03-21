classdef ProblemScalarEllipticTemplateExact < ProblemScalarEllipticTemplate
    %ProblemScalarEllipticTemplateExact - Scalar elliptic problem with exact solution
    %   ProblemScalarEllipticTemplateExact extends
    %   ProblemScalarEllipticTemplate with a known exact solution, from
    %   which it derives the Dirichlet and Neumann data; subclasses
    %   implement f, exactU, exactGradU and geom.
    %
    %   The functions are evaluated at the points X (N-by-d, one row per
    %   point). exactU returns an N-by-1 vector, exactGradU and exactSigma
    %   an N-by-d array, and g(OBJ,X,NORMALS), with the outer unit normals
    %   NORMALS (N-by-d), an N-by-1 vector. For a vector-valued solution
    %   (problemVectorEllipticSquareExact) f, u4Db and exactU return
    %   N-by-d arrays and exactGradU an N-by-d-by-d array with
    %   VAL(:,i,j) = d u_i/dx_j; exactSigma then is the flux grad(u) of
    %   the components only for the default alpha = I, and the inherited
    %   g assumes a scalar solution.
    %
    %   ProblemScalarEllipticTemplateExact functions:
    %       exactU      - Exact solution (abstract)
    %       exactGradU  - Exact gradient of the solution (abstract)
    %       u4Db        - Dirichlet data, taken from exactU
    %       g           - Neumann data alpha*grad(u)*n, from exactGradU
    %       exactSigma  - Diffusive flux alpha*grad(u) of the solution
    %
    %   See also ProblemScalarEllipticTemplate,
    %   problemScalarEllipticSquareExact

    % Copyright (c) 2026 by the fem4m authors

    methods(Abstract)
        val = exactU(obj,x);
        val = exactGradU(obj,x);
    end

    methods

        function val = u4Db(obj,x)
            val = obj.exactU(x);
        end

        function val = g(obj,x,normals4Nb)
            val = sum(obj.exactSigma(x).*normals4Nb,2);
        end

        function val = exactSigma(obj,x)
            val = matMul(obj.alpha(x),obj.exactGradU(x));
        end

    end
end