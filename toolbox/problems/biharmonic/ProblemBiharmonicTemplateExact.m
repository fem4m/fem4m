classdef ProblemBiharmonicTemplateExact < ProblemBiharmonicTemplate
    %ProblemBiharmonicTemplateExact - Biharmonic problem with known solution
    %   ProblemBiharmonicTemplateExact extends ProblemBiharmonicTemplate
    %   with abstract methods for an exact solution; subclasses provide
    %   it, and the boundary data are computed from it.
    %
    %   The exact solution is evaluated at the points X (N-by-d, one row
    %   per point). exactU and exactLaplaceU return an N-by-1 vector,
    %   exactGradU and exactDlaplaceU an N-by-d array, and exactD2 the
    %   second derivatives as an N-by-3 array [xx yy xy] in 2d or an
    %   N-by-6 array [xx yy zz xy xz yz] in 3d. exactHessianU returns the
    %   N-by-d-by-d array of the Hessians.
    %
    %   ProblemBiharmonicTemplateExact functions:
    %       u4Db           - Dirichlet data, taken from exactU
    %       g              - Normal derivative grad(u)*n, from exactGradU
    %       gradUt         - Tangential derivative grad(u)*t, from
    %                        exactGradU
    %       g2             - Second normal derivative n'*D2u*n, from
    %                        exactD2
    %       g3             - Normal derivative grad(Delta u)*n, from
    %                        exactDlaplaceU
    %       D2Ut           - Hessian times tangent D2u*t (N-by-d), from
    %                        exactHessianU
    %       laplaceU4Db    - Laplacian on the boundary, from exactLaplaceU
    %       exactU         - Exact solution (abstract)
    %       exactGradU     - Exact gradient of the solution (abstract)
    %       exactD2        - Exact second derivatives (abstract)
    %       exactLaplaceU  - Laplacian of the exact solution, from exactD2
    %       exactDlaplaceU - Gradient of the Laplacian (subclasses must
    %                        override it to use g3)
    %       exactHessianU  - Hessian of the exact solution, from exactD2
    %
    %   See also ProblemBiharmonicTemplate, problemBiharmonicCubeClampedExact

    % Copyright (c) 2026 by the fem4m authors

    methods(Abstract)
        val = exactU(x);
        val = exactGradU(x);
        val = exactD2(x);
    end


    methods

        function val = u4Db(obj,x)
            val = obj.exactU(x);
        end

        function val = g(obj,x,normals)
            val = sum(obj.exactGradU(x).*normals,2);
        end

        function val = gradUt(obj,x,tangents)
            val = sum(obj.exactGradU(x).*tangents,2);
        end

        function val = g2(obj,x,normals)
            d2Val = obj.exactD2(x);
            if size(x,2)==2
                normals = [normals.^2,2*prod(normals,2)];
            elseif size(x,2)==3
                normals = [normals.^2,2*prod(normals(:,[1 2]),2),2*prod(normals(:,[1 3]),2),2*prod(normals(:,[2 3]),2)];
            end
            val = sum(d2Val.*normals,2);
        end

        function val = g3(obj,x,normals)
            val = sum(obj.exactDlaplaceU(x).*normals,2);
        end

        function val = D2Ut(obj,x,tangents)
            val = matMul(obj.exactHessianU(x),tangents);
        end

        function val = laplaceU4Db(obj,x)
            val = obj.exactLaplaceU(x);
        end

        function val = exactLaplaceU(obj,x)
            d2 = obj.exactD2(x);
            val = sum(d2(:,1:size(x,2)),2);
        end

        function val = exactDlaplaceU(obj,~) %#ok<STOUT>
            error('fem4m:ProblemBiharmonicTemplateExact:notImplemented','%s: exactDlaplaceU (gradient of the Laplacian) not implemented, required by g3',class(obj));
        end

        function val = exactHessianU(obj,x)
            d2 = obj.exactD2(x);
            if size(x,2)==2
                val = zeros(size(d2,1),2,2);
                val(:,1,1) = d2(:,1);
                val(:,2,2) = d2(:,2);
                val(:,1,2) = d2(:,3);
                val(:,2,1) = d2(:,3);
            elseif size(x,2)==3
                val = zeros(size(d2,1),3,3);
                val(:,1,1) = d2(:,1);
                val(:,2,2) = d2(:,2);
                val(:,3,3) = d2(:,3);
                val(:,1,2) = d2(:,4);
                val(:,2,1) = d2(:,4);
                val(:,1,3) = d2(:,5);
                val(:,3,1) = d2(:,5);
                val(:,2,3) = d2(:,6);
                val(:,3,2) = d2(:,6);
            end
        end

    end

end