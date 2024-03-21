classdef ProblemSurfaceTemplate
    %ProblemSurfaceTemplate - Abstract surface problem
    %   ProblemSurfaceTemplate is an abstract template for problems on a
    %   surface in 3d, the zero level set of surfaceFunction. Subclasses
    %   define the level-set function, its gradient and its Hessian, and a
    %   property geom with the initial mesh; most problems with data also
    %   inherit a PDE template such as ProblemScalarEllipticTemplateExact.
    %   MeshSurface takes @p.surfaceFunction and @p.gradSurfaceFunction to
    %   project the mesh onto the surface.
    %
    %   The functions are evaluated at the points X (N-by-3, one row per
    %   point) and return an N-by-1 vector (surfaceFunction), an N-by-3
    %   array (gradSurfaceFunction, normal) or an N-by-3-by-3 array
    %   (d2SurfaceFunction with entry (i,j) = d_i d_j of the level-set
    %   function, normalProjector). curvature returns an N-by-2 array of
    %   the two principal curvatures, the eigenvalues of largest magnitude
    %   of P*H/|grad| (P the projection, H the Hessian, grad the gradient
    %   of the level-set function) sorted by decreasing magnitude; they are
    %   positive on a convex surface whose level-set gradient points
    %   outwards.
    %
    %   ProblemSurfaceTemplate functions:
    %       surfaceFunction     - Level-set function of the surface (abstract)
    %       gradSurfaceFunction - Gradient of the level-set function (abstract)
    %       d2SurfaceFunction   - Hessian of the level-set function (abstract)
    %       normal              - Unit normal (normalized level-set gradient)
    %       curvature           - Two principal curvatures of the surface
    %       normalProjector     - Projection I - n*n' onto the tangent plane
    %
    %   See also MeshSurface, problemSurfaceSphere, problemSurfaceEllipsoid,
    %   problemSurfaceTorus, ProblemScalarEllipticTemplateExact

    % Copyright (c) 2026 by the fem4m authors

    methods(Abstract)
        val = surfaceFunction(obj,x);
        val = gradSurfaceFunction(obj,x);
        val = d2SurfaceFunction(obj,x);
    end

    methods

        function val = normal(obj,x)
            gradVal = obj.gradSurfaceFunction(x);
            normGrad = vecnorm(gradVal,2,2);
            val = gradVal./normGrad;
        end

        function val = curvature(obj,x)
            gradVal = obj.gradSurfaceFunction(x);
            normGrad = vecnorm(gradVal,2,2);
            d2Val = obj.d2SurfaceFunction(x);
            val(:,1,1) = d2Val(:,1,1)./normGrad - gradVal(:,1).*sum(d2Val(:,1,:),3)./(normGrad.^3);
            val(:,1,2) = d2Val(:,1,2)./normGrad - gradVal(:,2).*sum(d2Val(:,1,:),3)./(normGrad.^3);
            val(:,1,3) = d2Val(:,1,3)./normGrad - gradVal(:,3).*sum(d2Val(:,1,:),3)./(normGrad.^3);
            val(:,2,1) = d2Val(:,2,1)./normGrad - gradVal(:,1).*sum(d2Val(:,2,:),3)./(normGrad.^3);
            val(:,2,2) = d2Val(:,2,2)./normGrad - gradVal(:,2).*sum(d2Val(:,2,:),3)./(normGrad.^3);
            val(:,2,3) = d2Val(:,2,3)./normGrad - gradVal(:,3).*sum(d2Val(:,2,:),3)./(normGrad.^3);
            val(:,3,1) = d2Val(:,3,1)./normGrad - gradVal(:,1).*sum(d2Val(:,3,:),3)./(normGrad.^3);
            val(:,3,2) = d2Val(:,3,2)./normGrad - gradVal(:,2).*sum(d2Val(:,3,:),3)./(normGrad.^3);
            val(:,3,3) = d2Val(:,3,3)./normGrad - gradVal(:,3).*sum(d2Val(:,3,:),3)./(normGrad.^3);
            %projection on tangent plane (P = I - n*n^T)
            normals = gradVal./normGrad;
            P = -matMul(normals,reshape(normals,size(normals,1),1,3));
            P(:,1,1) = P(:,1,1) + 1;
            P(:,2,2) = P(:,2,2) + 1;
            P(:,3,3) = P(:,3,3) + 1;
            A = permute(matMul(P,permute(val,[1 3 2])),[3 2 1]);
            val = zeros(size(A,3),size(A,1));
            for curElem = 1 : size(A,3)
                val(curElem,:) = eig(A(:,:,curElem));
            end
            val = sort(val,2,'descend','ComparisonMethod','abs');
            val = val(:,[1 2]);
        end

        function P = normalProjector(obj,x)
            %normalProjector - Projection I - n*n' onto the tangent plane
            %projection on tangent plane (P = I - n*n^T)
            normals = obj.normal(x);
            P = -matMul(normals,reshape(normals,size(normals,1),1,3));
            P(:,1,1) = P(:,1,1) + 1;
            P(:,2,2) = P(:,2,2) + 1;
            P(:,3,3) = P(:,3,3) + 1;
        end

    end

end