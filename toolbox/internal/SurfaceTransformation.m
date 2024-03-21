classdef SurfaceTransformation
    %SurfaceTransformation - Lifted surface data and geometric error coefficients
    %   OBJ = SurfaceTransformation(MESH,P) creates the coefficient
    %   functions A, B and f of the surface problem P on the surface
    %   mesh MESH (a MeshSurface), defined through the normal
    %   projection: f is the lifted right-hand side and B the
    %   coefficient of the geometric error estimator
    %   '(B(x)*grad(u))^2dx'. P is a ProblemSurfaceTemplate subclass
    %   with a right-hand side f.
    %
    %   SurfaceTransformation functions:
    %       A - Ratio of curvature-corrected distance factors
    %       B - Geometric error coefficient of the surface mesh
    %       f - Right-hand side lifted onto the surface mesh
    %
    %   SurfaceTransformation properties:
    %       mesh - Surface mesh
    %       p    - Surface problem with the right-hand side f
    %
    %   See also projectVertices2Surface, MeshSurface,
    %   ProblemSurfaceTemplate, EstimatorForm, LinearForm

    % Copyright (c) 2026 by the fem4m authors

    properties
        mesh  % Surface mesh
        p     % Surface problem with the right-hand side f
    end

    methods
        function obj = SurfaceTransformation(mesh,p)
            %SurfaceTransformation - Create the lifted surface coefficients
            obj.mesh = mesh;
            obj.p = p;
        end

        function val = A(obj,x,xref,~)
            %A - Ratio of curvature-corrected distance factors
            %   VAL = A(OBJ,X,XREF,PARTS) returns
            %   max(1-D*K1,1-D*K2)/(|dot(N,N_H)|*min(1-D*K1,1-D*K2))
            %   at the points X (nrElems-by-3, one point per element of MESH
            %   in element order) with the common reference point XREF
            %   (1-by-2), as the forms pass them to a coefficient function;
            %   VAL is nrElems-by-1. D is the signed distance of X to the
            %   surface, K1 and K2 the principal curvatures of the parallel
            %   surface through X, N the surface normal at the projection of
            %   X and N_H the element normal of MESH at XREF. PARTS is
            %   ignored.

            ax = obj.a(x);
            k = obj.eigenvalues(x,ax);
            d = obj.distance(x,ax);
            numerator = max(1-d.*k(:,1), 1-d.*k(:,2));
            v = obj.p.normal(ax);
            v_h = obj.mesh.surfaceNormal4e(xref);
            denominator = abs( dot(v,v_h,2) ) .* min(1-d.*k(:,1), 1-d.*k(:,2));
            val = numerator ./ denominator;
        end

        function val = B(obj,x,xref,~)
            %B - Geometric error coefficient of the surface mesh
            %   VAL = B(OBJ,X,XREF,PARTS) returns
            %   (|D*(K1-K2)| + |1-dot(N,N_H)|*(1+4*M))/MU at the points X
            %   (one row per point, see A), with M = max(1-D*K1,1-D*K2) and
            %   the area factor MU = (1-D*K1)*(1-D*K2)*dot(N,N_H). PARTS is
            %   ignored.

            ax = obj.a(x);
            k = obj.eigenvalues(x,ax);
            d = distance(obj,x,ax);
            v = obj.p.normal(ax);
            v_h = obj.mesh.surfaceNormal4e(xref);
            val = (1 ./ obj.mu(x,ax,xref)) .* ( vecnorm(d.*(k(:,1) - k(:,2)),2,2) + abs(1 - dot(v,v_h,2)).*(1 + 4*max(1-d.*k(:,1), 1-d.*k(:,2))) );
        end

        function val = f(obj,x,xref,~)
            %F - Right-hand side lifted onto the surface mesh
            %   VAL = F(OBJ,X,XREF,PARTS) returns MU times the right-hand
            %   side P.f at the projections AX of the points X onto the
            %   surface (one row per point, see A), with the area factor
            %   MU = (1-D*K1)*(1-D*K2)*dot(N,N_H). PARTS is ignored.

            ax = obj.a(x);
            val = obj.mu(x,ax,xref).*obj.p.f(ax);
        end
    end

    % Auxiliary functions
    methods (Access=private)

        function val = distance(obj,x,ax)
            %DISTANCE - Signed distance of the points to the surface
            val = sign(obj.p.surfaceFunction(x)) .* vecnorm(ax - x,2,2);
        end

        function k = eigenvalues(obj,x,ax)
            %EIGENVALUES - Principal curvatures of the parallel surface
            k = obj.p.curvature(ax);
            d = obj.distance(x,ax);
            k = k./(1 + d.*k);
        end

        function val = mu(obj,x,ax,xref)
            %MU - Area factor (1-d*k1)*(1-d*k2)*dot(n,n_h) of the lift
            k = obj.eigenvalues(x,ax);
            d = obj.distance(x,ax);
            val = (1 - d.*k(:,1)).*(1 - d.*k(:,2)).*dot(obj.p.normal(ax),obj.mesh.surfaceNormal4e(xref),2);
        end

        function val = a(obj,x)
            %A - Normal projection of the points onto the surface
            val = projectVertices2Surface(x, @(z)(obj.p.surfaceFunction(z)), @(z)(obj.p.gradSurfaceFunction(z)));
        end

    end
end