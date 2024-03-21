classdef MeshSurface < Mesh
    %MeshSurface - Triangular surface mesh embedded in 3d
    %   OBJ = MeshSurface(GEOM,SURFACEF,GRADSURFACEF) creates a surface
    %   mesh from the initial triangular mesh GEOM in 3d (a Geometry
    %   object, struct, or file name prefix, see Mesh) and projects its
    %   vertices onto the zero level set of the function handle
    %   SURFACEF with the gradient GRADSURFACEF (see
    %   projectVertices2Surface). SURFACEF(X) returns an N-by-1 vector
    %   and GRADSURFACEF(X) an N-by-3 array for the points X (N-by-3).
    %
    %   OBJ = MeshSurface(c4n,n4e,n4Db,n4Nb,n4Sb,SURFACEF,GRADSURFACEF)
    %   projects the vertices of the explicit initial mesh with the node
    %   coordinates c4n (N-by-3), the node ids n4e of the elements and
    %   the boundary node pairs n4Db, n4Nb and n4Sb (see Mesh) onto the
    %   surface.
    %
    %   OBJ = MeshSurface(...,DEGREE) also represents the surface by
    %   curved isoparametric elements: Lagrange elements (Pk) of degree
    %   DEGREE whose nodes are projected onto the surface. DEGREE = 1
    %   (the default, also for an empty DEGREE) uses the flat
    %   triangles.
    %
    %   Every refinement and coarsening projects the vertices onto the
    %   surface and, for DEGREE > 1, rebuilds the isoparametric map.
    %
    %   MeshSurface functions:
    %       point4e                 - Compute surface points of elements
    %       point4ed                - Compute surface points of edges
    %       area4e                  - Compute the surface areas of elements
    %       length4ed               - Compute the lengths of surface edges
    %       surfaceNormal4e         - Compute the unit surface normals
    %       surfaceGradNormal4e     - Compute gradient of the surface normal
    %       surfaceCurvature4e      - Compute the principal curvatures
    %       surfaceTotalCurvature4e - Compute sum of principal curvatures
    %       normal4ed               - Compute the outward co-normals of edges
    %       tangent4ed              - Compute the unit tangents of edges
    %       normalProjector         - Compute the tangent space projector
    %       transformation          - Compute the transposed pseudo-inverse
    %       gradTransformation      - Compute gradient of the transformation
    %       refine                  - Refine the mesh and reproject vertices
    %       plotSurfaceNormal       - Plot the surface normals as arrows
    %       plotCoNormal            - Plot the edge co-normals as arrows
    %       plotEnumeration         - Label nodes, elements, and edges with ids
    %       jacobian                - Compute the Jacobian of the embedding
    %       gradJacobian            - Compute second derivatives of embedding
    %       DJacobian               - Compute surface gradients of Jacobian
    %       jacobian4ed             - Compute the Jacobian along edges
    %
    %   MeshSurface properties:
    %       surfaceFunction     - Level-set function of the surface
    %       gradSurfaceFunction - Gradient of the level-set function
    %       pk                  - Pk element of the isoparametric map
    %       degree              - Degree of the surface interpolation
    %
    %   See also Mesh, Pk, projectVertices2Surface, ProblemSurfaceTemplate

    % Copyright (c) 2026 by the fem4m authors

    properties
        surfaceFunction     % level-set function of the surface
        gradSurfaceFunction % gradient of the level-set function
        pk                  % Pk element of the isoparametric map
        degree              % degree of the surface interpolation
    end

    methods

        function obj = MeshSurface(c4n,n4e,n4Db,n4Nb,n4Sb,surfaceF,gradSurfaceF,degree,varargin)
            %MeshSurface - Create a surface mesh object
            if nargin > 8 || (nargin > 5 && ~isa(surfaceF,'function_handle'))
                error('fem4m:Mesh:pkInterpolationRemoved', ...
                    ['MeshSurface takes no prolongation object anymore: use ', ...
                    'MeshSurface(c4n,n4e,n4Db,n4Nb,n4Sb,SURFACEF,GRADSURFACEF[,DEGREE]).']);
            end
            if nargin <= 4
                surfaceF = n4e;
                gradSurfaceF = n4Db;
                args{1} = c4n;
                if nargin == 4 && ~isempty(n4Nb)
                    degree = n4Nb;
                else
                    degree = 1;
                end
            else
                args = {c4n,n4e,n4Db,n4Nb,n4Sb};
                if nargin < 8 || isempty(degree)
                    degree = 1;
                end
            end
            obj@Mesh(args{:});
            obj.surfaceFunction = surfaceF;
            obj.gradSurfaceFunction = gradSurfaceF;
            obj.c4n = projectVertices2Surface(obj.c4n,obj.surfaceFunction,obj.gradSurfaceFunction);
            obj.degree = degree;
            if degree > 1
                obj.pk = Pk(degree,obj);
                obj.pk.c4dof = projectVertices2Surface(obj.pk.c4dof,obj.surfaceFunction,obj.gradSurfaceFunction);
            end
        end

        function val = point4e(obj,refp)
            %point4e - Compute surface points of elements
            %   VAL = point4e(OBJ,REFP) computes for every element the
            %   point (nrElems-by-3) with the reference coordinates REFP
            %   (1-by-2): on the flat triangle for DEGREE 1, else on the
            %   curved isoparametric element.
            if obj.degree == 1
                val = point4e@Mesh(obj,refp);
            else
                basis = obj.pk.monomBasis(refp(1),refp(2))';
                basis = obj.pk.refBasisCoeff*basis;
                val(:,1) = reshape(obj.pk.c4dof(obj.pk.dof4e(:),1),obj.nrElems,[])*basis;
                val(:,2) = reshape(obj.pk.c4dof(obj.pk.dof4e(:),2),obj.nrElems,[])*basis;
                val(:,3) = reshape(obj.pk.c4dof(obj.pk.dof4e(:),3),obj.nrElems,[])*basis;
            end
        end

        function val = point4ed(obj,refp)
            %point4ed - Compute surface points of edges
            %   VAL = point4ed(OBJ,REFP) computes for every edge the point
            %   (nrEdges-by-3) with the reference coordinate REFP along
            %   the edge, 0 at its first and 1 at its second node: on the
            %   straight edge for DEGREE 1, else on the curved edge of
            %   the element e4ed(:,1).
            if obj.degree == 1
                val = point4ed@Mesh(obj,refp);
            else
                basis = getEdgeValues4RefElemFunktion(refp(1),@(x)obj.pk.getBasis([],x,[]),obj);
                basis = basis(:,:,1);
                val(:,1) = sum(reshape(obj.pk.c4dof(obj.pk.dof4e(obj.e4ed(:,1),:),1),obj.nrEdges,[]).*basis,2);
                val(:,2) = sum(reshape(obj.pk.c4dof(obj.pk.dof4e(obj.e4ed(:,1),:),2),obj.nrEdges,[]).*basis,2);
                val(:,3) = sum(reshape(obj.pk.c4dof(obj.pk.dof4e(obj.e4ed(:,1),:),3),obj.nrEdges,[]).*basis,2);
            end
        end

        function area4e = area4e(obj,refp)
            %area4e - Compute the surface areas of elements
            %   AREA4E = area4e(OBJ) computes the areas (nrElems-by-1) of
            %   the flat triangles spanned by the vertices in 3d.
            %
            %   AREA4E = area4e(OBJ,REFP) computes for DEGREE > 1 half the
            %   area element |J1 x J2| of the curved elements at the
            %   reference point REFP (1-by-2), with the rows J1 and J2 of
            %   jacobian(OBJ,REFP); DEGREE 1 ignores REFP.
            % surface area in 3d
            if obj.degree == 1 || nargin < 2
                area4e = 1/2 * sqrt(sum((cross(obj.c4n(obj.n4e(:,3),:)-obj.c4n(obj.n4e(:,1),:),...
                    obj.c4n(obj.n4e(:,2),:)-obj.c4n(obj.n4e(:,1),:))).^2,2));
            else
                grad = permute(obj.jacobian(refp),[1 3 2]);
                area4e = 1/2 * sqrt(sum(cross(grad(:,:,1),grad(:,:,2)).^2,2));
            end
        end

        function val = length4ed(obj,refp)
            %length4ed - Compute the lengths of surface edges
            %   VAL = length4ed(OBJ) computes the lengths (nrEdges-by-1)
            %   of the straight edges in 3d.
            %
            %   VAL = length4ed(OBJ,REFP) computes for DEGREE > 1 the
            %   length element |jacobian4ed(OBJ,REFP)| of the curved
            %   edges at the reference coordinate REFP along the edge;
            %   DEGREE 1 ignores REFP.
            if obj.degree == 1 || nargin < 2
                val = length4ed@Mesh(obj);
            else
                grad = obj.jacobian4ed(refp);
                val = sqrt(sum(grad.^2,2));
            end
        end

        function normal4e = surfaceNormal4e(obj,refp)
            %surfaceNormal4e - Compute the unit surface normals
            %   NORMAL4E = surfaceNormal4e(OBJ) computes the unit normal
            %   (nrElems-by-3) of every flat triangle, (P2-P1)x(P3-P1)
            %   normalized for its nodes P1, P2 and P3; it points outwards
            %   when the nodes are ordered counterclockwise seen from
            %   outside.
            %
            %   NORMAL4E = surfaceNormal4e(OBJ,REFP) computes for
            %   DEGREE > 1 the unit normal of the curved elements at the
            %   reference point REFP (1-by-2), with the same orientation;
            %   DEGREE 1 ignores REFP.
            % unit normal (P2-P1)x(P3-P1), outward for a counterclockwise
            % node order seen from outside
            if obj.degree == 1 || nargin < 2
                normal4e = cross(obj.c4n(obj.n4e(:,3),:)-obj.c4n(obj.n4e(:,1),:),...
                    obj.c4n(obj.n4e(:,2),:)-obj.c4n(obj.n4e(:,1),:));
                normal4e = -normal4e./((sqrt(sum(normal4e.^2,2)))*ones(1,3));
            else
                grad = permute(obj.jacobian(refp),[1 3 2]);
                normal4e = cross(grad(:,:,1),grad(:,:,2));
                normal4e = normal4e./((sqrt(sum(normal4e.^2,2)))*ones(1,3));
            end
        end

        function gradNormal4e = surfaceGradNormal4e(obj,refp)
            %surfaceGradNormal4e - Compute gradient of the surface normal
            %   GRADNORMAL4E = surfaceGradNormal4e(OBJ,REFP) computes for
            %   every element the surface gradient (nrElems-by-3-by-3,
            %   entry (i,j) = d_j n_i) of the unit normal n of the curved
            %   element at the reference point REFP (1-by-2); it is zero
            %   for DEGREE 1.
            if obj.degree == 1
                gradNormal4e = zeros(obj.nrElems,3,3);
            else
                grad = permute(obj.jacobian(refp),[1 3 2]);
                normal4e = cross(grad(:,:,1),grad(:,:,2));
                normgrad = sqrt(sum(normal4e.^2,2));
                d2 = permute(obj.gradJacobian(refp),[1 3 2]);
                ct = obj.transformation(refp);
                tmp = (cross(d2(:,:,1),grad(:,:,2)) + cross(grad(:,:,1),d2(:,:,3)));
                gradNormal4e(:,:,1) = tmp./normgrad(:,[1 1 1]) - normal4e.*((sum(normal4e.*tmp,2)./(normgrad.^3))*ones(1,3));
                tmp = cross(grad(:,:,1),d2(:,:,2)) + cross(d2(:,:,3),grad(:,:,2));
                gradNormal4e(:,:,2) = tmp./normgrad(:,[1 1 1]) - normal4e.*((sum(normal4e.*tmp,2)./(normgrad.^3))*ones(1,3));
                gradNormal4e = matMul(gradNormal4e,permute(ct,[1 3 2]));
            end
        end

        function val = surfaceCurvature4e(obj,refp)
            %surfaceCurvature4e - Compute the principal curvatures
            %   VAL = surfaceCurvature4e(OBJ,REFP) computes the principal
            %   curvatures [k1,k2] (nrElems-by-2) of every element at the
            %   reference point REFP (1-by-2): the two eigenvalues of
            %   largest modulus of the shape operator, the surface
            %   gradient of the normal (surfaceGradNormal4e) times the
            %   tangent projector, sorted by decreasing modulus. They are
            %   positive on a convex surface with outward normals and
            %   vanish for DEGREE 1.
            if obj.degree == 1
                val = zeros(obj.nrElems,2);
            else
                normals = obj.surfaceNormal4e(refp);
                P = -matMul(normals,reshape(normals,size(normals,1),1,3));
                P(:,1,1) = P(:,1,1) + 1;
                P(:,2,2) = P(:,2,2) + 1;
                P(:,3,3) = P(:,3,3) + 1;
                gradNormal4e = obj.surfaceGradNormal4e(refp);
                A = permute(matMul(P,permute(gradNormal4e,[1 3 2])),[3 2 1]);
                val = zeros(size(A,3),size(A,1));
                for curElem = 1 : size(A,3)
                    val(curElem,:) = eig(A(:,:,curElem));
                end
                val= sort(val,2,'descend','ComparisonMethod','abs');
                val = val(:,[1 2]);
            end
        end

        function val = surfaceTotalCurvature4e(obj,refp)
            %surfaceTotalCurvature4e - Compute sum of principal curvatures
            %   VAL = surfaceTotalCurvature4e(OBJ,REFP) computes the sum
            %   k1+k2 (nrElems-by-1) of the principal curvatures of every
            %   element at the reference point REFP (1-by-2), the trace of
            %   surfaceGradNormal4e minus its normal-normal component; it
            %   vanishes for DEGREE 1.
            if obj.degree == 1
                val = zeros(obj.nrElems,1);
            else
                normals = obj.surfaceNormal4e(refp);
                gradNormal4e = obj.surfaceGradNormal4e(refp);
                val = gradNormal4e(:,1,1)+gradNormal4e(:,2,2)+gradNormal4e(:,3,3)...
                    - matMul(permute(normals,[1 3 2]),matMul(gradNormal4e,normals));
            end
        end

        function normal4ed = normal4ed(obj,refp)
            %normal4ed - Compute the outward co-normals of edges
            %   NORMAL4ED = normal4ed(OBJ) computes for every edge the
            %   unit co-normals (nrEdges-by-2-by-3) of its two triangles:
            %   NORMAL4ED(ED,K,:) lies in the plane of the triangle
            %   e4ed(ED,K), is orthogonal to the edge and points out of
            %   that triangle; it is zero for K = 2 on boundary edges.
            %
            %   NORMAL4ED = normal4ed(OBJ,REFP) computes for DEGREE > 1
            %   the co-normals of the curved elements at the reference
            %   coordinate REFP along the edge (see point4ed); DEGREE 1
            %   ignores REFP.
            % compute conormal for surface edges for T_1 and T_2
            if obj.degree == 1 || nargin < 2
                normal = obj.surfaceNormal4e;
                tangent = obj.tangent4ed;
                normal4ed(:,1,:) = cross(normal(obj.e4ed(:,1),:),-tangent);
                normal4ed(obj.e4ed(:,2)~=0,2,:) = cross(normal(obj.e4ed(obj.e4ed(:,2)~=0,2),:),tangent(obj.e4ed(:,2)~=0,:));
            else
                normal = getEdgeValues4RefElemFunktion(refp,@(pts)obj.surfaceNormal4e(pts),obj);
                tangent = obj.tangent4ed(refp);
                normal4ed(:,1,:) = cross(normal(:,:,1),-tangent);
                normal4ed(:,2,:) = cross(normal(:,:,2),tangent);
            end
        end

        function val = tangent4ed(obj,refp)
            %tangent4ed - Compute the unit tangents of edges
            %   VAL = tangent4ed(OBJ) computes the unit tangents
            %   (nrEdges-by-3) of the straight edges, pointing from the
            %   first to the second node in n4ed.
            %
            %   VAL = tangent4ed(OBJ,REFP) computes for DEGREE > 1 the
            %   normalized jacobian4ed(OBJ,REFP) of the curved edges at
            %   the reference coordinate REFP; DEGREE 1 ignores REFP.
            if obj.degree == 1 || nargin < 2
                val = tangent4ed@Mesh(obj);
            else
                grad = obj.jacobian4ed(refp);
                length4ed = sqrt(sum(grad.^2,2));
                val = grad./length4ed;
            end
        end

        function P = normalProjector(obj,refp)
            %normalProjector - Compute the tangent space projector
            %   P = normalProjector(OBJ,REFP) computes for every element
            %   the orthogonal projection I - n*n' (nrElems-by-3-by-3)
            %   onto the tangent plane, with the unit normal n of
            %   surfaceNormal4e(OBJ,REFP) at the reference point REFP
            %   (1-by-2).
            normals = obj.surfaceNormal4e(refp);
            P = -matMul(normals,permute(normals,[1 3 2]));
            P(:,1,1) = P(:,1,1) + 1;
            P(:,2,2) = P(:,2,2) + 1;
            P(:,3,3) = P(:,3,3) + 1;
        end

        function val = transformation(obj,refp)
            %TRANSFORMATION - Compute the transposed pseudo-inverse
            %   VAL = TRANSFORMATION(OBJ,REFP) computes for every element
            %   the transpose (nrElems-by-3-by-2) of the pseudo-inverse
            %   (A'*A)^(-1)*A' of the Jacobian A (3-by-2) of the map from
            %   the reference triangle to the element at the reference
            %   point REFP (1-by-2); VAL(e,:,:) times a reference
            %   gradient is the surface gradient. DEGREE 1 ignores REFP.
            if obj.degree == 1
                P1 = obj.c4n(obj.n4e(:,1),:);
                P2 = obj.c4n(obj.n4e(:,2),:);
                P3 = obj.c4n(obj.n4e(:,3),:);
                x2 = P2(:,1)-P1(:,1); y2 = P2(:,2)-P1(:,2); z2 = P2(:,3)-P1(:,3);
                x3 = P3(:,1)-P1(:,1); y3 = P3(:,2)-P1(:,2); z3 = P3(:,3)-P1(:,3);
            else
                %Jacobian matrix
                grad = obj.jacobian(refp);
                x2 = grad(:,1,1); y2 = grad(:,1,2); z2 = grad(:,1,3);
                x3 = grad(:,2,1); y3 = grad(:,2,2); z3 = grad(:,2,3);
            end
            %compute the pseudo-inverse via
            % (A^T*A)^(-1)*A^T
            d = (x2.^2.*y3.^2 + x2.^2.*z3.^2 - 2.*x2.*x3.*y2.*y3 - 2.*x2.*x3.*z2.*z3 + x3.^2.*y2.^2 + x3.^2.*z2.^2 + y2.^2.*z3.^2 - 2.*y2.*y3.*z2.*z3 + y3.^2.*z2.^2);
            val = [x2.*(x3.^2 + y3.^2 + z3.^2) - x3.*(x2.*x3 + y2.*y3 + z2.*z3), y2.*(x3.^2 + y3.^2 + z3.^2) - y3.*(x2.*x3 + y2.*y3 + z2.*z3), z2.*(x3.^2 + y3.^2 + z3.^2) - z3.*(x2.*x3 + y2.*y3 + z2.*z3),...
                x3.*(x2.^2 + y2.^2 + z2.^2) - x2.*(x2.*x3 + y2.*y3 + z2.*z3), y3.*(x2.^2 + y2.^2 + z2.^2) - y2.*(x2.*x3 + y2.*y3 + z2.*z3), z3.*(x2.^2 + y2.^2 + z2.^2) - z2.*(x2.*x3 + y2.*y3 + z2.*z3)]...
                ./d;
            val = reshape(val,obj.nrElems,3,2);
        end

        function val = gradTransformation(obj,refp)
            %gradTransformation - Compute gradient of the transformation
            %   VAL = gradTransformation(OBJ,REFP) computes for every
            %   element the derivatives (nrElems-by-2-by-3-by-2) of
            %   transformation(OBJ,REFP) with respect to the reference
            %   coordinates at the reference point REFP (1-by-2):
            %   VAL(e,k,:,:) is the derivative in the k-th reference
            %   direction. It is zero for DEGREE 1.
            if obj.degree == 1
                val = zeros(obj.nrElems,2,3,3);
            else
                % Jacobian matrix
                grad = obj.jacobian(refp);
                x2 = grad(:,[1 1],1); y2 = grad(:,[1 1],2); z2 = grad(:,[1 1],3);
                x3 = grad(:,[2 2],1); y3 = grad(:,[2 2],2); z3 = grad(:,[2 2],3);
                det = x2.^2.*y3.^2 + x2.^2.*z3.^2 - 2.*x2.*x3.*y2.*y3 - 2.*x2.*x3.*z2.*z3 + x3.^2.*y2.^2 + x3.^2.*z2.^2 + y2.^2.*z3.^2 - 2.*y2.*y3.*z2.*z3 + y3.^2.*z2.^2;
                %derivative of the entries of the Jacobian
                d2 = obj.gradJacobian(refp);
                dx2 = d2(:,[1 3],1); dy2 = d2(:,[1 3],2); dz2 = d2(:,[1 3],3);
                dx3 = d2(:,[3 2],1); dy3 = d2(:,[3 2],2); dz3 = d2(:,[3 2],3);
                %derivative of the pseudo-inverse
                %Product rule for the inverse of det
                ddet = 2*( dx2.*x2.*y3.^2 + dx2.*x2.*z3.^2 + dx3.*x3.*y2.^2 + dx3.*x3.*z2.^2 + dy2.*y2.*z3.^2 + dy3.*y3.*z2.^2 ...
                    + x2.^2.*dy3.*y3 + x2.^2.*z3.*dz3 + x3.^2.*y2.*dy2 + x3.^2.*z2.*dz2 + y2.^2.*z3.*dz3 + y3.^2.*z2.*dz2)...
                    -2*(dx2.*x3.*y2.*y3 + dx2.*x3.*z2.*z3 + dy2.*y3.*z2.*z3 ...
                    +x2.*dx3.*y2.*y3 + x2.*dx3.*z2.*z3 + y2.*dy3.*z2.*z3 ...
                    +x2.*x3.*dy2.*y3 + x2.*x3.*dz2.*z3 + y2.*y3.*dz2.*z3 ...
                    +x2.*x3.*y2.*dy3 + x2.*x3.*z2.*dz3 + y2.*y3.*z2.*dz3);
                ddet = -ddet./(det.^2);
                %Product rule
                r3 = (x3.^2 + y3.^2 + z3.^2)./det;
                r2 = (x2.^2 + y2.^2 + z2.^2)./det;
                dr3 = (dx3.*x3 + dy3.*y3 + dz3.*z3)./det;
                dr2 = (dx2.*x2 + dy2.*y2 + dz2.*z2)./det;
                m = (x2.*x3 + y2.*y3 + z2.*z3)./det;
                dm = (dx2.*x3 + dy2.*y3 + dz2.*z3)./det;
                md = (x2.*dx3 + y2.*dy3 + z2.*dz3)./det;
                ddr3 = (x3.^2 + y3.^2 + z3.^2).*ddet;
                ddr2 = (x2.^2 + y2.^2 + z2.^2).*ddet;
                ddm = (x2.*x3 + y2.*y3 + z2.*z3).*ddet;
                val = [dx2.*r3 - dx3.*m, dy2.*r3 - dy3.*m, dz2.*r3 - dz3.*m,...
                    dx3.*r2 - dx2.*m, dy3.*r2 - dy2.*m, dz3.*r2 - dz2.*m]...
                    +[x2.*dr3 - x3.*dm, y2.*dr3 - y3.*dm, z2.*dr3 - z3.*dm,...
                    x3.*dr2 - x2.*dm, y3.*dr2 - y2.*dm, z3.*dr2 - z2.*dm]...
                    +[x2.*dr3 - x3.*md, y2.*dr3 - y3.*md, z2.*dr3 - z3.*md,...
                    x3.*dr2 - x2.*md, y3.*dr2 - y2.*md, z3.*dr2 - z2.*md]...
                    +[x2.*ddr3 - x3.*ddm, y2.*ddr3 - y3.*ddm, z2.*ddr3 - z3.*ddm,...
                    x3.*ddr2 - x2.*ddm, y3.*ddr2 - y2.*ddm, z3.*ddr2 - z2.*ddm];
                val = reshape(val,obj.nrElems,2,3,2);
            end
        end

        function obj = refine(obj,varargin)
            %REFINE - Refine the mesh and reproject vertices
            %   OBJ = REFINE(OBJ) refines the mesh along the marked edges
            %   markedEdges by red-green-blue refinement, or uniformly
            %   with a warning when no edge is marked.
            %
            %   OBJ = REFINE(OBJ,MARKED4ED) marks the edges given by the
            %   0/1 or logical vector MARKED4ED (it replaces markedEdges)
            %   before refining; a vector of length nrElems marks the
            %   three edges of every selected element instead.
            %
            %   OBJ = REFINE(OBJ,'uniform',N), 'RGB' and 'NVB' work as in
            %   Mesh.refine. 'graded' and 'gradedred' measure the distances to
            %   the grading points and find them in the x-y coordinates only;
            %   'gradedbeta' raises an error on a surface mesh as soon as a node
            %   lies closer than R to a grading point (its shift has only two
            %   coordinates).
            %
            %   After every refinement round the vertices are projected
            %   onto the surface and, for DEGREE > 1, the isoparametric
            %   map is rebuilt with its nodes on the surface.
            obj = refine@Mesh(obj,varargin{:});
        end

        function obj = plotSurfaceNormal(obj)
            %plotSurfaceNormal - Plot the surface normals as arrows
            %   plotSurfaceNormal(OBJ) plots the mesh into the cleared
            %   figure together with an arrow for the unit normal at the
            %   image of the reference centroid of every element.
            normals = obj.surfaceNormal4e([1/3,1/3]);
            points = obj.point4e([1/3,1/3]);
            obj.plot;
            hold on;
            quiver3(points(:,1),points(:,2),points(:,3),normals(:,1),normals(:,2),normals(:,3));
        end

        function obj = plotCoNormal(obj)
            %plotCoNormal - Plot the edge co-normals as arrows
            %   plotCoNormal(OBJ) plots the mesh into the cleared figure
            %   together with arrows for both co-normals (normal4ed) at
            %   the midpoint of every edge.
            normals = permute(obj.normal4ed(1/2),[1 3 2]);
            points = obj.point4ed(1/2);
            obj.plot;
            hold on;
            quiver3(points(:,1),points(:,2),points(:,3),normals(:,1,1),normals(:,2,1),normals(:,3,1));
            quiver3(points(:,1),points(:,2),points(:,3),normals(:,1,2),normals(:,2,2),normals(:,3,2));
        end

        function obj = plotEnumeration(obj)
            %plotEnumeration - Label nodes, elements, and edges with ids
            %   plotEnumeration(OBJ) writes the ids of the nodes, of the
            %   elements (at the images of the reference centroids) and
            %   of the edges (at their midpoints) in different colors
            %   into the current 3d axes; it does not plot the mesh.
            mid4e = obj.point4e([1/3,1/3]);
            mid4ed = obj.point4ed(1/2);
            text(obj.c4n(:,1),obj.c4n(:,2),obj.c4n(:,3),num2str((1:size(obj.c4n,1))'),...
                'color',[0.26 0.49 0.19],'FontWeight','bold','FontSize',24,'HorizontalAlignment','center');
            text(mid4e(:,1),mid4e(:,2),mid4e(:,3),num2str((1:size(obj.n4e,1))'),...
                'color',[0.95 0.91 0],'FontWeight','bold','FontSize',24,'HorizontalAlignment','center');
            text(mid4ed(:,1),mid4ed(:,2),mid4ed(:,3),num2str((1:obj.nrEdges)'),...
                'color',[0.38 0.02 0.68],'FontWeight','bold','FontSize',24,'HorizontalAlignment','center');
        end

    end

    methods (Access = protected)
        function obj = finishRound(obj)
            %finishRound - Reproject the vertices after a change of the mesh
            %   OBJ = finishRound(OBJ) projects the vertices onto the
            %   surface and, for DEGREE > 1, rebuilds the isoparametric
            %   map of the current elements with its nodes on the surface.
            obj = finishRound@Mesh(obj);
            obj.c4n = projectVertices2Surface(obj.c4n,obj.surfaceFunction,obj.gradSurfaceFunction);
            if obj.degree > 1
                obj.pk = []; % no chain of meshes inside pk
                obj.pk = Pk(obj.degree,obj);
                obj.pk.c4dof = projectVertices2Surface(obj.pk.c4dof,obj.surfaceFunction,obj.gradSurfaceFunction);
            end
        end
    end

    methods(Access=public)

        function val = jacobian(obj,refp)
            %JACOBIAN - Compute the Jacobian of the embedding
            %   VAL = JACOBIAN(OBJ,REFP) computes for every element the
            %   transposed Jacobian (nrElems-by-2-by-3) of the map from
            %   the reference triangle to the element at the reference
            %   point REFP (1-by-2): VAL(e,k,:) is the derivative with
            %   respect to the k-th reference coordinate, P2-P1 and P3-P1
            %   for the nodes P1, P2, P3 and DEGREE 1, which ignores REFP.
            if obj.degree == 1
                P1 = obj.c4n(obj.n4e(:,1),:);
                P2 = obj.c4n(obj.n4e(:,2),:);
                P3 = obj.c4n(obj.n4e(:,3),:);
                x2 = P2(:,1)-P1(:,1); y2 = P2(:,2)-P1(:,2); z2 = P2(:,3)-P1(:,3);
                x3 = P3(:,1)-P1(:,1); y3 = P3(:,2)-P1(:,2); z3 = P3(:,3)-P1(:,3);
                val = zeros(obj.nrElems,2,3);
                val(:,1,1) = x2;
                val(:,1,2) = y2;
                val(:,1,3) = z2;
                val(:,2,1) = x3;
                val(:,2,2) = y3;
                val(:,2,3) = z3;
            else
                m = obj.pk.monomBasis(refp(1),refp(2)).';
                val = reshape(obj.pk.gradBasisCoefficients*m,obj.pk.nrBasisFunc,2);
                val = matMul(permute(val,[3 2 1]),reshape(obj.pk.c4dof(obj.pk.dof4e,:),obj.nrElems,obj.pk.nrBasisFunc,3));
            end
        end

        function val = gradJacobian(obj,refp)
            %gradJacobian - Compute second derivatives of embedding
            %   VAL = gradJacobian(OBJ,REFP) computes for every element
            %   the second derivatives (nrElems-by-3-by-3) of the map
            %   from the reference triangle at the reference point REFP
            %   (1-by-2): VAL(e,:,c) are the derivatives [xx yy xy] of
            %   the coordinate c with respect to the reference
            %   coordinates. They vanish for DEGREE 1.
            % 2d entries [dxx dyy dxy]*[x y z]
            if obj.degree == 1
                val = zeros(obj.nrElems,3,3);
            else
                m = obj.pk.monomBasis(refp(1),refp(2)).';
                val = reshape(obj.pk.d2BasisCoefficients*m,obj.pk.nrBasisFunc,3);
                val = matMul(permute(val,[3 2 1]),reshape(obj.pk.c4dof(obj.pk.dof4e,:),obj.nrElems,obj.pk.nrBasisFunc,3));
            end
        end

        function val = DJacobian(obj,refp)
            %DJacobian - Compute surface gradients of Jacobian
            %   VAL = DJacobian(OBJ,REFP) computes for every element the
            %   surface gradients (nrElems-by-2-by-3-by-3) of the entries
            %   of jacobian(OBJ,REFP) at the reference point REFP
            %   (1-by-2): VAL(e,k,c,:) is the gradient of the entry
            %   (k,c). They vanish for DEGREE 1.
            % surface gradients in 3d (last dimension)
            dJ = obj.gradJacobian(refp);
            dx2 = dJ(:,[1 3],1); dy2 = dJ(:,[1 3],2); dz2 = dJ(:,[1 3],3);
            dx3 = dJ(:,[3 2],1); dy3 = dJ(:,[3 2],2); dz3 = dJ(:,[3 2],3);
            T = obj.transformation(refp);
            val = zeros(obj.nrElems,2,3,3);
            val(:,1,1,:) = matMul(T,dx2);
            val(:,1,2,:) = matMul(T,dy2);
            val(:,1,3,:) = matMul(T,dz2);
            val(:,2,1,:) = matMul(T,dx3);
            val(:,2,2,:) = matMul(T,dy3);
            val(:,2,3,:) = matMul(T,dz3);
        end

        function val = jacobian4ed(obj,refp)
            %jacobian4ed - Compute the Jacobian along edges
            %   VAL = jacobian4ed(OBJ,REFP) computes for every edge the
            %   derivative (nrEdges-by-3) of the map from the reference
            %   edge [0,1] to the edge, from its first to its second node,
            %   at the reference coordinate REFP; it is P2-P1 for the
            %   nodes P1, P2 in n4ed and DEGREE 1, which ignores REFP.
            if obj.degree == 1
                P1 = obj.c4n(obj.n4ed(:,1),:);
                P2 = obj.c4n(obj.n4ed(:,2),:);
                x2 = P2(:,1)-P1(:,1); y2 = P2(:,2)-P1(:,2); z2 = P2(:,3)-P1(:,3);
                val = zeros(obj.nrEdges,3);
                val(:,1) = x2;
                val(:,2) = y2;
                val(:,3) = z2;
            else
                % derivatives along the reference edge (t,0) of the
                % basis functions of the edge dofs
                m = obj.pk.monomBasis(refp(1),0).';
                val = reshape(obj.pk.gradBasisCoefficients*m,obj.pk.nrBasisFunc,2);
                refEdgeIndx = [1,2,4:3:3*(obj.pk.degree)];
                val = val(refEdgeIndx,1);
                val = matMul(permute(val,[3 2 1]),reshape(obj.pk.c4dof(obj.pk.dof4ed,:),obj.nrEdges,2+(obj.pk.degree-1),3));
                val = permute(val,[1 3 2]);
            end
        end

    end

end