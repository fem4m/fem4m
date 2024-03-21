classdef DGPkFace < Pk3d
    %DGPkFace - Surface Pk space from the first-face basis of Pk3d
    %   OBJ = DGPkFace(DEGREE,MESH) creates the surface finite
    %   element space on the faces of the tetrahedral mesh MESH
    %   whose basis functions are the basis functions of the
    %   first face (z = 0) of the Pk3d element of polynomial degree
    %   DEGREE >= 1, i.e. the 2d Lagrange basis of degree DEGREE on
    %   each face. Each face has its own (DEGREE+1)*(DEGREE+2)/2
    %   degrees of freedom, numbered face by face in dof4f, with no
    %   continuity across the edges of the faces. The face functions
    %   return the values on both sides of each face (the second
    %   side is zero on boundary faces). Only the face functions
    %   listed below apply; the element routines inherited from Pk3d
    %   do not.
    %
    %   DGPkFace functions:
    %       getFaceBasis          - Evaluate the surface basis functions
    %       getFaceValue          - Evaluate the surface FE function
    %       getFaceGradBasis      - Evaluate surface gradients of the basis
    %       getFaceGradValue      - Evaluate surface gradients of the function
    %       getFaceCurlBasis      - Evaluate surface curls of the basis
    %       getFaceCurlValue      - Evaluate surface curls of the function
    %       getFaceTransformation - Build the face transformation
    %
    %   DGPkFace properties:
    %       firstFaceIndices      - Basis functions of the first face
    %       firstEdgeIndices      - Corner and first edge basis functions
    %       normal4f              - Face normals of the mesh
    %
    %   See also Pk3d, FEFace, DGDkFace, DGDkEdge
    %
    %   Note: The surface curl is curl_gamma(PHI) = grad_gamma(PHI) x n,
    %   which is the 2d curl (d_y PHI, -d_x PHI) for n = e_3 and the
    %   negative of n x grad_gamma(PHI) used in the stream function
    %   formulation of A. Reusken, Stream function formulation of
    %   surface Stokes equations, preprint RWTH Aachen, 2018.

    % Copyright (c) 2026 by the fem4m authors

    properties
        firstFaceIndices   % basis functions of the first face
        firstEdgeIndices   % corner and first edge basis functions
        normal4f           % face normals of the mesh
    end

    methods

        function obj = DGPkFace(degree,mesh)
            %DGPkFace - Create a surface Pk element from the first-face basis
            obj@Pk3d(degree,mesh);

            obj.nrBasisFunc = (degree+1)*(degree+2)/2;
            edgeindices = [4+(1:6:6*(degree-1));4+(2:6:6*(degree-1));4+(3:6:6*(degree-1))];
            obj.firstFaceIndices = [1:3,edgeindices(:)',...
                4+6*(degree-1)+(1:4:4*sum(1:(degree-2))) ];
            obj.firstEdgeIndices = [1,2,4+(1:6:6*(degree-1))];

            obj = dofManager(obj);
            obj.transformation = getFaceTransformation(obj);
            obj.normal4f = obj.mesh.normal4f;

        end

        function val = getFaceBasis(obj,~,x_ref,parts)
            %getFaceBasis - Evaluate the surface basis functions
            %   VAL = getFaceBasis(OBJ,PTS,X_REF,PARTS) evaluates
            %   the first-face Pk3d basis at the reference point X_REF
            %   (1-by-2) of the reference triangle, replicates it on
            %   both sides of each face (zero on the second side of
            %   boundary faces), and restricts it to the faces PARTS
            %   (numel(PARTS)-by-nrBasisFunc-by-2). PTS is ignored.
            val = obj.refBasisCoeff*obj.monomBasis(x_ref(:,1),x_ref(:,2),0)';
            %only first face (z=0)
            val = val(obj.firstFaceIndices,:);
            val = reshape((val(:)*ones(1,obj.mesh.nrFaces))',obj.mesh.nrFaces,[]);
            val = ones(obj.mesh.nrFaces,1,2).*val;
            val(obj.mesh.e4f(:,2)==0,:,2)=0;
            val = val(parts,:,:);
        end

        function val = getFaceValue(obj,pts,pts_ref,parts,c)
            %getFaceValue - Evaluate the surface FE function
            %   VAL = getFaceValue(OBJ,PTS,PTS_REF,PARTS,C)
            %   evaluates the surface finite element function with
            %   coefficients C on the faces PARTS at the reference
            %   point PTS_REF (1-by-2), one column per side of the
            %   face (numel(PARTS)-by-2).
            dof4e = obj.dof4f(parts,:);
            coeff4e = c(dof4e);
            basis4e = obj.getFaceBasis(pts,pts_ref,parts);
            val = permute(sum(basis4e.*coeff4e,2),[1 3 2]);
        end

        function val = getFaceGradBasis(obj,~,x_ref,parts)
            %getFaceGradBasis - Evaluate surface gradients of the basis
            %   VAL = getFaceGradBasis(OBJ,PTS,X_REF,PARTS)
            %   evaluates the surface gradients of the first-face
            %   Pk3d basis at the reference point X_REF (1-by-2),
            %   mapped to the faces and projected onto the tangent
            %   plane, on both sides of the faces PARTS
            %   (numel(PARTS)-by-nrBasisFunc-by-3-by-2, zero on the
            %   second side of boundary faces). PTS is ignored.
            m = obj.monomBasis(x_ref(:,1),x_ref(:,2),0).';
            grad = reshape(obj.gradBasisCoefficients*m,size(obj.refBasisCoeff,1),3);
            %only first face (z=0)
            grad = grad(obj.firstFaceIndices,:);
            grad = bsxfun(@times,ones(obj.mesh.nrFaces,1),permute(grad,[3 1 2]));
            ct = obj.transformation;
            val = matMul(grad,ct);
            %projection on tangent plane (P = I - n*n^T)
            normals = obj.normal4f;
            P = -matMul(normals,reshape(normals,size(normals,1),1,3));
            P(:,1,1) = P(:,1,1) + 1;
            P(:,2,2) = P(:,2,2) + 1;
            P(:,3,3) = P(:,3,3) + 1;
            val = permute(matMul(P,permute(val,[1 3 2])),[1 3 2]);
            val = ones(obj.mesh.nrFaces,1,1,2).*val;
            val(obj.mesh.e4f(:,2)==0,:,:,2)=0;
            val = val(parts,:,:,:);
        end

        function val = getFaceGradValue(obj,pts,pts_ref,parts,c)
            %getFaceGradValue - Evaluate surface gradients of the function
            %   VAL = getFaceGradValue(OBJ,PTS,PTS_REF,PARTS,C)
            %   evaluates the surface gradient of the finite
            %   element function with coefficients C on the faces
            %   PARTS at the reference point PTS_REF (1-by-2), on
            %   both sides of the faces (numel(PARTS)-by-3-by-2).
            dof4e = obj.dof4f(parts,:);
            coeff4e = c(dof4e);
            basisGrad4e = obj.getFaceGradBasis(pts,pts_ref,parts);
            val = permute(sum(basisGrad4e.*coeff4e,2),[1 3 4 2]);
        end

        function val = getFaceCurlBasis(obj,~,x_ref,parts)
            %getFaceCurlBasis - Evaluate surface curls of the basis
            %   VAL = getFaceCurlBasis(OBJ,PTS,X_REF,PARTS)
            %   evaluates the surface curl grad(PHI) x n =
            %   grad_gamma(PHI) x n of the first-face Pk3d basis
            %   at the reference point X_REF (1-by-2) on both sides
            %   of the faces PARTS (numel(PARTS)-by-nrBasisFunc-by-3-by-2,
            %   zero on the second side of boundary faces). PTS is
            %   ignored.
            grad = obj.getFaceGradBasis([],x_ref,parts);
            grad = grad(:,:,:,1);
            val = zeros(size(grad));
            normals = obj.normal4f(parts,:);
            val(:,:,1) = grad(:,:,2).*normals(:,3) - grad(:,:,3).*normals(:,2);
            val(:,:,2) = grad(:,:,3).*normals(:,1) - grad(:,:,1).*normals(:,3);
            val(:,:,3) = grad(:,:,1).*normals(:,2) - grad(:,:,2).*normals(:,1);
            val = ones(1,1,1,2).*val;
            val(obj.mesh.e4f(parts,2)==0,:,:,2)=0;
        end

        function val = getFaceCurlValue(obj,pts,pts_ref,parts,c)
            %getFaceCurlValue - Evaluate surface curls of the function
            %   VAL = getFaceCurlValue(OBJ,PTS,PTS_REF,PARTS,C)
            %   evaluates the surface curl grad(PHI) x n of the
            %   finite element function with coefficients C on
            %   the faces PARTS at the reference point PTS_REF
            %   (1-by-2), on both sides of the faces
            %   (numel(PARTS)-by-3-by-2).
            dof4e = obj.dof4f(parts,:);
            coeff4e = c(dof4e);
            basisCurl4e = obj.getFaceCurlBasis(pts,pts_ref,parts);
            val = permute(sum(basisCurl4e.*coeff4e,2),[1 3 4 2]);
        end

        function val = getFaceTransformation(obj)
            %getFaceTransformation - Build the face transformation
            %   VAL = getFaceTransformation(OBJ) computes for every
            %   face the inverse Jacobian of the affine map from the
            %   reference tetrahedron that maps the reference
            %   triangle onto the face (nodes in the order of n4f)
            %   and the fourth vertex to the face centroid moved by
            %   sqrt(area4f) along normal4f (nrFaces-by-3-by-3).
            mesh = obj.mesh;
            P1 = mesh.c4n(mesh.n4f(:,1),:);
            P2 = mesh.c4n(mesh.n4f(:,2),:);
            P3 = mesh.c4n(mesh.n4f(:,3),:);
            P4 = mesh.point4f([1/3,1/3,1/3])+mesh.normal4f.*(sqrt(mesh.area4f)*ones(1,3));
            x2 = P2(:,1)-P1(:,1); y2 = P2(:,2)-P1(:,2); z2 = P2(:,3)-P1(:,3);
            x3 = P3(:,1)-P1(:,1); y3 = P3(:,2)-P1(:,2); z3 = P3(:,3)-P1(:,3);
            x4 = P4(:,1)-P1(:,1); y4 = P4(:,2)-P1(:,2); z4 = P4(:,3)-P1(:,3);

            %compute matrix inverse
            d = x2.*y3.*z4 - x2.*z3.*y4 - y2.*x3.*z4 + y2.*z3.*x4 + z2.*x3.*y4 - z2.*y3.*x4;
            % [c11,c21,c31,c12,c22,c32,c13,c23,c33]
            val = [(y3.*z4-y4.*z3)./d,(-y2.*z4+y4.*z2)./d,(y2.*z3-y3.*z2)./d,...
                (-x3.*z4+x4.*z3)./d,(x2.*z4-x4.*z2)./d,(-x2.*z3+x3.*z2)./d,...
                (x3.*y4-x4.*y3)./d,(-x2.*y4+x4.*y2)./d,(x2.*y3-x3.*y2)./d];
            val = reshape(val,mesh.nrFaces,3,3);
        end
    end

    methods (Access = protected)

        function obj = dofManager(obj,mesh)
            %dofManager - Assemble the degrees of freedom
            %   OBJ = dofManager(OBJ) numbers the degrees of
            %   freedom face by face (dof4f, nrDoF), empties dof4e,
            %   and computes the coordinates of the dof points.
            %
            %   OBJ = dofManager(OBJ,MESH), as called by the Pk3d
            %   constructor, assembles the Pk3d degrees of freedom.
            if nargin > 1
                obj = dofManager@Pk3d(obj,mesh);
                return
            end
            mesh = obj.mesh;
            obj.nrDoF = obj.nrBasisFunc*mesh.nrFaces;
            obj.dof4f = reshape(1:obj.nrDoF,obj.nrBasisFunc,mesh.nrFaces)';
            obj.dof4e = [];
            obj.c4refDof = obj.c4refDof(obj.firstFaceIndices,:);
            obj.c4dof = zeros(obj.nrDoF,3);
            % nodal
            for k = 1 : 3
                obj.c4dof(obj.dof4f(:,k),:) = mesh.c4n(mesh.n4f(:,k),:);
            end
            % edge
            for k = 1 : obj.nrEdgePoints
                for j = 1 : 3
                    obj.c4dof(obj.dof4f(:,3+(k-1)*3+j),:) = ...
                        obj.local2global(obj.mesh.n4f,obj.c4refDof(3+(k-1)*3+j,:),mesh);
                end
            end
            % face points
            for k = 1 : obj.nrFacePoints
                obj.c4dof(obj.dof4f(:,3+3*obj.nrEdgePoints+k),:) = ...
                    obj.local2global(obj.mesh.n4f,obj.c4refDof(3+3*obj.nrEdgePoints+k,:),mesh);
            end
        end

    end
end