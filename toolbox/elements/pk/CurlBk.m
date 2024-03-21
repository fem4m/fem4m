classdef CurlBk < Bk
    %CurlBk - Curls of the cubic bubble functions of Bk
    %   OBJ = CurlBk(DEGREE,MESH) creates the divergence-free
    %   vector-valued element space of polynomial degree DEGREE >= 2 (a
    %   scalar or one value per element) on the triangular mesh MESH of
    %   the curls curl(phi) = (d_y phi, -d_x phi) of the cubic bubble
    %   basis functions phi of Bk(DEGREE+1,MESH).
    %
    %   CurlBk functions:
    %       getBasis     - Evaluate the basis functions
    %       getDivBasis  - Evaluate the divergence of the basis
    %       getGradBasis - Evaluate the gradients of the basis
    %
    %   See also Bk, CurlBRotk, Nd1k

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = CurlBk(degree,mesh)
            %CurlBk - Create a curl bubble finite element object
            if any(degree(:) < 2)
                error('fem4m:CurlBk:invalidDegree','CurlBk degree must be at least 2');
            end
            obj = obj@Bk(degree+1,mesh);
            obj.degree4e = obj.degree4e-1;
            obj.degree = obj.degree-1;
            obj.degree4DoF = obj.degree4DoF-1;
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,X,PARTS) evaluates the curls
            %   (d_y phi, -d_x phi) of the cubic bubble basis phi at the
            %   reference point X (1-by-2) on the physical elements.
            %   VAL is nrElems-by-nrBasisFunc-by-2. PTS and PARTS are
            %   ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,2);
            grad = permute(Bk.getGradElementBasis(obj.degree+1,x),[3 1 2]);
            c = obj.mesh.transformation;
            grad = matMul(grad,permute(c,[1 3 2]));
            % curl
            val(:,:,1) = grad(:,:,2);
            val(:,:,2) = -grad(:,:,1);
        end

        function val = getDivBasis(obj,~,~,~)
            %getDivBasis - Evaluate the divergence of the basis
            %   VAL = getDivBasis(OBJ,PTS,X,PARTS) returns the zero
            %   divergences of the divergence-free basis functions as
            %   an nrElems-by-nrBasisFunc array. PTS, X, and PARTS are
            %   ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc);
        end

        function val = getGradBasis(obj,~,x,~)
            %getGradBasis - Evaluate the gradients of the basis
            %   VAL = getGradBasis(OBJ,PTS,X,PARTS) evaluates the
            %   gradients of the curl bubble basis functions v at the
            %   reference point X (1-by-2), mapped to the physical
            %   elements. VAL is nrElems-by-nrBasisFunc-by-4 with the
            %   Jacobian entries [d_x v_1,d_x v_2,d_y v_1,d_y v_2]
            %   (column by column). PTS and PARTS are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,2,2);
            d2 = permute(Bk.getD2ElementBasis(obj.degree+1,x),[3 1 2]);
            c = obj.mesh.transformation(x);
            c11 = c(:,1,1);
            c12 = c(:,1,2);
            c21 = c(:,2,1);
            c22 = c(:,2,2);
            vald2 = zeros(obj.mesh.nrElems,obj.nrBasisFunc,3);
            for i = 1 : size(d2,2)
                z1 = d2(:,i,1);
                z2 = d2(:,i,2);
                z3 = d2(:,i,3);
                vald2(:,i,1) = z1.*c11.*c11+z2.*c12.*c12+2*z3.*c11.*c12;
                vald2(:,i,2) = z1.*c21.*c21+z2.*c22.*c22+2*z3.*c21.*c22;
                vald2(:,i,3) = z1.*c11.*c21+z2.*c12.*c22+  z3.*(c11.*c22+c21.*c12);
            end
            val(:,:,1,1) = vald2(:,:,3);
            val(:,:,1,2) = vald2(:,:,2);
            val(:,:,2,1) = -vald2(:,:,1);
            val(:,:,2,2) = -vald2(:,:,3);
            % flatten the Jacobian column by column
            val = val(:,:,:);
        end

    end

end