classdef CurlBRotk < DGDk
    %CurlBRotk - Matrix bubble element space with curl rows
    %   OBJ = CurlBRotk(DEGREE,MESH) creates the discontinuous
    %   matrix-valued element space of polynomial degree DEGREE >= 2, a
    %   scalar or one degree per element, on the triangular mesh MESH
    %   whose rows are the curls
    %   curl(b_T*d_x q) and curl(b_T*d_y q) of the cubic bubble b_T
    %   times the derivatives of the non-constant Dubiner polynomials q
    %   of degree DEGREE-1, with DEGREE*(DEGREE+1)/2 - 1 basis functions
    %   per element. FEMatrix(ROW1,ROW2) + CurlBRotk(DEGREE,MESH) appends
    %   these functions to a matrix-valued element with two row elements
    %   (FEMatrix.plus stores them in the property CurlBRotk).
    %
    %   CurlBRotk functions:
    %       getBasis    - Evaluate the basis functions
    %       getDivBasis - Evaluate the divergence of the basis
    %
    %   See also Bk, CurlBk, DGDk, FEMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods
        function obj = CurlBRotk(degree,mesh)
            %CurlBRotk - Create a matrix curl bubble finite element object
            if any(degree(:) < 2)
                error('fem4m:CurlBRotk:invalidDegree','CurlBRotk degree must be at least 2');
            end
            obj = obj@DGDk(degree-1,mesh);
            obj.degree4e = obj.degree4e+1;
            obj.degree = obj.degree+1;
            obj.degree4DoF = obj.degree4DoF+1;
            % the constant Dubiner polynomial has no gradient, so its
            % basis function vanishes and is removed
            obj.nrBasisFunc = obj.nrBasisFunc - 1;
            obj.nrBasis4e = obj.nrBasisFunc;
            obj.degree4DoF = obj.degree4DoF(2:end);
            obj.dof4e = reshape(1:mesh.nrElems*obj.nrBasisFunc,obj.nrBasisFunc,mesh.nrElems)';
            obj.allnrDoF = mesh.nrElems*obj.nrBasisFunc;
            % unused dofs of the elements of lower degree in the new
            % numbering
            obj.nonDoFNodes = sort(obj.dof4e(obj.degree4DoF > obj.degree4e));
            obj.nrDoF = obj.allnrDoF - numel(obj.nonDoFNodes);
            obj.freeNodes = setdiff(1:obj.allnrDoF,obj.nonDoFNodes);
            obj.fixedNodes = [];
        end

        function val = getBasis(obj,~,x,~)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,X,PARTS) evaluates the four
            %   components of the curl bubble basis at the
            %   reference point X (1-by-2) and scales them by the
            %   square root of the element area. VAL is
            %   nrElems-by-nrBasisFunc-by-4 with the components
            %   [s11,s21,s12,s22] (column by column) of the matrices
            %   s with the rows curl(b_T*d_x q) and curl(b_T*d_y q).
            %   PTS and PARTS are ignored.
            gradbasis = permute(DGDk.getGradElementBasis(obj.degree-1,x),[3 1 2]);
            % all Dubiner polynomials, the constant one is removed below
            nb = size(gradbasis,2);
            val = ones(obj.mesh.nrElems,nb,4);
            c = obj.mesh.transformation;
            gradbasis = matMul(gradbasis,permute(c,[1 3 2]));
            d2basis = permute(DGDk.getD2ElementBasis(obj.degree-1,x),[3 1 2]);
            c11 = c(:,1,1);
            c12 = c(:,1,2);
            c21 = c(:,2,1);
            c22 = c(:,2,2);
            d2 = zeros(obj.mesh.nrElems,nb,3);
            for i = 1 : size(d2basis,2)
                z1 = d2basis(:,i,1);
                z2 = d2basis(:,i,2);
                z3 = d2basis(:,i,3);
                d2(:,i,1) = z1.*c11.*c11+z2.*c12.*c12+2*z3.*c11.*c12;
                d2(:,i,2) = z1.*c21.*c21+z2.*c22.*c22+2*z3.*c21.*c22;
                d2(:,i,3) = z1.*c11.*c21+z2.*c12.*c22+  z3.*(c11.*c22+c21.*c12);
            end
            d2basis = d2;
            y = x(:,2)';
            x = x(:,1)';
            B(1,:,1) = (1-x-y).*x.*y;
            gradB(1,:,1) = -y.*(x+y-1.0)-x.*y;
            gradB(1,:,2) = -x.*(x+y-1.0)-x.*y;
            gradB = matMul(gradB,permute(c,[1 3 2]));
            val(:,:,1) = gradB(:,:,2).*gradbasis(:,:,1)+B.*d2basis(:,:,3);
            val(:,:,3) = -gradB(:,:,1).*gradbasis(:,:,1)-B.*d2basis(:,:,1);
            val(:,:,2) = gradB(:,:,2).*gradbasis(:,:,2)+B.*d2basis(:,:,2);
            val(:,:,4) = -gradB(:,:,1).*gradbasis(:,:,2)-B.*d2basis(:,:,3);

            val = val(:,2:end,:).*sqrt(obj.mesh.area4e);
        end

        function val = getDivBasis(obj,~,~,~)
            %getDivBasis - Evaluate the divergence of the basis
            %   VAL = getDivBasis(OBJ,PTS,X,PARTS) returns the zero
            %   row-wise divergences of the basis functions as an
            %   nrElems-by-nrBasisFunc-by-2 array. PTS, X, and PARTS
            %   are ignored.
            val = zeros(obj.mesh.nrElems,obj.nrBasisFunc,2);
        end

    end
end