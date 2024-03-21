classdef FESymMatrix < FEVector
    %FESymMatrix - Symmetric matrix-valued finite element object
    %   OBJ = FESymMatrix(FE11,FE22,FE12) builds a symmetric 2-by-2
    %   matrix field on a 2d mesh from three scalar elements, one per
    %   stored component [11 22 12].
    %
    %   OBJ = FESymMatrix(FE11,FE22,FE33,FE12,FE13,FE23) builds a
    %   symmetric 3-by-3 matrix field on a 3d mesh from six scalar
    %   elements, one per stored component [11 22 33 12 13 23].
    %
    %   Besides the routines of FEVector, FESymMatrix evaluates the row
    %   divergence, the double divergence and the trace (2d and 3d) and
    %   the row-wise rot and rot(rot) (2d only).
    %
    %   FESymMatrix functions:
    %       getDivBasis      - Evaluate the divergence of the basis
    %       getDivDivBasis   - Evaluate double divergences of the basis
    %       getRotBasis      - Evaluate the row-wise rot of the basis
    %       getRotRotBasis   - Evaluate the rot of the rot of the basis
    %       getDivValue      - Evaluate the row divergences of the function
    %       getTrValue       - Evaluate the trace of the function
    %
    %   See also FEVector, FEMatrix, FESkewSymMatrix, FormType

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = getDivBasis(obj,pts,pts_ref,parts)
            %getDivBasis - Evaluate the divergence of the basis functions
            %   VAL = getDivBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the row
            %   divergences sum_j d_j s_ij of the symmetric matrix basis
            %   functions s at the reference point PTS_REF (1-by-d) in the
            %   elements PARTS, one component per row of s.
            grad = obj.getGradBasis(pts,pts_ref,parts); % (:,:,component,[x y (z)])
            type = obj.symType();
            idx = reshape(type.full,type.dim,type.dim);
            val = zeros(size(grad,1),size(grad,2),type.dim);
            for i = 1 : type.dim
                for j = 1 : type.dim
                    val(:,:,i) = val(:,:,i) + grad(:,:,idx(i,j),j);
                end
            end
        end

        function val = getDivDivBasis(obj,pts,pts_ref,parts)
            %getDivDivBasis - Evaluate double divergences of the basis
            %   VAL = getDivDivBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   double divergence sum_ij d_i d_j s_ij of the symmetric
            %   matrix basis functions s at the reference point PTS_REF
            %   (1-by-d) in the elements PARTS. The second derivatives
            %   [xx yy xy] in 2d and [xx yy zz xy xz yz] in 3d are ordered
            %   like the stored components, so the double divergence is the
            %   sum over the stored components k of their k-th second
            %   derivative, the off-diagonal ones counted twice.
            d2 = obj.getD2Basis(pts,pts_ref,parts); % (:,:,component,second derivative)
            type = obj.symType();
            val = zeros(size(d2,1),size(d2,2));
            for k = type.diagonal
                val = val + d2(:,:,k,k);
            end
            for k = setdiff(1:type.nrComponents,type.diagonal)
                val = val + 2*d2(:,:,k,k);
            end
        end

        function val = getRotBasis(obj,pts,pts_ref,parts)
            %getRotBasis - Evaluate the row-wise rot of the basis functions
            %   VAL = getRotBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   row-wise rotations (d_x s12 - d_y s11, d_x s22 - d_y s12)
            %   of the 2d symmetric matrix basis functions s at the
            %   reference point PTS_REF (1-by-2) in the elements PARTS. On
            %   a 3d mesh it raises the error fem4m:FE:rotNotDefined.
            if obj.symType() ~= FormType.symmatrix2d
                obj.notDefined('rot');
            end
            grad = obj.getGradBasis(pts,pts_ref,parts); % (:,:,[11 22 12],[x y])
            val = cat(3,grad(:,:,3,1) - grad(:,:,1,2),grad(:,:,2,1) - grad(:,:,3,2));
        end

        function val = getRotRotBasis(obj,pts,pts_ref,parts)
            %getRotRotBasis - Evaluate the rot of the rot of the basis
            %   VAL = getRotRotBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   scalar rot(rot(s)) = d_yy s11 - 2 d_xy s12 + d_xx s22 of
            %   the 2d symmetric matrix basis functions s at the reference
            %   point PTS_REF (1-by-2) in the elements PARTS. On a 3d mesh
            %   it raises the error fem4m:FE:rotRotNotDefined.
            if obj.symType() ~= FormType.symmatrix2d
                obj.notDefined('rotRot');
            end
            d2 = obj.getD2Basis(pts,pts_ref,parts); % (:,:,[11 22 12],[xx yy xy])
            val = d2(:,:,1,2) + d2(:,:,2,1) - 2*d2(:,:,3,3);
        end

        function val = getDivValue(obj,pts,pts_ref,parts,c)
            %getDivValue - Evaluate the row divergences of the matrix function
            %   VAL = getDivValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the row
            %   divergences sum_j d_j s_ij of the symmetric matrix-valued
            %   finite element function s with coefficients C at the
            %   reference point PTS_REF (1-by-d) in the elements PARTS, one
            %   column per row of s.
            grad = obj.getGradValue(pts,pts_ref,parts,c); % (:,component,[x y (z)])
            type = obj.symType();
            idx = reshape(type.full,type.dim,type.dim);
            val = zeros(size(grad,1),type.dim);
            for i = 1 : type.dim
                for j = 1 : type.dim
                    val(:,i) = val(:,i) + grad(:,idx(i,j),j);
                end
            end
        end

        function val = getTrValue(obj,pts,pts_ref,parts,c)
            %getTrValue - Evaluate the trace of the matrix function
            %   VAL = getTrValue(OBJ,PTS,PTS_REF,PARTS,C) evaluates the trace
            %   of the symmetric matrix-valued finite element function with
            %   coefficients C at the reference point PTS_REF (1-by-d) in
            %   the elements PARTS, the sum of its diagonal components.
            val = obj.getValue(pts,pts_ref,parts,c);
            val = sum(val(:,obj.symType().diagonal),2);
        end
    end

    methods (Access = protected)

        function tf = isMatrixValued(~)
            %isMatrixValued - True for matrix-valued elements
            tf = true;
        end

        function type = symType(obj)
            %symType - FormType of the stored components
            %   TYPE = symType(OBJ) returns FormType.symmatrix2d or
            %   FormType.symmatrix3d, whose full and diagonal give the
            %   stored component of every matrix entry and of the
            %   diagonal entries.
            type = FormType.of(size(obj.fe,1),size(obj.mesh.c4n,2));
        end

    end
end