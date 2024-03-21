classdef FESkewSymMatrix < FEVector
    %FESkewSymMatrix - Skew-symmetric matrix-valued finite element object
    %   OBJ = FESkewSymMatrix(FE1) builds the skew-symmetric matrix field
    %   [0 -p; p 0] on a 2d mesh from the functions p of the scalar
    %   element FE1, stored column by column as [11 21 12 22] =
    %   [0 p -p 0]. Every basis function is divided by the square root
    %   of the area of the support of its degree of freedom (for DGDk
    %   the element area).
    %
    %   FESkewSymMatrix functions:
    %       getBasis     - Evaluate the skew-symmetric matrix basis
    %       prolongation - Prolongation matrix of the scaled basis
    %
    %   See also FEVector, FEMatrix, FESymMatrix

    % Copyright (c) 2026 by the fem4m authors

    methods

        function val = getBasis(obj,pts,pts_ref,parts)
            %getBasis - Evaluate the skew-symmetric matrix basis
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the basis
            %   functions p of the scalar element at the reference point
            %   PTS_REF (1-by-2) in the elements PARTS as the matrix fields
            %   [0 -p; p 0], stored as [0 p -p 0], each divided by the
            %   square root of the area of the support of its degree of
            %   freedom.
            basis = getBasis@FEVector(obj,pts,pts_ref,parts);
            val = zeros(size(basis,1),size(basis,2),4);
            val(:,:,2) = basis;
            val(:,:,3) = -basis;
            scaling = dofScaling(obj);
            val = val./sqrt(scaling(obj.dof4e));
        end

        function P = prolongation(obj,feOld)
            %PROLONGATION - Prolongation matrix of the scaled basis
            %   P = PROLONGATION(OBJ,FEOLD) returns the matrix that maps
            %   the coefficients of a function of FEOLD (an
            %   FESkewSymMatrix, its scalar element or a cell array of it)
            %   to the coefficients of the same function in OBJ. The
            %   component prolongation is rescaled by the mesh-dependent
            %   normalization of the basis.
            if ~isa(feOld,'FEVector')
                if ~iscell(feOld)
                    feOld = {feOld};
                end
                feOld = FESkewSymMatrix(feOld{:});
            end
            P = prolongation@FEVector(obj,feOld);
            sF = sqrt(dofScaling(obj));
            sC = 1./sqrt(dofScaling(feOld));
            sC(~isfinite(sC)) = 0;
            P = spdiags(sF,0,numel(sF),numel(sF))*P*spdiags(sC,0,numel(sC),numel(sC));
        end

    end

    methods (Access = protected)

        function scaling = dofScaling(obj)
            %dofScaling - Area of the support of every degree of freedom
            area4e = obj.mesh.area4e.*ones(1,size(obj.dof4e,2));
            scaling = accumarray(obj.dof4e(:),area4e(:),[obj.nrDoF 1]);
        end

        function tf = isMatrixValued(~)
            %isMatrixValued - True for matrix-valued elements
            tf = true;
        end

    end
end