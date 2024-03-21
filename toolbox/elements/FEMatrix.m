classdef FEMatrix < FEVector
    %FEMatrix - Matrix-valued finite element built from row or scalar elements
    %   OBJ = FEMatrix(ROW1,ROW2) builds a 2-by-2 matrix field on a 2d
    %   mesh from two vector-valued elements, one per row of the
    %   matrix, e.g. FEMatrix(BDMk(DEGREE,MESH),BDMk(DEGREE,MESH)) for
    %   stresses with H(div)-conforming rows.
    %
    %   OBJ = FEMatrix(FE11,FE21,FE12,FE22) builds a matrix field from
    %   one scalar element per entry, ordered column by column; in 3d it
    %   takes the nine elements FE11,FE21,FE31,FE12,...,FE33 (see
    %   DGPkMatrix).
    %
    %   Both store the values column by column, [11 21 12 22] in 2d, and
    %   manage the combined degrees of freedom as FEVector does. On a 2d
    %   mesh with a symmetry boundary (MESH.n4Sb), FEMatrix with two row
    %   elements that have edge degrees of freedom (dof4ed, e.g. RTk or
    %   BDMk) also fixes the edge degrees of freedom of ROW1 on the
    %   symmetry edges with fixed displacement component 2 and those of
    %   ROW2 on the edges with fixed component 1 (third column of n4Sb).
    %
    %   OBJ = FEMatrix(ROW1,ROW2) + CurlBRotk(DEGREE,MESH) appends a
    %   CurlBRotk enrichment, stored in the property CurlBRotk (see
    %   plus); FEMatrix(ROW1,ROW2) + FEMatrix(CurlBk(...),CurlBk(...))
    %   builds the Stenberg space.
    %
    %   FEMatrix functions:
    %       getBasis          - Evaluate the basis functions
    %       getDivBasis       - Evaluate the divergence of the basis
    %       getRotBasis       - Evaluate the row-wise rot of the basis
    %       getDivDivBasis    - Evaluate double divergences of the basis
    %       projectNbBoundary - Project data onto the Neumann boundary
    %       prolongation      - Prolongation matrix from a coarser space
    %       plus              - Merge two matrix-valued elements
    %
    %   FEMatrix properties:
    %       CurlBRotk - Optional CurlBRotk element appended to the basis
    %
    %   See also FEVector, FESymMatrix, FESkewSymMatrix, CurlBRotk

    % Copyright (c) 2026 by the fem4m authors

    properties
        CurlBRotk    % optional CurlBRotk element appended to the basis
    end

    methods

        function obj = FEMatrix(varargin)
            %FEMatrix - Create a matrix-valued finite element object
            obj = obj@FEVector(varargin{:});
            if isprop(obj.mesh,'n4Sb') && ~isempty(obj.mesh.n4Sb) && size(obj.mesh.c4n,2) == 2 ...
                    && numel(obj.fe) == 2 && isprop(obj.fe{1},'dof4ed')
                fixedNodesX = unique(obj.fe{1}.dof4ed(obj.mesh.ed4Sb(obj.mesh.n4Sb(:,3)==2),:));
                fixedNodesY = unique(obj.fe{1}.nrDoF + obj.fe{2}.dof4ed(obj.mesh.ed4Sb(obj.mesh.n4Sb(:,3)==1),:));
                obj.fixedNodes = [obj.fixedNodes(:);fixedNodesX(:);fixedNodesY(:)];
                obj.freeNodes = setdiff(1:obj.nrDoF,obj.fixedNodes)';
            end
        end

        function val = getBasis(obj,pts,pts_ref,parts)
            %getBasis - Evaluate the basis functions
            %   VAL = getBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the matrix
            %   basis functions at the reference point PTS_REF (1-by-d) in
            %   the elements PARTS, with the entries column by column along
            %   the third dimension ([11 21 12 22] in 2d), and appends the
            %   CurlBRotk basis if present.
            val = getBasis@FEVector(obj,pts,pts_ref,parts);
            if size(val,4) > 1
                val = reshape(val,size(val,1),size(val,2),[]);
            end
            if ~isempty(obj.CurlBRotk)
                val(:,end+(1:obj.CurlBRotk.nrBasisFunc),:) = obj.CurlBRotk.getBasis(pts,pts_ref,parts);
            end
        end

        function val = getDivBasis(obj,pts,pts_ref,parts)
            %getDivBasis - Evaluate the divergence of the basis functions
            %   VAL = getDivBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the row
            %   divergences of the matrix basis functions at the reference
            %   point PTS_REF (1-by-d) in the elements PARTS, VAL(:,J,I) the
            %   divergence of row I of the J-th basis function, from the
            %   divergences of the row elements or from the gradients of
            %   the scalar entries, and of the CurlBRotk basis if present.
            if isempty(parts) || obj.isRefElemBased
                nrRows = obj.mesh.nrElems;
            else
                nrRows = numel(parts);
            end
            val = zeros(nrRows,obj.nrBasisFunc,2);
            if size(obj.fe,1) == 2
                val(:,1:obj.fe{1}.nrBasisFunc,1) = obj.fe{1}.getDivBasis(pts,pts_ref,parts);
                val(:,obj.fe{1}.nrBasisFunc+(1:obj.fe{2}.nrBasisFunc),2) = obj.fe{2}.getDivBasis(pts,pts_ref,parts);
                indx = obj.fe{1}.nrBasisFunc+obj.fe{2}.nrBasisFunc;
                for curBFE = 2 : size(obj.fe,2)
                    val(:,indx+(1:obj.fe{1,curBFE}.nrBasisFunc),1) = obj.fe{1,curBFE}.getDivBasis(pts,pts_ref,parts);
                    val(:,indx+obj.fe{1,curBFE}.nrBasisFunc+(1:obj.fe{2,curBFE}.nrBasisFunc),2) = obj.fe{2,curBFE}.getDivBasis(pts,pts_ref,parts);
                    indx = indx + obj.fe{1,curBFE}.nrBasisFunc+obj.fe{2,curBFE}.nrBasisFunc;
                end
            else
                % components stored column-major, e.g. (:,:,[11 21 12 22],[x y])
                grad = obj.getGradBasis(pts,pts_ref,parts);
                dim = size(grad,4);
                val = zeros(nrRows,obj.nrBasisFunc,dim);
                for row = 1 : dim
                    for col = 1 : dim
                        val(:,1:size(grad,2),row) = val(:,1:size(grad,2),row) + grad(:,:,(col-1)*dim+row,col);
                    end
                end
            end
            if ~isempty(obj.CurlBRotk)
                val(:,end-obj.CurlBRotk.nrBasisFunc+1:end,:) = obj.CurlBRotk.getDivBasis(pts,pts_ref,parts);
            end
        end

        function val = getRotBasis(obj,pts,pts_ref,parts)
            %getRotBasis - Evaluate the row-wise rot of the basis functions
            %   VAL = getRotBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   row-wise rotations (d_x s12 - d_y s11, d_x s22 - d_y s21)
            %   of the 2d matrix basis functions s at the reference point
            %   PTS_REF (1-by-2) in the elements PARTS, from the rotations
            %   of the two row elements or from the gradients of the four
            %   scalar entries. On a 3d mesh, for other layouts and with a
            %   CurlBRotk element it raises an error.
            if isa(obj.mesh,'Mesh3d') || ~any(size(obj.fe,1) == [2 4])
                obj.notDefined('rot');
            elseif ~isempty(obj.CurlBRotk)
                error('fem4m:FEMatrix:rotNotSupported', ...
                    'FEMatrix: rot(u) of the CurlBRotk enrichment is not implemented.');
            end
            if isempty(parts) || obj.isRefElemBased
                nrRows = obj.mesh.nrElems;
            else
                nrRows = numel(parts);
            end
            if size(obj.fe,1) == 2
                val = zeros(nrRows,obj.nrBasisFunc,2);
                indx = 0;
                for curBFE = 1 : size(obj.fe,2)
                    for row = 1 : 2
                        nrBasisFunc = obj.fe{row,curBFE}.nrBasisFunc;
                        val(:,indx+(1:nrBasisFunc),row) = obj.fe{row,curBFE}.getRotBasis(pts,pts_ref,parts);
                        indx = indx + nrBasisFunc;
                    end
                end
            else
                grad = obj.getGradBasis(pts,pts_ref,parts); % (:,:,[11 21 12 22],[x y])
                val = cat(3,grad(:,:,3,1) - grad(:,:,1,2),grad(:,:,4,1) - grad(:,:,2,2));
            end
        end

        function val = getDivDivBasis(obj,pts,pts_ref,parts)
            %getDivDivBasis - Evaluate double divergences of the basis
            %   VAL = getDivDivBasis(OBJ,PTS,PTS_REF,PARTS) evaluates the
            %   double divergence sum_ij d_i d_j s_ij of the matrix basis
            %   functions s at the reference point PTS_REF (1-by-d) in the
            %   elements PARTS, from the second derivatives of the scalar
            %   entries. For row elements and with a CurlBRotk element it
            %   raises an error.
            dim = size(obj.mesh.c4n,2);
            if size(obj.fe,1) ~= dim^2 || ~isempty(obj.CurlBRotk)
                error('fem4m:FEMatrix:divDivNotSupported', ...
                    'FEMatrix: div(div(u)) is only implemented for scalar entry elements.');
            end
            d2 = obj.getD2Basis(pts,pts_ref,parts); % (:,:,[11 21 12 22],[xx yy xy])
            % the second derivative d_i d_j of the entry (i,j), ordered like
            % the stored components of a symmetric matrix
            deriv = FormType.of(dim*(dim+1)/2,dim).full;
            val = zeros(size(d2,1),size(d2,2));
            for k = 1 : dim^2
                val = val + d2(:,:,k,deriv(k));
            end
        end

        function val = projectNbBoundary(obj,g,gDegree)
            %projectNbBoundary - Project data onto the Neumann boundary
            %   VAL = projectNbBoundary(OBJ,G) projects the traction G
            %   row by row onto the Neumann boundary degrees of freedom of
            %   the two row elements (RTk or BDMk) and returns the nrDoF
            %   coefficients. G(X,NORMALS) gets the points X (N-by-2) and
            %   the outer unit normals NORMALS (N-by-2) and returns
            %   sigma*NORMALS as an N-by-2 array, whose column I is the
            %   normal trace of row I. Other layouts raise an error.
            %
            %   VAL = projectNbBoundary(OBJ,G,GDEGREE) also specifies the
            %   polynomial degree GDEGREE of G, which only sets the
            %   quadrature degree of the right-hand side (default: the
            %   degree of the row elements).
            if nargin < 3
                gDegree = [];
            end
            if size(obj.fe,1) ~= 2
                error('fem4m:FEMatrix:neumannNotImplemented','Neumann boundary not implemented');
            end
            val = zeros(obj.nrDoF,1);
            val(1:obj.fe{1}.nrDoF) = obj.fe{1}.projectNbBoundary(@(x,n)(sum(g(x,n).*[1,0],2)),gDegree);
            val(obj.fe{1}.nrDoF+(1:obj.fe{2}.nrDoF)) = obj.fe{2}.projectNbBoundary(@(x,n)(sum(g(x,n).*[0,1],2)),gDegree);
        end

        function P = prolongation(obj,feOld)
            %PROLONGATION - Prolongation matrix from a coarser space
            %   P = PROLONGATION(OBJ,FEOLD) returns the prolongation matrix
            %   of FEVector.prolongation. With a CurlBRotk element OBJ is a
            %   sum, and P is the local L2 projection onto the whole fine
            %   space, exact if the two spaces are nested, as for
            %   BDMk + CurlBRotk under uniform refinement. If OBJ has a
            %   CurlBRotk element, FEOLD must be an FEMatrix with one, and
            %   vice versa.
            hasOld = isa(feOld,'FEMatrix') && ~isempty(feOld.CurlBRotk);
            if isempty(obj.CurlBRotk) == hasOld
                error('fem4m:FEMatrix:curlBRotkMismatch', ...
                    'FEMatrix: either both or none of the spaces must have a CurlBRotk element.');
            end
            P = prolongation@FEVector(obj,feOld);
        end

        function obj = plus(obj1,obj2)
            %PLUS - Merge two matrix-valued finite element objects
            %   OBJ3 = PLUS(OBJ1,OBJ2), or OBJ3 = OBJ1 + OBJ2, merges two
            %   FEMatrix objects as FEVector.plus does, e.g. the Stenberg
            %   space FEMatrix(BDMk,BDMk) + FEMatrix(CurlBk,CurlBk).
            %
            %   If one of the inputs is a CurlBRotk element, it is stored
            %   in OBJ3.CurlBRotk and its basis is appended to that of the
            %   FEMatrix input; its degrees of freedom are shifted, and
            %   those of degree at most the degree of the FEMatrix (one
            %   less for RTk rows) are removed from the free set. Other
            %   inputs raise an error.
            if isprop(obj1,'fe') && isprop(obj2,'fe')
                obj = plus@FEVector(obj1,obj2);
            else
                if isprop(obj2,'fe') && isa(obj1,'CurlBRotk')
                    obj = obj2;
                    obj2 = obj1;
                    obj1 = obj;
                elseif isprop(obj1,'fe') && isa(obj2,'CurlBRotk')
                    obj = obj1;
                else
                    error('fem4m:FEMatrix:plusNotImplemented','not implemented');
                end
                obj.CurlBRotk = obj2;
                obj.degree = max(obj1.degree,obj2.degree);
                obj.dof4e = [obj1.dof4e,obj1.allnrDoF+obj2.dof4e];
                obj.nrDoF = obj1.nrDoF+obj2.nrDoF;
                obj.allnrDoF = obj1.allnrDoF+obj2.allnrDoF;
                obj.freeNodes = [obj1.freeNodes(:);obj1.allnrDoF+obj2.freeNodes(:)];
                if isa(obj1.fe{1},'RTk')
                    linearDependentDoFs = obj1.allnrDoF+obj2.dof4e(:,obj2.degree4DoF<=obj1.degree-1);
                else
                    linearDependentDoFs = obj1.allnrDoF+obj2.dof4e(:,obj2.degree4DoF<=obj1.degree);
                end
                obj.freeNodes = setdiff(obj.freeNodes,linearDependentDoFs);
                obj.fixedNodes = [obj1.fixedNodes(:);obj1.allnrDoF+obj2.fixedNodes(:)];
                obj.nonDoFNodes = [obj1.nonDoFNodes(:);obj1.allnrDoF+obj2.nonDoFNodes(:)];
                obj.nrBasisFunc = obj1.nrBasisFunc+obj2.nrBasisFunc;
            end
        end



    end

    methods (Access = protected)

        function tf = isMatrixValued(~)
            %isMatrixValued - True for matrix-valued elements
            tf = true;
        end

    end
end