classdef FormType
    %FormType - Value types of basis functions and coefficients in forms
    %   FormType names the kind of value an array of the form classes
    %   holds, e.g. FormType.vector2d for a vector field in 2d or
    %   FormType.symmatrix3d for a symmetric 3-by-3 matrix field stored
    %   by its 6 independent components (11,22,33,12,13,23). The members
    %   carry the data needed to work with the stored components, so
    %   that the form classes need not decide this from array sizes.
    %
    %   TYPE = FormType.OF(N,DIM) returns the type of values with N
    %   stored components in space dimension DIM.
    %
    %   Members compare with their names, e.g. TYPE == 'vector2d', and
    %   can be used in switch statements with case 'vector2d' or
    %   case {'matrix2d','matrix3d'}.
    %
    %   The enumeration members of FormType are:
    %       scalar      - scalar values (1 component)
    %       vector2d    - vectors in 2d (2 components)
    %       vector3d    - vectors in 3d (3 components)
    %       matrix2d    - 2-by-2 matrices (4 components)
    %       matrix3d    - 3-by-3 matrices (9 components)
    %       symmatrix2d - symmetric 2-by-2 matrices (3 components)
    %       symmatrix3d - symmetric 3-by-3 matrices (6 components)
    %       tensor      - anything else, e.g. gradients of matrix fields
    %
    %   FormType functions:
    %       of          - Type of values with N components in dimension DIM
    %       vector      - Vector type of the space dimension DIM
    %       isScalar    - True for scalar
    %       isVector    - True for vector2d and vector3d
    %       isMatrix    - True for full and symmetric matrices
    %       isSymmetric - True for symmatrix2d and symmatrix3d
    %       toMatrix    - Reshape stored components into DIM-by-DIM blocks
    %       symPart     - Symmetric part in the stored symmetric components
    %
    %   FormType properties:
    %       kind         - 'scalar', 'vector', 'matrix', 'symmatrix' or
    %                      'tensor'
    %       dim          - space dimension (0 for scalar and tensor)
    %       nrComponents - number of stored components (NaN for tensor)
    %       full         - stored component of each entry of the full
    %                      matrix (column-major)
    %       diagonal     - stored components of the diagonal entries
    %
    %   See also FormEval, FormBasis, FormValue

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = immutable)
        kind            % 'scalar', 'vector', 'matrix', 'symmatrix' or 'tensor'
        dim             % space dimension (0 for scalar and tensor)
        nrComponents    % number of stored components (NaN for tensor)
        full            % stored component of each entry of the full matrix (column-major)
        diagonal        % stored components of the diagonal entries
    end

    enumeration
        scalar      ('scalar',    0, 1,   1,                   1)
        vector2d    ('vector',    2, 2,   1:2,                 [])
        vector3d    ('vector',    3, 3,   1:3,                 [])
        matrix2d    ('matrix',    2, 4,   1:4,                 [1 4])
        matrix3d    ('matrix',    3, 9,   1:9,                 [1 5 9])
        symmatrix2d ('symmatrix', 2, 3,   [1 3 3 2],           [1 2])
        symmatrix3d ('symmatrix', 3, 6,   [1 4 5 4 2 6 5 6 3], [1 2 3])
        tensor      ('tensor',    0, NaN, [],                  [])
    end

    methods

        function obj = FormType(kind,dim,nrComponents,full,diagonal)
            %FormType - Create a member with its component data
            obj.kind = kind;
            obj.dim = dim;
            obj.nrComponents = nrComponents;
            obj.full = full;
            obj.diagonal = diagonal;
        end

        function tf = isScalar(obj)
            %isScalar - True for scalar
            %   TF = isScalar(TYPE) returns true if TYPE is
            %   FormType.scalar.
            tf = strcmp(obj.kind,'scalar');
        end

        function tf = isVector(obj)
            %isVector - True for vector2d and vector3d
            %   TF = isVector(TYPE) returns true if TYPE is
            %   FormType.vector2d or FormType.vector3d.
            tf = strcmp(obj.kind,'vector');
        end

        function tf = isMatrix(obj)
            %isMatrix - True for full and symmetric matrices
            %   TF = isMatrix(TYPE) returns true if TYPE is one of
            %   matrix2d, matrix3d, symmatrix2d or symmatrix3d.
            tf = strcmp(obj.kind,'matrix') || strcmp(obj.kind,'symmatrix');
        end

        function tf = isSymmetric(obj)
            %isSymmetric - True for symmatrix2d and symmatrix3d
            %   TF = isSymmetric(TYPE) returns true if TYPE is
            %   FormType.symmatrix2d or FormType.symmatrix3d.
            tf = strcmp(obj.kind,'symmatrix');
        end

        function val = toMatrix(obj,val,compDim)
            %toMatrix - Reshape stored components into DIM-by-DIM blocks
            %   VAL = toMatrix(TYPE,VAL,COMPDIM) expands the stored
            %   matrix components of VAL along dimension COMPDIM to all
            %   entries of the full matrix and reshapes them into the
            %   two dimensions COMPDIM and COMPDIM+1 of size TYPE.dim.
            %   The entry (i,j) is the component i+DIM*(j-1) of the full
            %   matrix. A TYPE that is no matrix type raises the error
            %   fem4m:FormType:notMatrix.
            if ~obj.isMatrix()
                error('fem4m:FormType:notMatrix', ...
                    'FormType: toMatrix needs a matrix type, not ''%s''',char(obj));
            end
            sz = size(val);
            sz(end+1:compDim) = 1;
            idx = repmat({':'},1,numel(sz));
            idx{compDim} = obj.full;
            val = reshape(val(idx{:}),[sz(1:compDim-1),obj.dim,obj.dim,sz(compDim+1:end)]);
        end

        function val = symPart(obj,val,compDim)
            %symPart - Symmetric part in the stored symmetric components
            %   VAL = symPart(TYPE,VAL,COMPDIM) returns the symmetric parts
            %   (M+M')/2 of the matrices M of TYPE stored along dimension
            %   COMPDIM of VAL as the components of the symmetric type of
            %   the same dimension, [11 22 12] in 2d and [11 22 33 12 13
            %   23] in 3d. Symmetric matrices are returned unchanged; a
            %   TYPE that is no matrix type raises the error
            %   fem4m:FormType:notMatrix.
            if obj.isSymmetric()
                return
            elseif ~obj.isMatrix()
                error('fem4m:FormType:notMatrix', ...
                    'FormType: symPart needs a matrix type, not ''%s''',char(obj));
            end
            sz = size(val);
            sz(end+1:compDim) = 1;
            M = obj.toMatrix(val,compDim);
            order = 1:ndims(M);
            order([compDim compDim+1]) = [compDim+1 compDim];
            M = reshape((M + permute(M,order))/2,[sz(1:compDim-1),obj.dim^2,sz(compDim+1:end)]);
            % first matrix entry of every stored symmetric component
            [~,entries] = unique(FormType.of(obj.dim*(obj.dim+1)/2,obj.dim).full);
            idx = repmat({':'},1,numel(sz));
            idx{compDim} = entries;
            val = M(idx{:});
        end

    end

    methods (Static)

        function type = vector(dim)
            %VECTOR - Vector type of the space dimension DIM
            %   TYPE = FormType.VECTOR(DIM) returns FormType.of(DIM,DIM),
            %   i.e. vector2d for DIM = 2 and vector3d for DIM = 3, e.g.
            %   the type of the columns of a matrix type.
            type = FormType.of(dim,dim);
        end

        function type = of(nrComponents,dim)
            %OF - Type of values with N components in dimension DIM
            %   TYPE = FormType.OF(N,DIM) returns the type of values
            %   with N stored components in space dimension DIM, e.g.
            %   symmatrix2d for N = 3 and DIM = 2 but vector3d for
            %   N = 3 and DIM = 3. N may also be the sizes of all
            %   component dimensions of an array; more than one
            %   nonsingleton component dimension gives tensor, as does
            %   any N that matches no other type.
            if ~isscalar(nrComponents)
                if nnz(nrComponents > 1) > 1
                    type = FormType.tensor;
                    return
                end
                nrComponents = prod(nrComponents);
            end
            type = FormType.tensor;
            if nrComponents == 1
                type = FormType.scalar;
            elseif dim == 2
                switch nrComponents
                    case 2
                        type = FormType.vector2d;
                    case 3
                        type = FormType.symmatrix2d;
                    case 4
                        type = FormType.matrix2d;
                end
            elseif dim == 3
                switch nrComponents
                    case 3
                        type = FormType.vector3d;
                    case 6
                        type = FormType.symmatrix3d;
                    case 9
                        type = FormType.matrix3d;
                end
            end
        end

    end
end