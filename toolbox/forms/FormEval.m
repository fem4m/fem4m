classdef FormEval
    %FormEval - Evaluate coefficients and functions in form strings
    %   FormEval provides static helpers used by the form classes,
    %   FormValue and the elements (FE, FEVector) to turn coefficients
    %   and source functions into array data on quadrature points.
    %   PTS, PTS_REF and PARTS are the arguments that integrate passes
    %   to an integrand: the points (one row per element or side of
    %   PARTS), the reference point (1-by-d, empty on polygonal meshes)
    %   and the element or side indices.
    %
    %   FormEval functions:
    %       evalCoefficient - Evaluate a coefficient as array data
    %       evalFunction    - Evaluate a source function on points
    %
    %   See also FormType, FormBasis, FormValue, Form

    % Copyright (c) 2026 by the fem4m authors

    methods(Static)

        function val = evalCoefficient(pts,pts_ref,parts,alpha,type,target)
            %evalCoefficient - Evaluate a coefficient as array data
            %   VAL = evalCoefficient(PTS,PTS_REF,PARTS,ALPHA,TYPE)
            %   evaluates the coefficient ALPHA that acts on values of
            %   the FormType TYPE: a scalar for FormType.scalar and
            %   FormType.tensor and otherwise an N-by-N matrix stored as
            %   1-by-N-by-N array, N = TYPE.nrComponents, e.g. 3-by-3
            %   for vector3d and symmatrix2d; for an ALPHA per element
            %   or per point VAL has one row per point instead of 1.
            %   Empty ALPHA gives the identity, a scalar ALPHA a multiple of
            %   the identity. A TYPE that is not a FormType raises the error
            %   fem4m:FormEval:invalidType.
            %
            %   VAL = evalCoefficient(PTS,PTS_REF,PARTS,ALPHA,TYPE,TARGET)
            %   with TARGET = FormType.scalar and a vector TYPE evaluates
            %   the vector coefficient that maps values of TYPE to
            %   scalars, e.g. b in b*grad(u)*v, as 1-by-N array (one row
            %   per point for an ALPHA per element or side; the value of a
            %   function handle unchanged) with default ones. TARGET = TYPE
            %   is the first syntax, and TYPE = FormType.tensor ignores
            %   TARGET; any other combination raises the error
            %   fem4m:FormEval:invalidType.
            %
            %   ALPHA may be empty, a scalar, the constant coefficient
            %   with N*N (or N) entries, a column of one scalar per
            %   element or side (indexed by PARTS), or a function handle
            %   of PTS or of PTS, PTS_REF and PARTS. A function handle
            %   that returns a column (one scalar per point) gives that
            %   multiple of the identity, except in the second syntax;
            %   other values are passed on as they are. Any other ALPHA raises
            %   the error fem4m:FormEval:invalidCoefficient.
            if ~isa(type,'FormType')
                error('fem4m:FormEval:invalidType', ...
                    'FormEval: the coefficient type must be a FormType, e.g. FormType.vector2d, not a number');
            end
            if nargin < 6
                target = type;
            end
            if type == FormType.tensor
                type = FormType.scalar;
                target = FormType.scalar;
            end
            n = type.nrComponents;
            if target == type
                shape = [1 n n];
                unit = reshape(eye(n),shape);
            elseif target == FormType.scalar && type.isVector()
                shape = [1 n];
                unit = ones(shape);
            else
                error('fem4m:FormEval:invalidType', ...
                    'FormEval: no coefficient maps %s to %s',char(type),char(target));
            end
            if isempty(alpha)
                val = unit;
            elseif isa(alpha,'numeric')
                if isscalar(alpha)
                    val = alpha*unit;
                elseif numel(alpha) == numel(unit)
                    val = reshape(alpha,shape);
                elseif numel(alpha)>=size(pts,1) && numel(size(alpha))==2 && size(alpha,2)==1
                    val = alpha(parts).*unit;
                else
                    error('fem4m:FormEval:invalidCoefficient','Error evaluation coefficient function');
                end
            elseif isa(alpha,'function_handle')
                if nargin(alpha) == 3
                    val = alpha(pts,pts_ref,parts);
                else
                    val = alpha(pts);
                end
                if size(val,1)==numel(val) && numel(shape)==3
                    val = val.*unit;
                end
            else
                error('fem4m:FormEval:invalidCoefficient','Error evaluation coefficient function');
            end
        end

        function val = evalFunction(pts,pts_ref,parts,f,normals)
            %evalFunction - Evaluate a source function on points
            %   VAL = evalFunction(PTS,PTS_REF,PARTS,F) evaluates the
            %   function handle or constant F at the points PTS. Empty F gives
            %   zeros(numel(PARTS),1), a numeric F the array
            %   F.*ones(numel(PARTS),1), i.e. the row F for every entry of
            %   PARTS. F is called as F(PTS) for 1 argument and as
            %   F(PTS,PTS_REF,PARTS) for 3 arguments; F with 2 or 4 arguments
            %   needs NORMALS (second syntax). Any other F (e.g. one with
            %   varargin) is called as F(PTS), and a failure of this call
            %   raises the error fem4m:FormEval:evaluationFailed.
            %
            %   VAL = evalFunction(PTS,PTS_REF,PARTS,F,NORMALS) also passes the
            %   side normals NORMALS (one row per entry of PARTS): F is called
            %   as F(PTS,NORMALS) for 2 and F(PTS,PTS_REF,PARTS,NORMALS) for 4
            %   arguments. NORMALS may be [] if F does not use it.
            if isempty(f)
                val = zeros(numel(parts),1);
            elseif isa(f,'numeric')
                val = f.*ones(numel(parts),1);
            elseif nargin(f) == 1
                val = f(pts);
            elseif nargin(f) == 2
                val = f(pts,normals);
            elseif nargin(f) == 3
                val = f(pts,pts_ref,parts);
            elseif nargin(f) == 4
                val = f(pts,pts_ref,parts,normals);
            else %fallback
                try
                    val = f(pts);
                catch
                    error('fem4m:FormEval:evaluationFailed','Error evaluating %s specify nargin>1 by stating @(x1,x2,...)func(x1,x2,...)',func2str(f));
                end
            end
        end

    end
end