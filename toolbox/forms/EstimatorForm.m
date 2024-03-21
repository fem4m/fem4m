classdef EstimatorForm < Form
    %EstimatorForm - Error estimator form parsed from a form string
    %   OBJ = EstimatorForm(STR) parses the estimator form STR in the
    %   discrete function u, e.g. "h^2*(f(x)+div(grad(u)))^2 dx", into
    %   form terms ready to be assembled into local error indicators.
    %
    %   OBJ = EstimatorForm(STR,varStr) uses the variable name varStr
    %   (a char vector) instead of u, or the two names of the cell
    %   array varStr = {NAME1,NAME2} instead of {'u','v'} for a form in
    %   two discrete functions (see assemble).
    %
    %   A function W(x) after the power of a term, as in
    %   "(exact(x)-u)^2 W(x) dx", weights the term with W (see
    %   compiler). The terms may use the domains dx, ds, dsi, dsd, dsn
    %   and dss of Form (not dsb) and the operators of FormValue.
    %
    %   EstimatorForm functions:
    %       assemble - Assemble the local error indicators
    %       compiler - Build the residual of an estimator term
    %
    %   See also ErrorForm, Form, BilinearForm, FormValue, FormEval

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = EstimatorForm(str,varStr)
            %EstimatorForm - Create an error estimator object
            if nargin < 2
                varStr = [];
            end
            obj@Form(str,varStr);
        end

        function eta = assemble(obj,c,fe,degree,varargin) %#ok<INUSD>
            %ASSEMBLE - Assemble the local error indicators
            %   ETA = ASSEMBLE(OBJ,C,FE) assembles the local error
            %   indicators of the discrete function with the coefficient
            %   vector C on the finite element FE. ETA is nrElems-by-1
            %   with ETA(T) = (sum of the contributions of all terms to
            %   T)^(1/p) for terms (...)^p, i.e. the indicator itself, not
            %   its p-th power. A side term contributes one half to each
            %   of the two neighboring elements of an interior side and
            %   in full to the element of a boundary side. Terms with
            %   different powers get the power of the first term, with
            %   the warning fem4m:EstimatorForm:powerMismatch.
            %
            %   ETA = ASSEMBLE(OBJ,C,FE,DEGREE) integrates polynomials of
            %   degree DEGREE exactly, by default twice the degree of FE.
            %
            %   ETA = ASSEMBLE(OBJ,C,FE,C2,FE2) also evaluates the
            %   operators of the second variable name of varStr (v by
            %   default) for the coefficients C2 on the element FE2, e.g.
            %   "(grad(u)-v)^2 dx" for the distance of grad(u) to a
            %   vector-valued discrete function v. The default degree is
            %   then twice the larger degree of FE and FE2. A form
            %   without the second name gives the warning
            %   fem4m:EstimatorForm:defaultVarStr.
            %
            %   ETA = ASSEMBLE(OBJ,C,FE,C2,FE2,DEGREE) integrates
            %   polynomials of degree DEGREE exactly.
            %
            %   ETA = ASSEMBLE(...,NAME,VALUE) substitutes VALUE for the
            %   coefficient or function NAME(x) of the form string, e.g.
            %   'f',F for "h^2*(f(x)+div(grad(u)))^2 dx", after any of
            %   the arguments above. NAME is a char vector. VALUE is a
            %   constant or a function handle (see FormEval); a
            %   coefficient without a value is the identity, a function
            %   without a value is zero, and a weight W(x) without a
            %   value leaves the term unweighted.
            fe2 = [];
            c2 = []; %#ok<NASGU>
            if numel(varargin)>=1
                if isobject(varargin{1})
                    fe2 = varargin{1};
                    c2 = degree; %#ok<NASGU>
                    if numel(varargin)>=2 && isa(varargin{2},'numeric') && isscalar(varargin{2})
                        degree=varargin{2};
                        varargin = varargin(3:end);
                    else
                        degree = 2*max(fe.degree,fe2.degree);
                        varargin = varargin(2:end);
                    end
                elseif isa(degree,'char')
                    varargin = [{degree},varargin];
                    degree = 2*fe.degree;
                end
            end
            if nargin < 4 || isempty(degree)
                degree = 2*fe.degree;
            end

            % sides: the edges of 2d, the faces of 3d meshes
            mesh = fe.mesh;
            [s4e,e4s,n4s,nrSides] = Form.sideOf(mesh,'s4e','e4s','n4s','nrSides');
            [s4Db,s4Nb,s4Sb] = Form.sideOf(mesh,'s4Db','s4Nb','s4Sb');
            h_S = Form.sideOf(mesh,'diam4s');
            degree4s = Form.sideOf(fe,'degree4s');
            if isempty(degree4s)
                degree4s = fe.degree.*ones(nrSides,1);
            end
            if isprop(fe, 'degree4e')
                degree4e = fe.degree4e;
            else
                degree4e = fe.degree*ones(fe.mesh.nrElems,1);
            end
            if ismethod(mesh,'diam4e')
                h_T = mesh.diam4e;
            else
                h_T = max(h_S(s4e),[],2);
            end

            % powers equal
            for k = 2 : obj.nrForms
                if obj.forms(k).pow~=obj.forms(k-1).pow
                    obj.forms(k).pow = obj.forms(k-1).pow;
                    warning('fem4m:EstimatorForm:powerMismatch','Setting power to be equal for all forms (=%s)',num2str(obj.forms(k).pow));
                end
            end

            %loop
            eta = zeros(fe.mesh.nrElems,1);
            className = class(obj);
            for k = 1 : obj.nrForms
                curForm = obj.forms(k);
                obj.checkOperators(curForm,[curForm.operator1,curForm.operator2],methods(FormValue),mesh);
                op = true;
                for i = 1 : numel(curForm.operator1)
                    op = op & all(any(strcmp(curForm.operator1{i},methods(FormValue)),1));
                end
                for i = 1 : numel(curForm.operator2)
                    op = op & all(any(strcmp(curForm.operator2{i},methods(FormValue)),1));
                end
                if op
                    if isempty(fe2)
                        str = ['@(pts,pts_ref,parts)',className,'.compiler(curForm.operator1,[],curForm.operatorsigns,curForm.weight,curForm.pow,pts,pts_ref,parts,{c},{fe}'];
                    else
                        str = ['@(pts,pts_ref,parts)',className,'.compiler(curForm.operator1,curForm.operator2,curForm.operatorsigns,curForm.weight,curForm.pow,pts,pts_ref,parts,{c,c2},{fe,fe2}'];
                    end
                else
                    if ~isempty(fe2)
                        str = ['@(pts,pts_ref,parts)ErrorForm.',curForm.kernel,'(pts,pts_ref,parts,c,fe,c2,fe2'];
                    else
                        str = ['@(pts,pts_ref,parts)',className,'.',curForm.kernel,'(pts,pts_ref,parts,c,fe'];
                    end
                end
                for curCoeff = 1 : numel(curForm.coefficients)
                    func = '[]';
                    for l = 1 : numel(varargin)/2
                        if strcmp(varargin{2*(l-1)+1},curForm.coefficients{curCoeff})
                            func = ['varargin{',int2str(2*(l-1)+2),'}'];
                        end
                    end
                    str = [str,',',func]; %#ok<AGROW>
                end
                for curFunc = 1 : numel(curForm.function)
                    func = '[]';
                    for l = 1 : numel(varargin)/2
                        if strcmp(varargin{2*(l-1)+1},curForm.function{curFunc})
                            func = ['varargin{',int2str(2*(l-1)+2),'}'];
                        end
                    end
                    str = [str,',',func]; %#ok<AGROW>
                end
                str = [str,')']; %#ok<AGROW>
                handle = eval(str);
                if strcmp(curForm.domain,'dx')
                    eta4e = integrate(handle,fe.mesh.n4e,1:fe.mesh.nrElems,fe.mesh,degree);
                    eta4e = h_T.^(str2double(curForm.meshsize_power)).*degree4e.^(str2double(curForm.degree_power)).*eta4e;
                    eta = eta + eta4e;
                elseif any(strcmp(curForm.domain,{'ds','dsi','dsd','dsn','dss'}))
                    if strcmp(curForm.domain,'ds')
                        indices = 1:nrSides;
                    elseif strcmp(curForm.domain,'dsi')
                        indices = setdiff(1:nrSides,[s4Db;s4Nb;s4Sb]);
                    elseif strcmp(curForm.domain,'dsd')
                        indices = s4Db;
                    elseif strcmp(curForm.domain,'dsn')
                        indices = s4Nb;
                    elseif strcmp(curForm.domain,'dss')
                        indices = s4Sb;
                    end
                    eta4s = integrate(handle,n4s(indices,:),indices,mesh,degree);
                    eta4s = .5*h_S(indices).^(str2double(curForm.meshsize_power)).*degree4s(indices).^(str2double(curForm.degree_power)).*eta4s;
                    eta4s(e4s(indices,2)==0) = 2*eta4s(e4s(indices,2)==0);
                    eta = eta + accumarray(e4s(indices,1),eta4s,size(eta));
                    if any(e4s(indices,2)~=0)
                        eta = eta + accumarray(e4s(indices(e4s(indices,2)~=0),2),eta4s(e4s(indices,2)~=0),size(eta));
                    end
                else
                    error('fem4m:EstimatorForm:unknownDomain','Form error: Form is not (yet) implemented: %s',curForm.domain);
                end
            end
            % return the element vector with the entries \eta_T (not \eta_T^p)
            eta = eta.^(1/curForm.pow);
        end

    end

    methods(Static)

        function val = compiler(op1,op2,sign,weight,pow,pts,pts_ref,parts,c,fe,varargin)
            %COMPILER - Build the residual of an estimator term
            %   VAL = COMPILER(OP1,OP2,SIGN,WEIGHT,POW,PTS,PTS_REF,PARTS,C,FE,ALPHA1,...,F)
            %   builds the integrand of an estimator or error term in the
            %   elements or sides PARTS at the points PTS (one row per
            %   part) and the reference point PTS_REF (a row of
            %   reference coordinates of the element or side, empty for
            %   the polygons of a MeshPoly), as integrate calls it.
            %   The FormValue operators of the cell array OP1 act on the
            %   coefficients C{1} on the element FE{1}, those of OP2
            %   (empty for a form in one function) on C{2} on FE{2}. Each
            %   value of OP1 is multiplied by its coefficient ALPHAi (one
            %   per operator of OP1, empty for the identity), all values
            %   get the signs SIGN{i} and are summed, and the function F
            %   is added. VAL is the column |val|^POW, summed over the
            %   components. WEIGHT is false.
            %
            %   VAL = COMPILER(...,ALPHA1,...) without F uses F = 0.
            %
            %   VAL = COMPILER(...,W) with WEIGHT true weights the term
            %   with the last argument W (after F, if any): VAL =
            %   val'*W*val, POW must be 2. W is evaluated by
            %   FormEval.evalCoefficient (a constant, a column of one
            %   value per part, or a function handle). A constant 3-by-3
            %   (2d) or 6-by-6 (3d) W acts on the scaled Voigt notation
            %   of matrix values.
            if numel(fe)==2 && isempty(op2)
                warning('fem4m:EstimatorForm:defaultVarStr','Specify variable names different from {''u'',''v''}');
            end
            nOp1 = numel(op1);
            nOp2 = numel(op2);
            op = [op1,op2];
            if weight
                weight = varargin{end};
                varargin(end) = [];
            else
                weight = [];
            end
            if numel(varargin) <= nOp1
                f = [];
                alpha = varargin(1:end);
            else
                f = varargin{end};
                alpha = varargin(1:end-1);
            end
            feIdx = [ones(1,nOp1),2*ones(1,nOp2)];
            if nOp2 > 0
                alpha{end+1:end+nOp2} = [];
            end
            curNormals4s = Form.sideOf(fe{1}.mesh,'normal4s');
            curNormals4s = curNormals4s(parts,:);
            fVal = FormEval.evalFunction(pts,pts_ref,parts,f,curNormals4s);
            for i = 1 : numel(op)
                curOp = op{i};
                if nargin(str2func(['FormValue.',curOp]))==6
                    % operators that take the coefficient as argument
                    curUVal = FormValue.(curOp)(pts,pts_ref,parts,c{feIdx(i)},fe{feIdx(i)},alpha{i});
                else
                    % operators without a coefficient argument
                    curUVal = FormValue.(curOp)(pts,pts_ref,parts,c{feIdx(i)},fe{feIdx(i)});
                    % eval coefficient: a matrix acting on the values, a
                    % scalar for scalars and tensors (e.g. vector Laplace)
                    type = FormType.of(size(curUVal,[2 3]),size(pts,2));
                    curAlpha = FormEval.evalCoefficient(pts,pts_ref,parts,alpha{i},type);
                    % multiply coefficient
                    if numel(curAlpha)==size(curAlpha,1)
                        curUVal = curAlpha.*curUVal;
                    elseif all(size(curAlpha,[1 2 3])==[size(curUVal,[1 2]),1])
                        curUVal = sum(curAlpha.*curUVal,2);
                    else
                        curUVal = matMul(curAlpha,curUVal);
                    end
                end
                type = FormType.of(size(curUVal,[2 3]),size(pts,2));
                if type.isSymmetric() % sym to full matrix
                    curUVal = curUVal(:,type.full);
                end
                % sub/add value
                if i==1
                    uVal = (-1).^(strcmp(sign{i},'-')).*curUVal;
                else
                    uVal = uVal + (-1).^(strcmp(sign{i},'-')).*curUVal;
                end
            end
            % f+u
            fType = FormType.of(size(fVal,2),size(pts,2));
            if fType.isSymmetric() % sym to full matrix
                fVal = fVal(:,fType.full);
            elseif size(fVal,2)==size(pts,2) && size(fVal,3)==size(pts,2) && size(uVal,2)==size(pts,2)^2 % dim x dim matrix
                fVal = fVal(:,:);
            end
            val = fVal + uVal;
            % return |val|^p, or val'*W*val for a weight W
            if ~isempty(weight)
                if pow~=2;error('fem4m:EstimatorForm:weightedPowerNotImplemented','weighted p-norm not implemented');end
                % a weight of the size of a symmetric matrix acts on the
                % (scaled) Voigt notation of the full matrix values
                type = FormType.of(size(val,2),size(pts,2));
                if type == FormType.matrix2d && isequal(size(weight),[3 3])
                    type = FormType.symmatrix2d;
                    val = [val(:,1),val(:,4),sqrt(2)/2*(val(:,2)+val(:,3))]; % a{1,1},a{2,2},sqrt(2)*a{1,2}
                elseif type == FormType.matrix3d && isequal(size(weight),[6 6])
                    type = FormType.symmatrix3d;
                    % a{1,1},a{2,2},a{3,3},sqrt(2)*a{1,2},sqrt(2)*a{1,3},sqrt(2)*a{2,3}
                    val = [val(:,1),val(:,5),val(:,9),sqrt(2)/2*(val(:,2)+val(:,4)),sqrt(2)/2*(val(:,3)+val(:,7)),sqrt(2)/2*(val(:,6)+val(:,8))];
                end
                weight = FormEval.evalCoefficient(pts,pts_ref,parts,weight,type);
                if numel(weight)==size(weight,1)
                    weightVal = weight.*val;
                else
                    weightVal = matMul(weight,val);
                end
                val = sum(weightVal.*val,2);
            else
                val = sum(abs(val).^pow,[2,3]);
            end
        end

    end
end