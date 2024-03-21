classdef ErrorForm < EstimatorForm
    %ErrorForm - Global error form parsed from a form string
    %   OBJ = ErrorForm(STR) parses the error form STR in the discrete
    %   function u, e.g. "(exactGrad(x)-grad(u))^2 dx", into form
    %   terms ready to be assembled into a single global error value.
    %
    %   OBJ = ErrorForm(STR,varStr) uses the variable name varStr (a
    %   char vector) instead of u, or the two names of the cell array
    %   varStr = {NAME1,NAME2} instead of {'u','v'} for a form in two
    %   discrete functions (see assemble).
    %
    %   The terms are those of EstimatorForm and may also use the
    %   domain dsb (whole boundary).
    %
    %   ErrorForm functions:
    %       assemble - Assemble the global error value
    %
    %   See also EstimatorForm, Form, BilinearForm, FormValue

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = ErrorForm(str,varStr)
            %ErrorForm - Create a global error form object
            if nargin < 2
                varStr = [];
            end
            obj@EstimatorForm(str,varStr);
        end

        function err = assemble(obj,c,fe,degree,varargin) %#ok<INUSD>
            %ASSEMBLE - Assemble the global error value
            %   ERR = ASSEMBLE(OBJ,C,FE) assembles the global error of
            %   the discrete function with the coefficient vector C on
            %   the finite element FE: ERR = (sum of all terms)^(1/p) for
            %   terms (...)^p, e.g. the L2 norm of the error for
            %   "(exact(x)-u)^2 dx". Side terms are integrated over the
            %   sides of their domain. Terms with different powers get
            %   the power of the first term, with the warning
            %   fem4m:ErrorForm:powerMismatch.
            %
            %   ERR = ASSEMBLE(OBJ,C,FE,DEGREE) integrates polynomials of
            %   degree DEGREE exactly, by default twice the degree of FE.
            %
            %   ERR = ASSEMBLE(OBJ,C,FE,C2,FE2) also evaluates the
            %   operators of the second variable name of varStr (v by
            %   default) for the coefficients C2 on the element FE2, e.g.
            %   "(u-v)^2 dx" for the distance of two discrete functions.
            %   The default degree is then twice the larger degree of FE
            %   and FE2. A form without the second name gives the warning
            %   fem4m:EstimatorForm:defaultVarStr.
            %
            %   ERR = ASSEMBLE(OBJ,C,FE,C2,FE2,DEGREE) integrates
            %   polynomials of degree DEGREE exactly.
            %
            %   ERR = ASSEMBLE(...,NAME,VALUE) substitutes VALUE for the
            %   coefficient or function NAME(x) of the form string, e.g.
            %   'exactGrad',GRADU for the form of the class help, after
            %   any of the arguments above. NAME is a char vector. VALUE
            %   is a constant or a function handle (see FormEval); a
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
                    degree = 2*fe.degree ;
                end
            end
            if nargin < 4 || isempty(degree)
                degree = 2*fe.degree;
            end

            % sides: the edges of 2d, the faces of 3d meshes
            mesh = fe.mesh;
            [s4e,n4s,nrSides] = Form.sideOf(mesh,'s4e','n4s','nrSides');
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
                    warning('fem4m:ErrorForm:powerMismatch','Setting power to be equal for all forms (=%s)',num2str(obj.forms(k).pow));
                end
            end

            % loop
            err = 0;
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
                        str = ['@(pts,pts_ref,parts)',className,'.',curForm.kernel,'(pts,pts_ref,parts,c,fe,c2,fe2'];
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
                    curError = integrate(handle,fe.mesh.n4e,1:fe.mesh.nrElems,fe.mesh,degree);
                    curError = h_T.^(str2double(curForm.meshsize_power)).*degree4e.^(str2double(curForm.degree_power)).*curError;
                elseif any(strcmp(curForm.domain,{'ds','dsi','dsd','dsn','dss','dsb'}))
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
                    elseif strcmp(curForm.domain,'dsb')
                        indices = [s4Db;s4Nb;s4Sb];
                    end
                    curError = integrate(handle,n4s(indices,:),indices,mesh,degree);
                    curError = h_S(indices).^(str2double(curForm.meshsize_power)).*degree4s(indices).^(str2double(curForm.degree_power)).*curError;
                else
                    error('fem4m:ErrorForm:unknownDomain','FormError: Domain ''%s'' not implemented',curForm.domain);
                end
                err = err + sum(curError);
            end
            err = err.^(1/curForm.pow);
        end

    end
end