classdef LinearForm < Form
    %LinearForm - Linear form parsed from a form string
    %   OBJ = LinearForm(STR) parses the linear form STR in the test
    %   function v, e.g. "f(x)*v dx" or "g(x)*[v] ds", into form
    %   terms ready to be assembled.
    %
    %   OBJ = LinearForm(STR,varStr) uses the test variable name TEST
    %   of the cell array varStr = {TRIAL,TEST} instead of v, e.g.
    %   {'sigma','tau'} for "f(x)*(tau*n) dsd"; TRIAL must differ from
    %   TEST. A char varStr names both variables, and such a form
    %   cannot be assembled.
    %
    %   Each term holds at most one function NAME(x). The terms may use
    %   the domains dx, ds, dsi, dsd, dsn and dss of Form (not dsb) and
    %   the operators of FormBasis.
    %
    %   LinearForm functions:
    %       assemble - Assemble the load vector
    %       compiler - Build the basis evaluated against a function
    %
    %   See also Form, BilinearForm, FormBasis, FormEval

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = LinearForm(str,varStr)
            %LinearForm - Create a linear form object
            if nargin < 2
                varStr = [];
            end
            obj@Form(str,varStr);
        end

        function [rhs,localRhs] = assemble(obj,fe,degree,subdomain,varargin)
            %ASSEMBLE - Assemble the load vector
            %   RHS = ASSEMBLE(OBJ,FE) assembles the linear form into
            %   the global load vector RHS, using FE as test space. RHS
            %   has nrDoF(FE) entries (allnrDoF for elements that have
            %   it).
            %
            %   RHS = ASSEMBLE(OBJ,FE,DEGREE) integrates polynomials of
            %   degree DEGREE exactly, by default twice the degree of FE
            %   (plus the mesh degree minus one on curved meshes).
            %   DEGREE must not be empty: unlike in BilinearForm, an
            %   empty DEGREE is not replaced by the default.
            %
            %   RHS = ASSEMBLE(OBJ,FE,DEGREE,SUBDOMAIN) assembles only
            %   on the elements (dx terms) or sides (ds terms) given in
            %   SUBDOMAIN; for side terms SUBDOMAIN replaces the sides of
            %   the domain of the term (dsi, dsd, ...).
            %
            %   RHS = ASSEMBLE(...,NAME,VALUE) substitutes VALUE for the
            %   function NAME(x) of the form string, e.g. 'f',F for
            %   "f(x)*v dx", after any of the arguments above. NAME is a
            %   char vector. VALUE is a constant or a function handle
            %   (see FormEval.evalFunction); a function without a value
            %   is zero.
            %
            %   [RHS,localRhs] = ASSEMBLE(...) also returns the local
            %   element or side vectors in localRhs, one row per element
            %   or side assembled: indexed by element or side number for
            %   dx and ds terms without SUBDOMAIN, else by the position in
            %   SUBDOMAIN or in the list of sides of the domain (dsi,
            %   dsd, ...). All terms must then have the same domain (else
            %   the error fem4m:LinearForm:mixedLocalDomains) and local
            %   vectors of the same size, e.g. side terms must not mix
            %   jumps or averages with one-sided operators.
            if nargin < 3
                degree = 2*fe.degree;
                if isprop(fe.mesh,'degree')
                    degree = degree + fe.mesh.degree-1;
                end
            end
            if nargin < 4
                subdomain = [];
            end
            if isa(subdomain,'char')
                varargin = [{subdomain},varargin];
                subdomain = [];
            end
            if isa(degree,'char')
                varargin = [{degree},{subdomain},varargin];
                degree = 2*fe.degree;
                if isprop(fe.mesh,'degree')
                    degree = degree + fe.mesh.degree-1;
                end
                subdomain = [];
            end

            % sides: the edges of 2d, the faces of 3d meshes
            mesh = fe.mesh;
            [s4e,e4s,n4s,nrSides] = Form.sideOf(mesh,'s4e','e4s','n4s','nrSides');
            [s4Db,s4Nb,s4Sb] = Form.sideOf(mesh,'s4Db','s4Nb','s4Sb');
            h_S = Form.sideOf(mesh,'diam4s');
            [fedof4s,degree4s] = Form.sideOf(fe,'dof4s','degree4s');
            if isempty(degree4s)
                degree4s = fe.degree.*ones(nrSides,1);
            end
            if isprop(fe,'degree4e')
                degree4e = fe.degree4e;
            else
                degree4e = fe.degree.*ones(mesh.nrElems,1);
            end
            if ismethod(mesh,'diam4e')
                h_T = mesh.diam4e;
            else
                h_T = max(h_S(s4e),[],2);
            end

            % edge/face - space dofs of both neighboring elements; the
            % second element of boundary sides is replaced by the first,
            % whose side 2 basis values vanish there
            e4s(e4s==0) = 1;
            if isprop(fe, 'dof4e') && ~isempty(fe.dof4e)
                dof4e = fe.dof4e;
                dof4s = [dof4e(e4s(:,1),:),dof4e(e4s(:,2),:)];
            else
                dof4e = [];
                dof4s = [fedof4s,fedof4s];
            end

            % hp-adaptivity
            if isprop(fe, 'allnrDoF')
                rhs = zeros(fe.allnrDoF,1);
            else
                rhs = zeros(fe.nrDoF,1);
            end
            % local vectors need all terms on the same domain
            if nargout > 1 && numel(unique({obj.forms.domain})) > 1
                error('fem4m:LinearForm:mixedLocalDomains', ...
                    'Local vectors (second output) need all terms on the same domain.');
            end
            localRhs = [];

            % loop over forms
            className = class(obj);
            for k = 1 : obj.nrForms
                curForm = obj.forms(k);
                obj.checkOperators(curForm,curForm.operator2,methods(FormBasis),mesh);
                if all(any(strcmp(curForm.operator2,methods(FormBasis)),1))
                    str = ['@(pts,pts_ref,parts)',className,'.compiler(curForm.operator2,pts,pts_ref,parts,fe'];
                else
                    str = ['@(pts,pts_ref,parts)',className,'.',curForm.kernel,'(pts,pts_ref,parts,fe'];
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
                    if isempty(subdomain)
                        indices = 1:fe.mesh.nrElems;
                    else
                        indices = subdomain;
                    end
                    b = integrate(handle,fe.mesh.n4e(indices,:),indices,fe.mesh,degree);
                    b = b.*(h_T(indices).^(str2double(curForm.meshsize_power)).*degree4e(indices).^(str2double(curForm.degree_power)));
                    dofs = dof4e(indices,:);
                elseif any(strcmp(curForm.domain,{'ds','dsi','dsd','dsn','dss'}))
                    if ~isempty(subdomain)
                        indices = subdomain;
                    elseif strcmp(curForm.domain,'ds')
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
                    b = integrate(handle,n4s(indices,:),indices,mesh,degree);
                    b = b.*(h_S(indices).^(str2double(curForm.meshsize_power)).*degree4s(indices).^(str2double(curForm.degree_power)));
                    dofs = dof4s(indices,:);
                    % operators without jump or average (e.g. f(x)*v ds)
                    % only see the first neighboring element of each side
                    if ~iscell(dofs) && size(b,2) < size(dofs,2)
                        dofs = dofs(:,1:size(b,2));
                    end
                else
                    error('fem4m:LinearForm:unknownDomain','FormError: Not implemented');
                end
                if isa(dofs,'cell')
                    dofs = cell2mat(dofs);
                end
                if strcmp(curForm.sign,'+')
                    rhs = rhs + accumarray(dofs(:),b(:),[size(rhs,1),1]);
                    if nargout > 1
                        if isempty(localRhs)
                            localRhs = b;
                        else
                            localRhs = localRhs + b;
                        end
                    end
                else
                    rhs = rhs - accumarray(dofs(:),b(:),[size(rhs,1),1]);
                    if nargout > 1
                        if isempty(localRhs)
                            localRhs = -b;
                        else
                            localRhs = localRhs - b;
                        end
                    end
                end
            end
        end
    end

    methods(Static)

        % compiler
        function val = compiler(op,pts,pts_ref,parts,fe,f)
            %COMPILER - Build the basis evaluated against a function
            %   VAL = COMPILER(OP,PTS,PTS_REF,PARTS,FE,F) builds the
            %   products of the test basis functions of the FormBasis
            %   operator OP{1} on the element FE with the function F (see
            %   FormEval.evalFunction; for side terms an F with 2 or 4
            %   arguments also gets the side normals of PARTS), in the
            %   elements or sides PARTS at the points PTS (one row per
            %   part) and the reference point PTS_REF (a row of
            %   reference coordinates of the element or side, empty for
            %   the polygons of a MeshPoly), as integrate calls the
            %   integrand, for load terms like "f(x)*v dx".
            %   VAL(:,b) belongs to the test basis function b.
            %
            %   VAL = COMPILER(OP,PTS,PTS_REF,PARTS,FE) uses F = 1, for
            %   terms without a function, like "v dx".
            if nargin < 6
                f = 1;
            end
            basis = FormBasis.(op{1})(pts,pts_ref,parts,fe);
            curNormals4s = Form.sideOf(fe.mesh,'normal4s');
            curNormals4s = curNormals4s(parts,:);
            evalF = FormEval.evalFunction(pts,pts_ref,parts,f,curNormals4s);
            basis = basis(:,:,:);
            evalF = evalF(:,:);
            % all entries of symmetric matrices
            type = FormBasis.typeOf(basis,fe);
            if type.isSymmetric()
                basis = basis(:,:,type.full);
            end
            fType = FormValue.typeOf(evalF,fe);
            if fType.isSymmetric()
                evalF = evalF(:,fType.full);
            end
            % multiply f : v
            val = matMul(basis,evalF);
        end

    end
end