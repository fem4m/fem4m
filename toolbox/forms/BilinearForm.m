classdef BilinearForm < Form
    %BilinearForm - Bilinear form parsed from a form string
    %   OBJ = BilinearForm(STR) parses the bilinear form STR in the
    %   trial function u and the test function v, e.g.
    %   "grad(u)*grad(v) dx", into form terms ready to be assembled.
    %
    %   OBJ = BilinearForm(STR,varStr) uses the trial and test
    %   variable names of the cell array varStr = {TRIAL,TEST} instead
    %   of {'u','v'}, e.g. {'sigma','tau'}. The two names must differ;
    %   a char varStr names both variables, and such a form cannot be
    %   assembled.
    %
    %   The terms may use all domains of Form (dsb included) and the
    %   operators of FormBasis.
    %
    %   BilinearForm functions:
    %       assemble - Assemble the system matrix
    %       compiler - Build the product of two basis operators
    %
    %   See also Form, LinearForm, FormBasis, FormEval

    % Copyright (c) 2026 by the fem4m authors

    methods

        function obj = BilinearForm(str,varStr)
            %BilinearForm - Create a bilinear form object
            if nargin < 2
                varStr = [];
            end
            obj@Form(str,varStr);
        end

        function [A,localA] = assemble(obj,fe1,fe2,degree,subdomain,varargin)
            %ASSEMBLE - Assemble the system matrix
            %   A = ASSEMBLE(OBJ,FE1) assembles the bilinear form into
            %   the sparse global matrix A, using FE1 for trial and
            %   test space.
            %
            %   A = ASSEMBLE(OBJ,FE1,FE2) assembles the form with trial
            %   space FE1 and test space FE2, e.g. for mixed
            %   formulations. A is nrDoF(FE2)-by-nrDoF(FE1), with
            %   allnrDoF instead of nrDoF for elements that have it
            %   (hp-adaptivity).
            %
            %   A = ASSEMBLE(OBJ,FE1,FE2,DEGREE) integrates polynomials
            %   of degree DEGREE exactly, by default the sum of the
            %   degrees of FE1 and FE2 (plus the mesh degree minus one on
            %   curved meshes).
            %
            %   A = ASSEMBLE(OBJ,FE1,FE2,DEGREE,SUBDOMAIN) assembles only
            %   on the elements (dx terms) or sides (ds terms) given in
            %   SUBDOMAIN; for side terms SUBDOMAIN replaces the sides of
            %   the domain of the term (dsi, dsd, ...). Empty FE2, DEGREE
            %   or SUBDOMAIN take the defaults.
            %
            %   A = ASSEMBLE(OBJ,FE1,DEGREE,SUBDOMAIN) is the same with
            %   FE2 = FE1. With only one numeric argument after FE1, a
            %   scalar is the degree and a vector the subdomain, so
            %   ASSEMBLE(OBJ,FE1,SUBDOMAIN) also works; with two, the
            %   first is the degree, unless it is empty: then it is read
            %   as an empty FE2 and the second is the degree. So a
            %   SUBDOMAIN with a single entry must follow a DEGREE, as in
            %   ASSEMBLE(OBJ,FE1,DEGREE,SUBDOMAIN), or
            %   ASSEMBLE(OBJ,FE1,[],[],SUBDOMAIN) for the default degree.
            %
            %   A = ASSEMBLE(...,NAME,VALUE) substitutes VALUE for the
            %   coefficient NAME(x) in front of the trial operator in the
            %   form string, e.g. 'alpha',ALPHA for
            %   "alpha(x)*grad(u)*grad(v) dx", after any of the arguments
            %   above. NAME is a char vector. VALUE is a constant, a
            %   column of one value per element or side, or a function
            %   handle (see FormEval.evalCoefficient); a coefficient
            %   without a value is the identity (ones for a vector
            %   coefficient b as in "b(x)*grad(u)*v dx").
            %
            %   [A,localA] = ASSEMBLE(...) also returns the local
            %   element or side matrices in localA, indexed by element
            %   or side number: localA(K,a,b) belongs to the trial basis
            %   function a and the test basis function b of the element
            %   or side K, and is zero for the elements or sides not
            %   assembled. All terms must then be element (dx) terms or
            %   all side (ds) terms (else the error
            %   fem4m:BilinearForm:mixedLocalDomains), with local matrices
            %   of the same size, e.g. side terms must not mix jumps or
            %   averages with one-sided operators.

            % parse input
            if nargin < 3 || isempty(fe2)
                fe2 = fe1;
            end
            if nargin < 4
                degree = [];
            end
            if nargin < 5
                subdomain = [];
            end
            if isa(fe2,'char')
                if nargin == 3
                    varargin = [{fe2},varargin];
                elseif nargin < 5
                    varargin = [{fe2},{degree},varargin];
                else
                    varargin = [{fe2},{degree},{subdomain},varargin];
                end
                fe2 = fe1;
                degree = [];
                subdomain = [];
            end
            if isa(fe2,'numeric') && nargin < 4
                if isscalar(fe2)
                    degree = fe2;
                    subdomain = [];
                else
                    degree = [];
                    subdomain = fe2;
                end
                fe2 = fe1;
            end
            if isa(fe2,'numeric') && isa(degree,'char')
                varargin = [{degree},{subdomain},varargin];
                if isscalar(fe2)
                    degree = fe2;
                    subdomain = [];
                else
                    degree = [];
                    subdomain = fe2;
                end
                fe2 = fe1;
            end
            if isa(fe2,'numeric') && isa(degree,'numeric')
                if nargin > 4
                    varargin = [{subdomain},varargin];
                end
                subdomain = degree;
                degree = fe2;
                fe2 = fe1;
            end
            if isa(degree,'char')
                varargin = [{degree},{subdomain},varargin];
                degree = [];
                subdomain = [];
            end
            if isa(subdomain,'char')
                varargin = [{subdomain},varargin];
                subdomain = [];
            end
            if isempty(degree)
                degree = fe1.degree + fe2.degree;
                if isprop(fe1.mesh,'degree')
                    degree = degree + fe1.mesh.degree-1;
                end
            end

            % sides: the edges of 2d, the faces of 3d meshes
            mesh = fe1.mesh;
            [s4e,e4s,n4s,nrSides] = Form.sideOf(mesh,'s4e','e4s','n4s','nrSides');
            [s4Db,s4Nb,s4Sb] = Form.sideOf(mesh,'s4Db','s4Nb','s4Sb');
            h_S = Form.sideOf(mesh,'diam4s');
            [fe1dof4s,degree4s] = Form.sideOf(fe1,'dof4s','degree4s');
            fe2dof4s = Form.sideOf(fe2,'dof4s');
            if isempty(degree4s)
                degree4s = fe1.degree.*ones(nrSides,1);
            end
            if isprop(fe1, 'degree4e')
                degree4e = fe1.degree4e;
            else
                degree4e = fe1.degree*ones(mesh.nrElems,1);
            end
            if ismethod(mesh,'diam4e')
                h_T = mesh.diam4e;
            else
                h_T = max(h_S(s4e),[],2);
            end

            % edge/face - space dofs
            e4s(e4s==0) = 1;
            if isprop(fe1, 'dof4e') && ~isempty(fe1.dof4e)
                dof4e1 = fe1.dof4e;
                dof4s1 = [dof4e1(e4s(:,1),:),dof4e1(e4s(:,2),:)];
            else
                dof4e1 = [];
                dof4s1 = [fe1dof4s,fe1dof4s];
            end
            if isprop(fe2, 'dof4e') && ~isempty(fe2.dof4e)
                dof4e2 = fe2.dof4e;
                dof4s2 = [dof4e2(e4s(:,1),:),dof4e2(e4s(:,2),:)];
            else
                dof4e2 = [];
                dof4s2 = [fe2dof4s,fe2dof4s];
            end

            % hp-adaptivity
            if isprop(fe1, 'allnrDoF') && isprop(fe2, 'allnrDoF')
                A = sparse(fe2.allnrDoF,fe1.allnrDoF);
            elseif isprop(fe1, 'allnrDoF')
                A = sparse(fe2.nrDoF,fe1.allnrDoF);
            elseif isprop(fe2, 'allnrDoF')
                A = sparse(fe2.allnrDoF,fe1.nrDoF);
            else
                A = sparse(fe2.nrDoF,fe1.nrDoF);
            end

            % local matrices need one kind of domain: element or side terms
            if nargout > 1
                isDx = strcmp({obj.forms.domain},'dx');
                if any(isDx) && ~all(isDx)
                    error('fem4m:BilinearForm:mixedLocalDomains', ...
                        'Local matrices (second output) need all terms on elements (dx) or all on sides (ds...).');
                end
            end

            % loop over forms
            localA = [];
            className = class(obj);
            for k = 1 : obj.nrForms
                curForm = obj.forms(k);
                obj.checkOperators(curForm,[curForm.operator1,curForm.operator2],methods(FormBasis),mesh);
                if all(any(strcmp(curForm.operator1,methods(FormBasis)),1)) && all(any(strcmp(curForm.operator2,methods(FormBasis)),1))
                    func = ['@(pts,pts_ref,parts)',className,'.compiler([curForm.operator1,curForm.operator2],pts,pts_ref,parts,{fe1,fe2}'];
                else
                    func = ['@(pts,pts_ref,parts)',className,'.',curForm.kernel,'(pts,pts_ref,parts,fe1,fe2'];
                end
                for curCoeff = 1 : numel(curForm.coefficients)
                    coeff = '[]';
                    for l = 1 : numel(varargin)/2
                        if strcmp(varargin{2*(l-1)+1},curForm.coefficients{curCoeff})
                            coeff = ['varargin{',int2str(2*(l-1)+2),'}'];
                        end
                    end
                    func = [func,',',coeff]; %#ok<AGROW>
                end
                func = [func,')']; %#ok<AGROW>
                handle = eval(func);
                if strcmp(curForm.domain,'dx')
                    if isempty(subdomain)
                        indices = 1:fe1.mesh.nrElems;
                    else
                        indices = subdomain;
                    end
                    S = integrate(handle,mesh.n4e(indices,:),indices,mesh,degree);
                    S = S.*(h_T(indices).^(str2double(curForm.meshsize_power)).*degree4e(indices).^(str2double(curForm.degree_power)));
                    [I,J] = localDoFtoGlobalDoF(dof4e1(indices,:),dof4e2(indices,:));
                elseif any(strcmp(curForm.domain,{'ds','dsi','dsd','dsn','dss','dsb'}))
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
                    elseif strcmp(curForm.domain,'dsb')
                        indices = [s4Db;s4Nb;s4Sb];
                    end
                    S = integrate(handle,n4s(indices,:),indices,mesh,degree);
                    S = S.*(h_S(indices).^(str2double(curForm.meshsize_power)).*degree4s(indices).^(str2double(curForm.degree_power)));
                    dofs1 = dof4s1(indices,:);
                    dofs2 = dof4s2(indices,:);
                    % operators without jump or average (e.g. u*v ds) only
                    % see the first neighboring element of each side
                    if ~iscell(dofs1) && size(S,2) < size(dofs1,2)
                        dofs1 = dofs1(:,1:size(S,2));
                    end
                    if ~iscell(dofs2) && size(S,3) < size(dofs2,2)
                        dofs2 = dofs2(:,1:size(S,3));
                    end
                    [I,J] = localDoFtoGlobalDoF(dofs1,dofs2);
                else
                    error('fem4m:BilinearForm:unknownDomain','FormError: Domain ''%s'' not implemented',curForm.domain);
                end
                if nargout > 1 && isempty(localA)
                    sz = size(S);
                    if strcmp(curForm.domain,'dx')
                        sz(1) = fe1.mesh.nrElems;
                    else
                        sz(1) = nrSides;
                    end
                    localA = zeros(sz);
                end
                if strcmp(curForm.sign,'+')
                    A = A + sparse(I(:),J(:),S(:),size(A,1),size(A,2));
                    if nargout > 1
                        localA(indices,:,:,:,:) = localA(indices,:,:,:,:) + S; %#ok<AGROW>
                    end
                else
                    A = A - sparse(I(:),J(:),S(:),size(A,1),size(A,2));
                    if nargout > 1
                        localA(indices,:,:,:,:) = localA(indices,:,:,:,:) - S; %#ok<AGROW>
                    end
                end
            end
        end

    end

    methods(Static)

        %% compiler
        function val = compiler(op,pts,pts_ref,parts,fe,alpha)
            %COMPILER - Build the product of two basis operators
            %   VAL = COMPILER(OP,PTS,PTS_REF,PARTS,FE,ALPHA) builds the
            %   products of the trial basis functions of the FormBasis
            %   operator OP{1} on the element FE{1}, multiplied by the
            %   coefficient ALPHA (see FormEval.evalCoefficient; an empty
            %   ALPHA is the identity), with the test basis functions of
            %   the operator OP{2} on the element FE{2}, in the elements or
            %   sides PARTS at the points PTS (one row per part) and the
            %   reference point PTS_REF (a row of reference coordinates
            %   of the element or side, empty for the polygons of a
            %   MeshPoly), as integrate calls the integrand.
            %   VAL(:,a,b) belongs to the trial basis function a and the
            %   test basis function b.
            basis1 = FormBasis.(op{1})(pts,pts_ref,parts,fe{1});
            basis2 = FormBasis.(op{2})(pts,pts_ref,parts,fe{2});
            type1 = FormBasis.typeOf(basis1,fe{1});
            type2 = FormBasis.typeOf(basis2,fe{2});
            % eval coefficient: vector*scalar (e.g. b*grad(u)*v) needs a
            % vector, all other types a matrix acting on type1
            if type1.isVector() && type2 == FormType.scalar
                curAlpha = FormEval.evalCoefficient(pts,pts_ref,parts,alpha,type1,FormType.scalar);
            else
                curAlpha = FormEval.evalCoefficient(pts,pts_ref,parts,alpha,type1);
            end
            % multiply coefficient
            if numel(curAlpha)==size(curAlpha,1)
                basis1 = curAlpha.*basis1;
            else
                basis1 = matMul(basis1,curAlpha);
            end
            % multiply u : v with all entries of symmetric matrices
            if type1.isSymmetric()
                basis1 = basis1(:,:,type1.full);
            end
            if type2.isSymmetric()
                basis2 = basis2(:,:,type2.full);
            end
            val = matMul(basis1(:,:,:),permute(basis2(:,:,:),[1 3 2]));
        end

    end
end