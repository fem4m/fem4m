classdef Form < handle
    %FORM - Base class for variational forms parsed from strings
    %   OBJ = FORM(STR) parses the form string STR (char or string),
    %   e.g. "grad(u)*grad(v) dx", into the struct array forms of form
    %   terms with operators, domains, and scaling powers. The trial
    %   variable is u and the test variable v.
    %
    %   OBJ = FORM(STR,varStr) uses the variable names of the cell
    %   varStr = {TRIAL,TEST} instead of u and v; a char varStr names
    %   both variables, as does a cell with one name, and an empty
    %   varStr gives {'u','v'}; a varStr of another class raises the
    %   error fem4m:Form:invalidVarStr.
    %
    %   Form is abstract: BilinearForm, LinearForm and EstimatorForm
    %   (and through it ErrorForm) call this constructor with their STR
    %   and varStr.
    %
    %   Every term of STR ends with its integration domain: dx
    %   (elements), ds (all sides, i.e. edges in 2d and faces in 3d),
    %   dsi (interior), dsd (Dirichlet), dsn (Neumann), dss (symmetry)
    %   or dsb (whole boundary, BilinearForm and ErrorForm only). A term
    %   may start with h^K* and p^K*, which multiply it by the integer
    %   power K of the local mesh size and of the polynomial degree.
    %
    %   FORM functions:
    %       assemble       - Assemble the form (abstract)
    %       checkDomains   - Check the integration domains of all terms
    %       checkOperators - Check that the operators of a form term exist
    %       maskKeywords   - Hide operator and function names in a string
    %       unmaskKeywords - Restore the names hidden by maskKeywords
    %       getNormal4ed   - Edge normals of both neighboring elements
    %       getTangent4ed  - Edge tangents of both neighboring elements
    %       sideOf         - Kind and data of the sides of a mesh
    %       directions     - Normals or tangents of both neighboring elements
    %       checkSides     - Check that an operator is defined on the sides
    %       getT1          - Keep the first entry along the last dimension
    %
    %   FORM properties:
    %       forms         - struct array of parsed form terms
    %       nrForms       - number of parsed form terms
    %       varStr        - trial and test variable names
    %
    %   See also BilinearForm, LinearForm, EstimatorForm, ErrorForm

    % Copyright (c) 2026 by the fem4m authors

    properties
        forms                  % struct array of parsed form terms
        nrForms                % number of parsed form terms
        varStr                 % trial and test variable names
    end

    methods (Abstract)
        %ASSEMBLE - Assemble the form
        %   ASSEMBLE(OBJ,...) is implemented by BilinearForm, LinearForm,
        %   EstimatorForm and ErrorForm with their own arguments.
        assemble(obj);
    end

    methods

        function obj = Form(str,varStr)
            %FORM - Parse a form string into form terms
            if nargin < 2 || isempty(varStr)
                obj.varStr = {'u','v'};
            elseif isa(varStr,'char')
                obj.varStr{1} = varStr;
            elseif isa(varStr,'cell')
                obj.varStr = varStr;
            else
                error('fem4m:Form:invalidVarStr','Form: varStr not recognized')
            end
            if numel(obj.varStr)<2
                obj.varStr{2} = obj.varStr{1};
            end
            str = char(str);
            formStr = str;
            str(regexp(str,'[^a-zA-Z0-9\^\+\-\.\(\)\*\:\{\[\}\]]')) = '';
            if isempty(str) || str(1)~='+' && str(1)~='-'
                str = ['+',str];
            end
            [str,keywords] = Form.maskKeywords(str,obj.varStr);
            if strcmpi(obj.varStr{1},obj.varStr{2})
                str = strrep(str,obj.varStr{1},'$');
            elseif ~isempty(strfind(obj.varStr{2},obj.varStr{1}))
                str = strrep(str,obj.varStr{2},'#');
                str = strrep(str,obj.varStr{1},'$');
            else
                str = strrep(str,obj.varStr{1},'$');
                str = strrep(str,obj.varStr{2},'#');
            end
            str = Form.unmaskKeywords(str,keywords);

            [obj.forms,termStart,termEnd] = regexpi(str,'([\+\-])(h\^\-?[0-9\.]+\*?\.?)?(p\^\-?[0-9\.]+\*\.?)?(.*?)(d[xsf][indsb]?)','tokens','start','end');
            Form.checkDomains(formStr,str,obj.forms,termStart,termEnd);
            obj.nrForms = numel(obj.forms);
            for curForm = 1 : obj.nrForms
                % power of h/p
                for k = 2 : 3
                    str = obj.forms{curForm}{k};
                    str(regexpi(str,'[^0-9\-]')) = '';
                    if isempty(str);str='0';end
                    obj.forms{curForm}{k} = str;
                end
                % extract coefficients (of u -> $)
                str = obj.forms{curForm}{4};
                coefficients = regexpi(str,['([a-z]*[0-9]?\(x\))?[\*\.\:]?(?:[a-z,0-9]*[\(\{\[])*','\$','[\)\]\}]*'],'tokens');
                obj.forms{curForm}{6} = [];
                for k = 1 : numel(coefficients)
                    coeff_str = coefficients{k}{1};
                    if ~isempty(coeff_str)
                        str = strrep(str,coeff_str,'');
                        obj.forms{curForm}{6}{k} = strrep(coeff_str,'(x)','');
                    else
                        obj.forms{curForm}{6}{k} = [];
                    end
                end
                % extract function
                func = regexpi(str,'([a-z]*[0-9]?\(x\))','tokens');
                obj.forms{curForm}{7} = [];
                for k = 1 : numel(func)
                    func_str = regexpi(func{k}{1},'([a-z]*[0-9]?\(x\))','tokens');
                    str = strrep(str,func_str{1}{1},'func');
                    obj.forms{curForm}{7}{k} = strrep(func_str{1}{1},'(x)','');
                end
                str = lower(str);
                str = strrep(str,'^','_pow');
                str = strrep(str,'+','_plus_');
                str = strrep(str,'-','_minus_');
                str = strrep(str,'{','average');
                str = strrep(str,'[','jump');
                str(regexp(str,'[^a-zA-Z0-9_$#]')) = '';
                str = strrep(str,'__','_');
                if str(end)=='_';str(end)='';end
                if str(1)=='_';str(1)='';end
                if ~isempty(strfind(obj.forms{curForm}{5},'x'))
                    str = [str,'_T']; %#ok<AGROW>
                else
                    str = [str,'_S']; %#ok<AGROW>
                end
                str = strrep(str,'func_T','_weight_T');
                str = strrep(str,'func_S','_weight_S');
                str = strrep(str,'tr','spur');

                obj.forms{curForm}{8} = [];
                suffix = str(end);
                op1str = regexpi(str,'(?:func)?(?:(_plus_)|(_minus_))?([a-z0-9]*\$[xn]?[n,t]*)','tokens');
                if ~isempty(op1str)
                    % the jump/average of the first operator applies to
                    % all operators, e.g. [a-b] or [a]-[b]
                    if contains(op1str{1}{1},'jump')
                        prefix = 'jump';
                    elseif contains(op1str{1}{1},'average')
                        prefix = 'average';
                    else
                        prefix = '';
                    end
                    for k = 1 : numel(op1str)
                        if ~isempty(prefix)
                            op1str{k}{1} = strrep(op1str{k}{1},prefix,'');
                        end
                        op1str{k}{1} = strrep(op1str{k}{1},'spur','tr');
                        obj.forms{curForm}{8}{k} = [prefix,strrep(op1str{k}{1},'$','w'),'_',suffix];
                    end
                end
                obj.forms{curForm}{9} = [];
                op2str = regexpi(str,'(?:func)?(?:[a-z0-9]*\$[xn]?[n,t]*)?([a-z0-9]*#[xn]?[n,t]*)','tokens');
                if ~isempty(op2str)
                    if contains(op2str{1}{1},'jump')
                        prefix = 'jump';
                    elseif contains(op2str{1}{1},'average')
                        prefix = 'average';
                    else
                        prefix = '';
                    end
                    for k = 1 : numel(op2str)
                        if ~isempty(prefix)
                            op2str{k}{1} = strrep(op2str{k}{1},prefix,'');
                        end
                        op2str{k}{1} = strrep(op2str{k}{1},'spur','tr');
                        obj.forms{curForm}{9}{k} = [prefix,strrep(op2str{k}{1},'#','w'),'_',suffix];
                    end
                end
                obj.forms{curForm}{10} = [];
                signstr = regexpi(str,'((?:_plus_)|(?:_minus_))*','tokens');
                if numel(str)<4 || ~strcmpi(str(1:4),'func')
                    if ~isempty(signstr)
                        signstr{2:end+1} = signstr{:};
                    end
                    signstr{1} = {'_plus_'};
                end
                for k = 1 : numel(signstr)
                    obj.forms{curForm}{10}{k} = strrep(strrep(signstr{k}{1},'_plus_','+'),'_minus_','-');
                end
                powstr = regexpi(str,'pow([1-9]?)','tokens');
                if isempty(powstr); powstr = {{'1'}};end
                obj.forms{curForm}{11} = str2double(powstr{1}{1});
                obj.forms{curForm}{12} = contains(str,'weight');
                str = strrep(str,'$','u');
                str = strrep(str,'#','v');
                str = strrep(str,'spur','tr');
                obj.forms{curForm}{13} = str;
                obj.forms{curForm} = cell2struct(obj.forms{curForm},...
                    {'sign','meshsize_power','degree_power','integrand','domain','coefficients','function','operator1','operator2','operatorsigns','pow','weight','kernel'},...
                    2);
            end
            obj.forms = cell2mat(obj.forms);
        end

        function checkOperators(obj,curForm,ops,known,mesh)
            %checkOperators - Check that the operators of a form term exist
            %   checkOperators(OBJ,CURFORM,OPS,KNOWN,MESH) raises the error
            %   fem4m:<Class>:unknownOperator if one of the operator names
            %   OPS of the form term CURFORM (an entry of OBJ.forms) is not
            %   in the list KNOWN of FormBasis or FormValue methods and OBJ
            %   has no kernel method for the term, with a hint for rot on
            %   the 3d mesh MESH.
            if isempty(ops)
                return
            end
            unknown = ops(~ismember(ops,known));
            if isempty(unknown) || any(strcmp(curForm.kernel,methods(obj)))
                return
            end
            hint = '';
            if isa(mesh,'Mesh3d') && any(contains(unknown,'rot'))
                hint = '; rot is the 2d rotation, use curl in 3d';
            end
            term = strrep(strrep(curForm.integrand,'$',obj.varStr{1}),'#',obj.varStr{2});
            error(['fem4m:',class(obj),':unknownOperator'], ...
                '%s: unknown operator ''%s'' in the term ''%s %s''%s.', ...
                class(obj),unknown{1},term,curForm.domain,hint);
        end
    end

    %% helper functions
    methods (Static)

        function checkDomains(formStr,str,forms,termStart,termEnd)
            %checkDomains - Check the integration domains of all terms
            %   checkDomains(FORMSTR,STR,FORMS,TERMSTART,TERMEND) raises
            %   the error fem4m:Form:missingDomain if a term of the form
            %   string FORMSTR has no integration domain and
            %   fem4m:Form:unknownDomain if it has an unknown one. STR is
            %   the preprocessed form string, FORMS are its parsed terms
            %   and TERMSTART and TERMEND the positions of their first and
            %   last characters in STR. The terms must cover STR and each
            %   must end with one of the domains
            %       dx   elements
            %       ds   all sides, i.e. edges in 2d and faces in 3d
            %       dsi  interior sides
            %       dsd  Dirichlet boundary
            %       dsn  Neumann boundary
            %       dss  symmetry boundary
            %       dsb  whole boundary (BilinearForm and ErrorForm only)
            domains = {'dx','ds','dsi','dsd','dsn','dss','dsb'};
            nextStart = [termStart(2:end),numel(str)+1];
            if isempty(forms) || termStart(1) ~= 1
                error('fem4m:Form:missingDomain', ...
                    'Form: term without integration domain (e.g. dx, ds) in ''%s''',formStr);
            end
            for k = 1 : numel(forms)
                domain = forms{k}{5};
                rest = regexp(str(termEnd(k)+1:end),'^[a-zA-Z0-9]*','match','once');
                if ~isempty(rest) || ~any(strcmpi(domain,domains))
                    hint = '';
                    if startsWith(domain,'df','IgnoreCase',true)
                        hint = ['; the sides are ds, dsi, dsd, dsn, dss and dsb, ', ...
                            'also for the faces in 3d'];
                    end
                    error('fem4m:Form:unknownDomain', ...
                        'Form: unknown integration domain ''%s'' in ''%s''%s',[domain,rest],formStr,hint);
                end
                if termEnd(k)+1 ~= nextStart(k)
                    error('fem4m:Form:missingDomain', ...
                        'Form: term without integration domain (e.g. dx, ds) in ''%s''',formStr);
                end
            end
        end

        function [str,keywords] = maskKeywords(str,varStr)
            %maskKeywords - Hide operator and function names in a string
            %   [STR,KEYWORDS] = maskKeywords(STR,VARSTR) replaces every
            %   function or coefficient token NAME(x) and every run of
            %   operator keywords in front of a parenthesis or of one of
            %   the variable names VARSTR, e.g. curl in curl(u) or curlu,
            %   or rotrot in rotrot(u), by a placeholder character, so that
            %   substituting the variable names cannot change their
            %   letters. KEYWORDS holds the replaced tokens in order of
            %   appearance for unmaskKeywords.
            operators = {'biharmonic','laplace','curl','grad','skew','rot', ...
                'div','dev','eps','sym','tr','d2','d3'};
            next = strjoin([{'\('},regexptranslate('escape',varStr)],'|');
            pattern = ['[a-zA-Z]+[0-9]?\(x\)|(?:',strjoin(operators,'|'),')+(?=',next,')'];
            [keywords,rest] = regexpi(str,pattern,'match','split');
            str = rest{1};
            for k = 1 : numel(keywords)
                str = [str,char(57343+k),rest{k+1}]; %#ok<AGROW>
            end
        end

        function str = unmaskKeywords(str,keywords)
            %unmaskKeywords - Restore the names hidden by maskKeywords
            %   STR = unmaskKeywords(STR,KEYWORDS) replaces the placeholder
            %   characters of maskKeywords by the tokens KEYWORDS.
            for k = 1 : numel(keywords)
                str = strrep(str,char(57343+k),keywords{k});
            end
        end

        function val = getNormal4ed(pts_ref,parts,mesh)
            %getNormal4ed - Edge normals of both neighboring elements
            %   VAL = getNormal4ed(PTS_REF,PARTS,MESH) returns the unit
            %   normals of the edges PARTS of MESH as a
            %   numel(PARTS)-by-DIM-by-2 array, DIM the number of
            %   coordinates: the first slice is normal4ed, the second its
            %   negative, for the second neighboring element. For a
            %   MeshSurface the slices are the co-normals
            %   normal4ed(PTS_REF) in the two neighboring triangles,
            %   and the second slice is zero on boundary edges; other
            %   meshes ignore PTS_REF.
            if isa(mesh,'MeshSurface')
                val = permute(mesh.normal4ed(pts_ref),[1 3 2]);
            else
                normals = mesh.normal4ed;
                val = normals;
                val(:,:,2) = -normals;
            end
            val = val(parts,:,:);
        end

        function val = getTangent4ed(pts_ref,parts,mesh)
            %getTangent4ed - Edge tangents of both neighboring elements
            %   VAL = getTangent4ed(PTS_REF,PARTS,MESH) returns the unit
            %   tangents of the edges PARTS of MESH as a
            %   numel(PARTS)-by-DIM-by-2 array, DIM the number of
            %   coordinates: the first slice is tangent4ed (for a
            %   MeshSurface tangent4ed(PTS_REF), other meshes ignore
            %   PTS_REF), the second its negative, for the second
            %   neighboring element.
            if isa(mesh,'MeshSurface')
                t = mesh.tangent4ed(pts_ref);
            else
                t = mesh.tangent4ed;
            end
            val = t;
            val(:,:,2) = -t;
            val = val(parts,:,:);
        end

        function varargout = sideOf(obj,varargin)
            %sideOf - Kind and data of the sides of a mesh
            %   SIDE = sideOf(MESH) returns 'Face' for the faces of a 3d
            %   mesh (Mesh3d) and 'Edge' for the edges of all other
            %   meshes, i.e. 2d, polygonal and surface meshes. The side
            %   operators (suffix _S) and the ds domains act on these
            %   sides.
            %
            %   SIDE = sideOf(FE) returns the kind of the sides of
            %   FE.mesh.
            %
            %   [V1,V2,...] = sideOf(MESH,NAME1,NAME2,...) returns the
            %   data NAME1, NAME2, ... of these sides, the edge data or
            %   the face data of Mesh3d, and computes only the requested
            %   data. NAME is the edge or face name with s for ed or f,
            %   except for nrSides and diam4s:
            %       'n4s'      - n4ed, n4f: nodes of the sides
            %       'e4s'      - e4ed, e4f: neighboring elements, 0 for
            %                    the second one of boundary sides
            %       's4e'      - ed4e, f4e: sides of the elements
            %       'nrSides'  - nrEdges, nrFaces: number of sides
            %       's4Db'     - ed4Db, f4Db: Dirichlet boundary sides
            %       's4Nb'     - ed4Nb, f4Nb: Neumann boundary sides
            %       's4Sb'     - ed4Sb, f4Sb: symmetry boundary sides
            %       'normal4s' - normal4ed, normal4f: unit normals (the
            %                    co-normals of MeshSurface)
            %       'diam4s'   - length4ed, for faces the longest edge
            %                    max(length4ed(ed4f),[],2)
            %
            %   [V1,V2,...] = sideOf(FE,NAME1,NAME2,...) returns the same
            %   data of FE.mesh and the side data of the element FE (an
            %   FE, FEEdge or FEFace), [] if FE has none (e.g. DGDk has no
            %   dof4ed):
            %       'dof4s'    - dof4ed, dof4f of FE
            %       'degree4s' - degree4ed, degree4f of FE
            %   An unknown NAME, or an element name for a mesh, raises the
            %   error fem4m:Form:unknownSideData.
            mesh = obj;
            isFE = isa(obj,'FE') || isa(obj,'FEEdge') || isa(obj,'FEFace');
            if isFE
                mesh = obj.mesh;
            end
            isFace = isa(mesh,'Mesh3d');
            if nargin == 1
                if isFace
                    varargout = {'Face'};
                else
                    varargout = {'Edge'};
                end
                return
            end
            % the side names and the edge and face names they stand for
            meshNames = {'n4s','n4ed','n4f'; 'e4s','e4ed','e4f'; ...
                's4e','ed4e','f4e'; 'nrSides','nrEdges','nrFaces'; ...
                's4Db','ed4Db','f4Db'; 's4Nb','ed4Nb','f4Nb'; ...
                's4Sb','ed4Sb','f4Sb'; 'normal4s','normal4ed','normal4f'; ...
                'diam4s','length4ed','length4ed'};
            feNames = {'dof4s','dof4ed','dof4f'; ...
                'degree4s','degree4ed','degree4f'};
            varargout = cell(1,numel(varargin));
            for k = 1 : numel(varargin)
                isMeshName = strcmp(varargin{k},meshNames(:,1));
                isFEName = strcmp(varargin{k},feNames(:,1));
                if any(isMeshName)
                    varargout{k} = mesh.(meshNames{isMeshName,2+isFace});
                elseif any(isFEName) && isFE
                    name = feNames{isFEName,2+isFace};
                    if isprop(obj,name)
                        varargout{k} = obj.(name);
                    end
                else
                    error('fem4m:Form:unknownSideData', ...
                        'Form: no side data ''%s'' of %s', ...
                        varargin{k},class(obj));
                end
                % the diameter of a face is its longest edge
                if isFace && strcmp(varargin{k},'diam4s')
                    varargout{k} = max(varargout{k}(mesh.ed4f),[],2);
                end
            end
        end

        function dirs = directions(letters,pts_ref,parts,mesh)
            %DIRECTIONS - Normals or tangents of both neighboring elements
            %   DIRS = DIRECTIONS(LETTERS,PTS_REF,PARTS,MESH) returns for
            %   every letter of LETTERS the unit normals ('n') or tangents
            %   ('t', edges only) of the sides PARTS of MESH (see sideOf)
            %   as numel(PARTS)-by-DIM-by-2 arrays in the cell DIRS, the
            %   second slice for the second neighboring element (see
            %   getNormal4ed, getTangent4ed; face normals are normal4f and
            %   its negative). A repeated letter reuses the array; 't' on
            %   faces raises the error fem4m:Form:noFaceTangent.
            dirs = cell(1,numel(letters));
            isEdge = strcmp(Form.sideOf(mesh),'Edge');
            for k = 1 : numel(letters)
                prev = find(letters(1:k-1) == letters(k),1);
                if ~isempty(prev)
                    dirs{k} = dirs{prev};
                elseif letters(k) == 't' && isEdge
                    dirs{k} = Form.getTangent4ed(pts_ref,parts,mesh);
                elseif letters(k) == 't'
                    error('fem4m:Form:noFaceTangent', ...
                        'Form: faces have no tangent direction');
                elseif isEdge
                    dirs{k} = Form.getNormal4ed(pts_ref,parts,mesh);
                else
                    normals = mesh.normal4f;
                    normals = normals(parts,:);
                    dirs{k} = cat(3,normals,-normals);
                end
            end
        end

        function checkSides(mesh,where,className,opName)
            %checkSides - Check that an operator is defined on the sides
            %   checkSides(MESH,WHERE,CLASSNAME,OPNAME) raises the error
            %   fem4m:<CLASSNAME>:sideDimension unless the sides of MESH
            %   are the edges of a 2d mesh (WHERE = 'planarEdges') or the
            %   faces of a 3d mesh (WHERE = 'faces'), the only sides on
            %   which the side operator OPNAME is defined.
            if strcmp(where,'faces')
                ok = strcmp(Form.sideOf(mesh),'Face');
                sides = 'the faces of 3d meshes';
            else
                ok = strcmp(Form.sideOf(mesh),'Edge') && size(mesh.c4n,2) == 2;
                sides = 'the edges of 2d meshes';
            end
            if ~ok
                error(['fem4m:',className,':sideDimension'], ...
                    '%s: %s is only defined on %s',className,opName,sides);
            end
        end

        function val = getT1(val)
            %getT1 - Keep the first entry along the last dimension
            %   VAL = getT1(VAL) keeps only the first entry along the last
            %   dimension of size(VAL) (trailing singleton dimensions do
            %   not count) of a 2- to 5-dimensional array VAL, e.g. the
            %   values or directions of the first neighboring element of
            %   the sides from an array that holds both in its last
            %   dimension. Arrays with more dimensions are unchanged.
            indx = numel(size(val));
            if indx==5
                val = val(:,:,:,:,1);
            elseif indx==4
                val = val(:,:,:,1);
            elseif indx==3
                val = val(:,:,1);
            elseif indx==2
                val = val(:,1);
            end
        end

    end
end