classdef fem4mResult
    %fem4mResult - File name and contents of level-wise results
    %   R = fem4mResult() creates the description of a results file with
    %   the eight empty slots of the file name
    %   problem_fe_refine_estimator_solver_param1_param2_maxNrDoF that
    %   plotApp reads, and the default labels of the slots. All slots
    %   have to be set before fileName or save is called.
    %
    %   R = fem4mResult(Name=VALUE) also sets the slots and properties
    %   given by name, e.g. fem4mResult(problem='Lshape',fe='Pk'). A
    %   label struct given by label=LABEL only replaces the labels of the
    %   slots it contains.
    %
    %   Numbers are written with p as decimal point and m as minus sign,
    %   e.g. 0.3 as 0p3 and -1e-5 as m1em05, so that the file names
    %   contain no '.', which MATLAB, plotApp and mat2data would take for
    %   the start of the file extension.
    %
    %   Example: Save the errors of an adaptive loop
    %       r = fem4mResult(problem='Lshape',fe='Pk',refine='adaptive', ...
    %           estimator='residual',solver='direct',param1=2, ...
    %           param2=0,maxNrDoF=1e4,label=struct('param1','degree'));
    %       file = r.save(tempdir,nrDoF=[3 21],L2Error=[0.1 0.02])
    %
    %   fem4mResult functions:
    %       fileName     - Name of the results file without extension
    %       save         - Save the level-wise results to a MAT-file
    %       encodeNumber - Number as part of a file name
    %       decodeNumber - Number of a part of a file name
    %       displayParts - Parts of file names as shown by plotApp
    %
    %   fem4mResult properties:
    %       Slots     - Names of the file name slots in their order
    %       problem   - Problem or domain name (slot 1)
    %       fe        - Finite element space name (slot 2)
    %       refine    - Refinement strategy (slot 3)
    %       estimator - Error estimator name (slot 4)
    %       solver    - Solver name (slot 5)
    %       param1    - First parameter (slot 6)
    %       param2    - Second parameter (slot 7)
    %       maxNrDoF  - Maximal number of degrees of freedom (slot 8)
    %       label     - Struct of the slot labels shown by plotApp
    %
    %   See also fem4mPath, mat2data, fem4mModuleTemplate

    % Copyright (c) 2026 by the fem4m authors

    properties (Constant)
        Slots = {'problem','fe','refine','estimator','solver', ...
            'param1','param2','maxNrDoF'} % Names of the file name slots in their order
    end

    properties
        problem    % Problem or domain name (slot 1)
        fe         % Finite element space name (slot 2)
        refine     % Refinement strategy (slot 3)
        estimator  % Error estimator name (slot 4)
        solver     % Solver name (slot 5)
        param1     % First parameter (slot 6)
        param2     % Second parameter (slot 7)
        maxNrDoF   % Maximal number of degrees of freedom (slot 8)
        label     = struct('problem','problem','fe','fe','refine','refine', ...
            'estimator','estimator','solver','solver','param1','', ...
            'param2','','maxNrDoF','maxNrDoF') % Struct of the slot labels shown by plotApp
    end

    methods

        function r = fem4mResult(opts)
            %fem4mResult - Create the description of a results file
            arguments
                opts.?fem4mResult
            end
            names = fieldnames(opts);
            for j = 1 : numel(names)
                if strcmp(names{j},'label')
                    labels = fieldnames(opts.label);
                    for k = 1 : numel(labels)
                        if ~ismember(labels{k},fem4mResult.Slots)
                            error('fem4m:fem4mResult:unknownSlot', ...
                                'Unknown slot ''%s'' in label.',labels{k});
                        end
                        r.label.(labels{k}) = opts.label.(labels{k});
                    end
                else
                    r.(names{j}) = opts.(names{j});
                end
            end
        end

        function name = fileName(r)
            %fileName - Name of the results file without extension
            %   NAME = fileName(R) joins the values of the eight slots of R
            %   by underscores. Numbers are written by encodeNumber. A slot
            %   must not be empty, a number must be a scalar and a text
            %   must not contain '_' or '.', which separate the parts of
            %   the names in plotApp and mat2data.
            parts = cell(1,numel(fem4mResult.Slots));
            for j = 1 : numel(parts)
                slot = fem4mResult.Slots{j};
                value = r.(slot);
                if isempty(value) || (isstring(value) && strlength(value) == 0)
                    error('fem4m:fem4mResult:emptySlot', ...
                        'The slot %s of the file name is empty.',slot);
                end
                if isnumeric(value) || islogical(value)
                    if ~isscalar(value)
                        error('fem4m:fem4mResult:invalidSlotValue', ...
                            'The number of the slot %s must be a scalar.',slot);
                    end
                    parts{j} = fem4mResult.encodeNumber(value);
                else
                    parts{j} = char(value);
                    if contains(parts{j},{'_','.'})
                        error('fem4m:fem4mResult:invalidSlotValue', ...
                            'The value ''%s'' of the slot %s contains ''_'' or ''.''.', ...
                            parts{j},slot);
                    end
                end
            end
            name = strjoin(parts,'_');
        end

        function file = save(r,folder,varargin)
            %SAVE - Save the level-wise results to a MAT-file
            %   FILE = SAVE(R,FOLDER,LEVELS) saves the struct p to the file
            %   fileName(R).mat in FOLDER, creating FOLDER if needed, and
            %   returns the full file name. LEVELS is a scalar struct of
            %   numeric or logical vectors with one entry per level. It
            %   must contain the field nrDoF, and no field name may
            %   contain '_', which separates the parts of the names in
            %   mat2data. p.level(k).NAME is LEVELS.NAME(k). The fields
            %   p.problem, ..., p.maxNrDoF hold the values of the slots
            %   and p.label is the 1-by-8 cell of the slot labels ('' for
            %   a slot without a field in R.label).
            %
            %   FILE = SAVE(R,FOLDER,Name1=VALUE1,Name2=VALUE2,...) gives
            %   the vectors of the levels as name-value arguments.
            if isscalar(varargin) && isstruct(varargin{1})
                levels = varargin{1};
            elseif mod(numel(varargin),2) == 0 && ...
                    all(cellfun(@(s)ischar(s) || (isstring(s) && isscalar(s)),varargin(1:2:end)))
                names = cellfun(@char,varargin(1:2:end),'UniformOutput',false);
                levels = cell2struct(varargin(2:2:end)',names',1);
            else
                error('fem4m:fem4mResult:invalidLevels', ...
                    'The levels must be a struct or name-value arguments.');
            end
            if ~isscalar(levels) || ~isfield(levels,'nrDoF')
                error('fem4m:fem4mResult:missingNrDoF', ...
                    'The levels must contain the field nrDoF.');
            end
            names = fieldnames(levels);
            nrLevels = numel(levels.nrDoF);
            for j = 1 : numel(names)
                if contains(names{j},'_')
                    error('fem4m:fem4mResult:invalidLevelName', ...
                        'The level field name %s must not contain ''_''.',names{j});
                end
                value = levels.(names{j});
                if ~(isnumeric(value) || islogical(value)) || ~isvector(value) ...
                        || numel(value) ~= nrLevels
                    error('fem4m:fem4mResult:levelSizeMismatch', ...
                        'The level field %s must be a numeric vector with %d entries.', ...
                        names{j},nrLevels);
                end
            end
            labels = fieldnames(r.label);
            unknown = setdiff(labels,fem4mResult.Slots);
            if ~isempty(unknown)
                error('fem4m:fem4mResult:unknownSlot', ...
                    'Unknown slot ''%s'' in label.',unknown{1});
            end

            p = struct();
            for j = 1 : numel(fem4mResult.Slots)
                p.(fem4mResult.Slots{j}) = r.(fem4mResult.Slots{j});
            end
            p.label = cell(1,numel(fem4mResult.Slots));
            for j = 1 : numel(fem4mResult.Slots)
                if isfield(r.label,fem4mResult.Slots{j})
                    p.label{j} = char(r.label.(fem4mResult.Slots{j}));
                else
                    p.label{j} = '';
                end
            end
            for k = 1 : nrLevels
                for j = 1 : numel(names)
                    p.level(k).(names{j}) = levels.(names{j})(k);
                end
            end

            if ~isfolder(folder)
                mkdir(folder);
            end
            file = fullfile(folder,[r.fileName(),'.mat']);
            save(file,'p');
        end

    end

    methods (Static)

        function part = encodeNumber(x)
            %encodeNumber - Number as part of a file name
            %   PART = fem4mResult.encodeNumber(X) writes the number X with
            %   15 significant digits and replaces the decimal point by p
            %   and the minus sign by m, e.g. 0.3 by 0p3, -1 by m1 and 1e-5
            %   by 1em05. The plus sign of an exponent is dropped, e.g.
            %   1e+20 is written as 1e20. Integers with at most 15 digits
            %   keep all their digits.
            part = sprintf('%.15g',double(x));
            part = strrep(strrep(strrep(part,'+',''),'.','p'),'-','m');
        end

        function x = decodeNumber(part)
            %decodeNumber - Number of a part of a file name
            %   X = fem4mResult.decodeNumber(PART) returns the number that
            %   encodeNumber writes as PART, and NaN if PART is no such
            %   number, e.g. for the texts p1 or 00.
            part = char(part);
            x = str2double(strrep(strrep(part,'p','.'),'m','-'));
            if isnan(x) || ~isreal(x) || ~strcmp(fem4mResult.encodeNumber(x),part)
                x = NaN;
            end
        end

        function values = displayParts(parts)
            %displayParts - Parts of file names as shown by plotApp
            %   VALUES = fem4mResult.displayParts(PARTS) returns the texts
            %   that plotApp shows for the N-by-8 cell PARTS of the parts
            %   of N file names. The numbers of param1 and param2 are
            %   shown with decimal point and minus sign, e.g. 0p3 as 0.3,
            %   and the trailing zeros of maxNrDoF as exponent, e.g.
            %   150000 as 15e4. Different numbers give different texts.
            if size(parts,2) ~= numel(fem4mResult.Slots)
                error('fem4m:fem4mResult:invalidParts', ...
                    'The parts must have %d columns, one per slot.',numel(fem4mResult.Slots));
            end
            values = parts;
            for k = 1 : size(parts,1)
                for j = 6 : 7
                    x = fem4mResult.decodeNumber(parts{k,j});
                    if ~isnan(x)
                        values{k,j} = sprintf('%.15g',x);
                    end
                end
                part = char(parts{k,8});
                digits = regexprep(part,'0+$','');
                if ~isempty(regexp(part,'^\d+$','once')) && ~isempty(digits) ...
                        && numel(digits) < numel(part)
                    values{k,8} = sprintf('%se%d',digits,numel(part)-numel(digits));
                end
            end
        end

    end

end