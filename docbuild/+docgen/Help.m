classdef Help
    %Help - Read the help text of MATLAB source files
    %   Help provides static methods that extract the help block of a
    %   class or function file and split it into the parts shown on a
    %   reference page: description, member lists, examples, See also
    %   names, notes and the call syntax.
    %
    %   Help functions:
    %       block         - Kind, first code line and help lines of a file
    %       fileSummary   - Kind and summary line of a source file
    %       commentBlock  - Comment lines starting at a line
    %       splitHelp     - Description, member lists and examples
    %       seeAlso       - See also names and the notes after them
    %       findFunction  - Signature and help of a function in a file
    %       classHeader   - Attributes, name and superclasses of a classdef
    %       callLines     - Calls of a name found in help lines
    %       callRe        - Regular expression of a call of a name
    %       normCall      - Call in reference-page style
    %       lowerCaps     - Lower-case an all-caps argument name
    %       argTokens     - Argument names used in the calls of a name
    %       uncomment     - Line without its leading '%' characters
    %       functionLines - Lines that may start a function definition
    %       stripComment  - Code of a line without its comment
    %
    %   See also docgen.RefGenerator, docgen.Text

    % Copyright (c) 2026 by the fem4m authors

    properties (Constant)
        % heading of a member list, e.g. '%   Pk functions:' (any case)
        HeadRe = '^\s*%\s*(?:[\w.]+\s+)?(functions|properties|methods|events)\s*:'
        % entry of a member list, e.g. '%       dof4e - Local dofs'
        EntryRe = '^\s*%\s+([A-Za-z_]\w*)\s*-\s+(.*)'
    end

    methods (Static)

        function [kind, sig, hb, lines] = block(file)
            %block - Kind, first code line and help lines of a file
            %   [KIND,SIG,HB,LINES] = docgen.Help.block(FILE) returns the
            %   KIND ('class', 'function', 'script' or 'comment-only'),
            %   the first code line SIG, the help lines HB (with '%', the
            %   comment block right after the first code line; all comment
            %   lines of a script) and all LINES of FILE.
            lines = docgen.Text.readLines(file);
            [kind, first] = docgen.Help.codeStart(lines);
            switch kind
                case 'comment-only'
                    sig = '';
                    hb = cell(1, 0);
                case 'script'
                    sig = '(script)';
                    hb = lines(cellfun(@(l) startsWith(strtrim(l), ...
                        '%'), lines));
                otherwise
                    sig = strtrim(lines{first});
                    hb = docgen.Help.commentBlock(lines, first + 1);
            end
        end

        function [kind, h1] = fileSummary(lines)
            %fileSummary - Kind and summary line of a source file
            %   [KIND,H1] = docgen.Help.fileSummary(LINES) returns KIND
            %   'class', 'function', 'script' or 'comment-only' and the
            %   first comment line H1 after the classdef line (within 5
            %   lines) or function line (within 9 lines), without '%'.
            [kind, first] = docgen.Help.codeStart(lines);
            h1 = '';
            span = 0;
            if strcmp(kind, 'class')
                span = 5;
            elseif strcmp(kind, 'function')
                span = 9;
            end
            for j = first+1:min(first+span, numel(lines))
                raw = strtrim(lines{j});
                if startsWith(raw, '%')
                    h1 = strtrim(docgen.Help.uncomment(raw));
                    return
                end
            end
        end

        function hb = commentBlock(lines, j)
            %commentBlock - Comment lines starting at line J
            %   HB = docgen.Help.commentBlock(LINES,J) collects the comment
            %   lines from line J on (stripped), skipping empty lines
            %   before the first comment; it stops at the first empty line
            %   after a comment or at a code line.
            hb = cell(1, 0);
            while j <= numel(lines)
                raw = strtrim(lines{j});
                if startsWith(raw, '%')
                    hb{end+1} = raw; %#ok<AGROW>
                elseif isempty(raw)
                    if ~isempty(hb)
                        break
                    end
                else
                    break
                end
                j = j + 1;
            end
        end

        function [desc, lists, examples] = splitHelp(hb)
            %splitHelp - Description, member lists and examples
            %   [DESC,LISTS,EXAMPLES] = docgen.Help.splitHelp(HB) splits
            %   the help lines HB up to 'See also': DESC holds the lines
            %   (without '%') outside of member lists and examples, LISTS
            %   (a dictionary, read with LISTS{NAME}) maps member names to
            %   their text in the member lists
            %   ('<Name> functions:', 'properties:', ...) and EXAMPLES
            %   holds the lines after an 'Examples:' heading.
            desc = cell(1, 0);
            examples = cell(1, 0);
            lists = configureDictionary("string", "cell");
            inList = false;
            inExamples = false;
            cur = '';
            for k = 1:numel(hb)
                ln = hb{k};
                if docgen.Help.isSeeAlso(ln)
                    break
                end
                if ~isempty(regexpi(ln, docgen.Help.HeadRe, 'once'))
                    inList = true;
                    inExamples = false;
                    cur = '';
                    continue
                end
                if ~isempty(regexpi(ln, '^\s*%\s*Examples?\s*:', 'once'))
                    inExamples = true;
                    inList = false;
                    cur = '';
                    continue
                end
                if inList
                    tok = regexp(ln, docgen.Help.EntryRe, 'tokens', 'once');
                    if ~isempty(tok)
                        if ~isKey(lists, tok{1})
                            lists{tok{1}} = strtrim(tok{2});
                        end
                        cur = tok{1};
                    elseif ~isempty(cur) && ~isempty(strtrim(ln))
                        cont = strtrim(docgen.Help.uncomment(ln));
                        if ~isempty(cont)
                            lists{cur} = [lists{cur} ' ' cont];
                        end
                    end
                    continue
                end
                if inExamples
                    t = deblank(docgen.Help.uncomment(ln));
                    if ~isempty(t)
                        examples{end+1} = t; %#ok<AGROW>
                    end
                    continue
                end
                if startsWith(strtrim(ln), '%')
                    desc{end+1} = deblank(docgen.Help.uncomment(ln)); %#ok<AGROW>
                else
                    desc{end+1} = ''; %#ok<AGROW>
                end
            end
        end

        function [names, notes] = seeAlso(hb)
            %seeAlso - See also names and the notes after them
            %   [NAMES,NOTES] = docgen.Help.seeAlso(HB) returns the names
            %   listed after 'See also' (the list ends at the first empty
            %   comment line) and the paragraphs that follow it up to the
            %   copyright line, e.g. 'Note: ...'.
            names = cell(1, 0);
            notes = cell(1, 0);
            state = "before";
            for k = 1:numel(hb)
                ln = hb{k};
                t = strtrim(docgen.Help.uncomment(strtrim(ln)));
                if state == "before"
                    tok = regexpi(strtrim(ln), '^%\s*See also\s+(.*)', ...
                        'tokens', 'once');
                    if isempty(tok)
                        continue
                    end
                    state = "list";
                    t = strtrim(tok{1});
                end
                if state == "list"
                    if isempty(t)
                        state = "after";
                        continue
                    end
                    for n = regexp(t, '[,;]\s*', 'split')
                        name = strtrim(regexprep(strtrim(n{1}), '\.+$', ''));
                        if ~isempty(regexp(name, '^[A-Za-z_][\w.]*$', 'once'))
                            names{end+1} = name; %#ok<AGROW>
                        end
                    end
                    continue
                end
                if startsWith(lower(t), 'copyright')
                    break
                end
                if isempty(t)
                    if ~isempty(notes) && ~isempty(notes{end})
                        notes{end+1} = ''; %#ok<AGROW>
                    end
                    continue
                end
                if ~isempty(notes) && ~isempty(notes{end})
                    notes{end} = [notes{end} ' ' t];
                else
                    if ~isempty(notes)
                        notes(end) = [];
                    end
                    notes{end+1} = t; %#ok<AGROW>
                end
            end
            notes = notes(~cellfun(@isempty, notes));
        end

        function info = findFunction(lines, fname)
            %findFunction - Signature and help of a function in a file
            %   INFO = docgen.Help.findFunction(LINES,FNAME) finds the
            %   first line 'function ... = FNAME(...)' in LINES and returns
            %   a struct with fields outputs, inputs and help (the comment
            %   block after it), or [] if there is none.
            info = [];
            id = '[A-Za-z_][\w.]*';
            for i = docgen.Help.functionLines(lines)
                code = docgen.Help.stripComment(lines{i});
                % no captures inside other groups: MATLAB drops them
                tok = regexp(code, ['^\s*function\s+(\[?[\w,\s~]+\]?\s*=\s*)?' ...
                    '(' id ')\s*(\([^)]*\))?\s*(\([^)]*\))?\s*$'], 'tokens', 'once');
                if ~isempty(tok) && strcmp(tok{2}, fname)
                    if ~isempty(tok{3})
                        inputs = tok{3}(2:end-1);
                    elseif ~isempty(tok{4})
                        inputs = tok{4}(2:end-1);
                    else
                        inputs = '';
                    end
                    outputs = regexprep(tok{1}, '^\[|\]?\s*=\s*$', '');
                    info = struct('outputs', strtrim(outputs), ...
                        'inputs', strtrim(inputs), ...
                        'help', {docgen.Help.commentBlock(lines, i + 1)});
                    return
                end
            end
        end

        function [attrs, name, supers] = classHeader(lines)
            %classHeader - Attributes, name and superclasses of a classdef
            %   [ATTRS,NAME,SUPERS] = docgen.Help.classHeader(LINES) parses
            %   the classdef line: ATTRS is the attribute text, NAME the
            %   class name and SUPERS a cell array of superclass names.
            attrs = '';
            name = '';
            supers = cell(1, 0);
            for i = 1:numel(lines)
                code = strtrim(docgen.Help.stripComment(lines{i}));
                if startsWith(code, 'classdef')
                    tok = regexp(code, ['^classdef\s*(\([^)]*\)\s*)?' ...
                        '([A-Za-z_]\w*)(\s*<\s*.*)?\s*$'], 'tokens', 'once');
                    if ~isempty(tok)
                        attrs = regexprep(tok{1}, '^\(|\)\s*$', '');
                        name = tok{2};
                        superText = regexprep(tok{3}, '^\s*<\s*', '');
                        if ~isempty(superText)
                            supers = cellfun(@strtrim, ...
                                strsplit(superText, '&'), 'UniformOutput', false);
                            supers = supers(~cellfun(@isempty, supers));
                        end
                    end
                    return
                end
            end
        end

        function calls = callLines(hb, name, display)
            %callLines - Calls of a name found in help lines
            %   CALLS = docgen.Help.callLines(HB,NAME,DISPLAY) returns the
            %   calls like 'OBJ = PK(DEGREE,MESH)' of NAME (matched
            %   ignoring case) in the help lines HB before 'See also', in
            %   reference-page style with NAME written as DISPLAY.
            if nargin < 3 || isempty(display)
                display = name;
            end
            calls = cell(1, 0);
            rx = docgen.Help.callRe(name);
            for k = 1:numel(hb)
                ln = hb{k};
                if docgen.Help.isSeeAlso(ln)
                    break
                end
                if startsWith(strtrim(ln), '%')
                    t = strtrim(docgen.Help.uncomment(strtrim(ln)));
                else
                    t = ln;
                end
                toks = regexpi(t, rx, 'tokens');
                for m = 1:numel(toks)
                    out = regexprep(toks{m}{1}, '\s*=\s*$', '');
                    c = docgen.Help.normCall(out, toks{m}{3}, display);
                    if ~any(strcmp(calls, c))
                        calls{end+1} = c; %#ok<AGROW>
                    end
                end
            end
        end

        function c = normCall(out, args, name)
            %normCall - Call in reference-page style
            %   C = docgen.Help.normCall(OUT,ARGS,NAME) returns
            %   'OUT = NAME(ARGS)' (or 'NAME(ARGS)' for empty OUT) with
            %   all-caps argument names in lower case.
            args = docgen.Help.lowerAllCaps(args);
            if ~isempty(out)
                c = [docgen.Help.lowerAllCaps(out) ' = ' name args];
            else
                c = [name args];
            end
        end

        function s = lowerCaps(tok)
            %lowerCaps - Lower-case an all-caps argument name
            %   'DEGREE' -> 'degree'; mixed-case names and single letters
            %   stay as they are.
            if numel(tok) > 1 && docgen.Text.isUpper(tok)
                s = lower(tok);
            else
                s = tok;
            end
        end

        function toks = argTokens(hb, name)
            %argTokens - Argument names used in the calls of a name
            %   TOKS = docgen.Help.argTokens(HB,NAME) returns the argument
            %   and output names as written in the help text (e.g.
            %   'DEGREE', 'varStr') of all calls of NAME in HB.
            toks = cell(1, 0);
            rx = docgen.Help.callRe(name);
            for k = 1:numel(hb)
                ln = hb{k};
                if docgen.Help.isSeeAlso(ln)
                    break
                end
                found = regexpi(docgen.Help.uncomment(ln), rx, 'tokens');
                for m = 1:numel(found)
                    for part = found{m}([1 3])  % outputs and arguments
                        p = regexprep(part{1}, '"[^"]*"|''[^'']*''', '');
                        toks = [toks regexp(p, '[A-Za-z_]\w*', 'match')]; %#ok<AGROW>
                    end
                end
            end
            toks = unique(toks);
        end

        function s = uncomment(s)
            %uncomment - Line without its leading '%' characters
            s = regexprep(s, '^%+', '');
        end

        function rx = callRe(name)
            %callRe - Regular expression of a call of a name
            %   RX = docgen.Help.callRe(NAME) matches calls like
            %   'obj = Pk(degree,mesh)' or 'Pk(degree,mesh)' of NAME; its
            %   tokens are the output part ('obj = ' or ''), NAME and the
            %   argument list with parentheses. NAME must not follow a word
            %   character or '.'.
            rx = ['((?:\[[^\]]*\]|[A-Za-z_]\w*)\s*=\s*)?(?<![\w.])(' ...
                regexptranslate('escape', name) ')\s*(\([^)]*\))'];
        end

        function idx = functionLines(lines)
            %functionLines - Lines that may start a function definition
            %   IDX = docgen.Help.functionLines(LINES) returns the indices
            %   (row vector) of the lines starting with 'function', so the
            %   callers strip the comments of these lines only.
            idx = find(~cellfun(@isempty, regexp(lines, '^\s*function\>', 'once')));
            idx = reshape(idx, 1, []);
        end

        function code = stripComment(line)
            %stripComment - Code of a line without its comment
            %   CODE = docgen.Help.stripComment(LINE) removes the comment
            %   of LINE; a '%' inside a 'char' or "string" literal is
            %   kept. A quote right after a name, number, closing bracket,
            %   '.' or another quote is a transpose (A', x.'), not the
            %   start of a literal.
            quote = '';  % quote character of the open literal, '' outside
            i = 1;
            n = numel(line);
            while i <= n
                c = line(i);
                if ~isempty(quote)
                    if c == quote
                        if i < n && line(i+1) == quote
                            i = i + 1;  % doubled quote inside the literal
                        else
                            quote = '';
                        end
                    end
                elseif c == '"' || (c == '''' && (i == 1 || ...
                        isempty(regexp(line(i-1), '[\w)\]}.''"]', 'once'))))
                    quote = c;
                elseif c == '%'
                    break
                end
                i = i + 1;
            end
            code = line(1:i-1);
        end

    end

    methods (Static, Access = private)

        function [kind, first] = codeStart(lines)
            % kind of a file from its first line that is neither empty nor
            % a comment (FIRST, 0 if there is none)
            first = find(~cellfun(@(l) isempty(strtrim(l)) || ...
                startsWith(strtrim(l), '%'), lines), 1);
            if isempty(first)
                first = 0;
                kind = 'comment-only';
            elseif ~isempty(regexp(lines{first}, '^\s*classdef\>', 'once'))
                kind = 'class';
            elseif ~isempty(regexp(lines{first}, '^\s*function\>', 'once'))
                kind = 'function';
            else
                kind = 'script';
            end
        end

        function tf = isSeeAlso(ln)
            tf = ~isempty(regexpi(ln, '^\s*%\s*See also\>', 'once'));
        end

        function s = lowerAllCaps(s)
            [words, gaps] = regexp(s, '[A-Za-z_]\w*', ...
                'match', 'split');
            words = cellfun(@docgen.Help.lowerCaps, words, 'UniformOutput', false);
            s = [gaps; [words {''}]];
            s = [s{:}];
        end

    end
end