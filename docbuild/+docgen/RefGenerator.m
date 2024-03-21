classdef RefGenerator < handle
    %RefGenerator - Generate the HTML reference pages of the toolbox
    %   OBJ = docgen.RefGenerator(CTX) prepares the generation of the
    %   reference pages for the classes and functions of the toolbox
    %   described by the docgen.Context CTX (page tree, examples, name
    %   lookup), to be written to CTX.DocFolder. The toolbox must be on the
    %   path, since class members are read with meta.class.
    %
    %   Pages follow the MathWorks layout. Functions: Syntax, Description,
    %   Examples, Tips, See Also. Classes: Description, Creation (Syntax
    %   and Description), Properties, Methods, Examples, Tips, Method
    %   Details, See Also.
    %
    %   RefGenerator functions:
    %       run        - Generate the pages of some or all folders
    %       render     - HTML of the reference page of one source file
    %       printNotes - Print what may need a look in the help texts
    %
    %   RefGenerator properties:
    %       Context      - Folders and toolbox data (docgen.Context)
    %       EmptyDesc    - Members without any description
    %       SyntaxLog    - Notes on generated Syntax sections
    %
    %   See also genRef, docgen.Context, docgen.Help, docgen.classMembers

    % Copyright (c) 2026 by the fem4m authors

    properties
        EmptyDesc = cell(0, 3)  % members without any description: src, kind, name
        SyntaxLog = cell(0, 2)  % notes on generated Syntax sections: src, message
    end

    properties (SetAccess = private)
        Context     % folders and toolbox data (docgen.Context)
    end

    properties (Access = private)
        SubMap      % class name -> direct subclasses
        Classes     % source file -> cached class information
    end

    properties (Constant)
        MaxRelated = 6  % example pages listed on a reference page
    end

    methods

        function obj = RefGenerator(ctx)
            %RefGenerator - Prepare the generation of the reference pages
            obj.Context = ctx;
            obj.Classes = configureDictionary("string", "cell");
        end

        function pages = run(obj, folders)
            %run - Generate the pages of some or all folders
            %   PAGES = run(OBJ) writes the reference pages of all folders
            %   (docgen.pageTree); run(OBJ,FOLDERS) only those of the
            %   folders in the cell array FOLDERS (e.g. {'forms'}). PAGES
            %   lists the written files. The links are checked by
            %   checkDocs once all pages are written.
            rows = obj.Context.Pages;
            if nargin < 2 || isempty(folders)
                folders = unique(rows(:, 4))';  % sorted
            end
            obj.buildSubMap(rows);
            pages = cell(1, 0);
            for f = 1:numel(folders)
                folder = folders{f};
                count = 0;
                empty = cell(1, 0);
                for r = find(strcmp(rows(:, 4), folder))'
                    [main, page, src] = rows{r, [5 1 2]};
                    html = obj.render(fullfile(obj.Context.Toolbox, src), main, folder);
                    outDir = fullfile(obj.Context.DocFolder, main);
                    if ~isfolder(outDir)
                        mkdir(outDir);
                    end
                    file = fullfile(outDir, page);
                    docgen.Text.writeText(file, html);
                    pages{end+1} = file; %#ok<AGROW>
                    count = count + 1;
                    if contains(html, '(no help text)')
                        empty{end+1} = page; %#ok<AGROW>
                    end
                end
                fprintf('folder %s: %d pages, empty desc: [%s]\n', folder, ...
                    count, strjoin(empty, ', '));
            end
        end

        function html = render(obj, file, pageDir, folder)
            %render - HTML of the reference page of one source file
            %   HTML = render(OBJ,FILE,PAGEDIR,FOLDER) returns the page of
            %   the class or function in FILE for the page folder PAGEDIR
            %   (a main folder such as 'elements'); FOLDER is the source
            %   folder relative to the toolbox (e.g. 'elements/mixed').
            pg = obj.pageInfo(file, pageDir, folder);
            E = @docgen.Text.escape;
            parts = {['<h1>' pg.stem '</h1>'], ...
                ['<div class="purpose_line"><span class="purpose">' E(pg.h1) ...
                '</span><span class="reftype">' refType(pg) '</span></div>']};
            details = {};
            if pg.isClass
                [body, details] = obj.classBody(pg);
            else
                body = functionBody(pg);
            end
            parts = [parts body obj.examplesSection(pg) obj.tipsSection(pg) ...
                details obj.seeAlsoSection(pg)];
            html = docgen.page(obj.Context, [pg.stem ' - ' pg.h1], crumbs(pg), ...
                parts, Root='../', BodyClass='refpage');
        end

        function printNotes(obj)
            %printNotes - Print what may need a look in the help texts
            %   printNotes(OBJ) prints the class members without any
            %   description, the ambiguous names (not linked) and the
            %   notes on generated Syntax sections, except the common
            %   case of a class without an explicit constructor.
            if isempty(obj.EmptyDesc)
                fprintf('all members have descriptions\n');
            else
                rows = unique(strcat(obj.EmptyDesc(:, 1), {': '}, ...
                    obj.EmptyDesc(:, 2), {' '}, obj.EmptyDesc(:, 3)));
                fprintf('members without description: %d\n', numel(rows));
                fprintf('  %s\n', rows{:});
            end
            amb = sort(obj.Context.Lookup.Ambiguous);
            if isempty(amb)
                fprintf('no ambiguous names\n');
            else
                fprintf('ambiguous names (not linked): %s\n', strjoin(amb, ', '));
            end
            notes = obj.SyntaxLog(~startsWith(obj.SyntaxLog(:, 2), ...
                'no explicit constructor'), :);
            if ~isempty(notes)
                notes = unique(strcat(notes(:, 1), {': '}, notes(:, 2)), 'stable');
                fprintf('Syntax notes: %d\n', numel(notes));
                fprintf('  %s\n', notes{:});
            end
        end

    end

    methods (Access = private)

        function pg = pageInfo(obj, file, pageDir, folder)
            % everything the sections of a reference page are made of:
            % the help text split into its parts, the class information
            % and the calls of the Syntax block
            T = docgen.Text;
            [pg.kind, pg.sig, pg.hb, pg.lines] = docgen.Help.block(file);
            pg.srcRel = obj.relSrc(file);
            pg.stem = T.stem(file);
            pg.pageDir = pageDir;
            pg.folder = folder;
            [pg.h1, h1Lines] = summaryLine(pg.hb, pg.stem);

            % ---- description lines without summary line ----
            [desc, pg.lists, pg.examples] = docgen.Help.splitHelp(pg.hb);
            while ~isempty(desc) && isempty(strtrim(desc{1}))
                desc(1) = [];
            end
            if ~isempty(desc) && (startsWith(lower(strtrim(desc{1})), lower(pg.stem)) ...
                    || strcmp(strtrim(desc{1}), pg.h1))
                desc = desc(h1Lines+1:end);
            end
            while ~isempty(desc) && isempty(strtrim(desc{end}))
                desc(end) = [];
            end

            pg.isClass = strcmp(pg.kind, 'class');
            pg.abstract = false;
            pg.cname = pg.stem;
            pg.qname = pg.stem;
            pg.cls = [];
            if pg.isClass
                [~, headerName] = docgen.Help.classHeader(pg.lines);
                if ~isempty(headerName)
                    pg.cname = headerName;
                end
                pg.qname = docgen.qualifiedName(pg.srcRel, pg.cname);
                pg.cls = obj.classInfo(pg.srcRel, pg.qname, pg.lines);
                pg.abstract = pg.cls.abstract;
            end
            [pg.calls, pg.staticOnly] = obj.syntaxCalls(pg);
            pg.args = docgen.Help.argTokens(pg.hb, pg.cname);
            [pg.intro, pg.items, pg.ids] = obj.describe(desc, pg.cname, pg.qname, ...
                pg.args, pageDir, pg.srcRel, 'syntax_');
        end

        function [calls, staticOnly] = syntaxCalls(obj, pg)
            % calls shown in the Syntax block, from the help text, else
            % from the constructor or function signature; STATICONLY is
            % true for a class with only static methods
            calls = cell(1, 0);
            staticOnly = false;
            if ~pg.isClass
                if strcmp(pg.kind, 'function')
                    calls = docgen.Help.callLines(pg.hb, pg.stem);
                end
                if isempty(calls)
                    tok = regexp(strtrim(docgen.Help.stripComment(pg.sig)), ...
                        '^function\s+(.*?\s*=\s*)?(\w+)\s*(\(.*\))?', 'tokens', 'once');
                    if isempty(tok)
                        calls = {pg.sig};
                    else
                        out = regexprep(tok{1}, '\s*=\s*$', '');
                        calls = {docgen.Help.normCall(out, tok{3}, tok{2})};
                    end
                end
                return
            end
            [cname, qname, srcRel] = deal(pg.cname, pg.qname, pg.srcRel);
            ctor = docgen.Help.findFunction(pg.lines, cname);
            staticOnly = ~isempty(pg.cls.methods) && all([pg.cls.methods.static]) ...
                && isempty(ctor);
            if staticOnly || pg.abstract
                return
            end
            calls = docgen.Help.callLines(pg.hb, cname, qname);
            if isempty(ctor)
                % no own constructor: the superclass constructor is used
                for c = obj.inheritedCalls(pg.lines, qname, {})
                    if ~any(strcmp(calls, c{1}))
                        calls{end+1} = c{1}; %#ok<AGROW>
                    end
                end
            end
            if ~isempty(calls)
                return
            end
            if isempty(ctor)
                obj.SyntaxLog(end+1, :) = {srcRel, ['no explicit ' ...
                    'constructor; showing `obj = ' qname '()`']};
                calls = {['obj = ' qname '()']};
                return
            end
            calls = docgen.Help.callLines(ctor.help, cname, qname);
            if isempty(calls)
                inputs = strsplit(ctor.inputs, ',');
                if ~isempty(ctor.inputs) && all(ismember( ...
                        cellfun(@strtrim, inputs, 'UniformOutput', false), ...
                        {'varargin', 'varargout'}))
                    obj.SyntaxLog(end+1, :) = {srcRel, ['constructor ' ...
                        'takes only `' ctor.inputs '`; help does not ' ...
                        'describe the arguments; showing call from signature']};
                end
                out = 'obj';
                if ~isempty(ctor.outputs)
                    outs = strsplit(ctor.outputs, ',');
                    out = strtrim(outs{1});
                end
                calls = {docgen.Help.normCall(out, ['(' ctor.inputs ')'], qname)};
            end
        end

        function [parts, details] = classBody(obj, pg)
            % MathWorks class page: Description, Creation (Syntax +
            % Description), Subclasses, Properties, Methods, Events; DETAILS
            % holds the Method Details section
            E = @docgen.Text.escape;
            cls = pg.cls;
            [pageDir, srcRel] = deal(pg.pageDir, pg.srcRel);
            parts = {};
            if ~isempty(cls.super)
                parts{end+1} = ['<p class="class_meta"><span class="label">' ...
                    'Superclasses:</span> ' ...
                    obj.Context.Lookup.linkNames(cls.super, pageDir, true) '</p>'];
            end
            subclasses = {};
            if pg.abstract
                subclasses = obj.subclassTable(pg.stem, srcRel, pageDir);
                if isempty(subclasses)  % no section to link to
                    subLink = 'subclasses';
                else
                    subLink = '<a href="#subclasses">subclasses</a>';
                end
                parts{end+1} = ['<div class="note">This is an abstract class and ' ...
                    'cannot be instantiated directly. Use one of its ' ...
                    subLink '.</div>'];
            end

            noCreation = pg.abstract || pg.staticOnly;
            if ~isempty(pg.intro) || (~isempty(pg.items) && noCreation)
                parts{end+1} = '<h2 id="description">Description</h2>';
                parts{end+1} = pg.intro;
                if noCreation  % no Creation section to hold the items
                    parts = [parts pg.items];
                end
            elseif isempty(pg.h1) && isempty(pg.items)
                parts{end+1} = ['<h2 id="description">Description</h2>' ...
                    '<p><em>(no help text)</em></p>'];
            end
            if pg.staticOnly
                parts{end+1} = ['<h2 id="creation">Creation</h2><p><code>' ...
                    E(pg.qname) '</code> has only static methods; call them as ' ...
                    '<code>' E(pg.qname) '.<var>method</var>(...)</code>.</p>'];
            elseif ~pg.abstract
                parts{end+1} = '<h2 id="creation">Creation</h2>';
                parts{end+1} = ['<h3>Syntax</h3>' syntaxHtml(pg.calls, pg.ids)];
                if ~isempty(pg.items)
                    parts{end+1} = '<h3>Description</h3>';
                    parts = [parts pg.items];
                end
            end
            parts = [parts subclasses];

            inheritNote = obj.inheritNote(cls, pageDir, srcRel);
            if ~isempty(cls.props)
                parts = [parts obj.propertyTable(pg)];
                if ~isempty(inheritNote)
                    parts{end+1} = inheritNote;
                    inheritNote = '';
                end
            end
            detailMap = obj.methodDetails(pg.lines, cls, pageDir, srcRel);
            if ~isempty(cls.methods)
                parts = [parts obj.methodTable(pg, detailMap)];
            end
            if ~isempty(inheritNote)
                parts{end+1} = inheritNote;
            end
            if ~isempty(cls.events)
                parts{end+1} = ['<h2>Events</h2><table class="members"><tr>' ...
                    '<th>Name</th><th>Description</th></tr>'];
                for ev = cls.events
                    parts{end+1} = ['<tr><td><code>' E(ev.name) '</code></td><td>' ...
                        E(docgen.stripName(ev.name, ev.desc)) '</td></tr>']; %#ok<AGROW>
                end
                parts{end+1} = '</table>';
            end

            details = {};
            if numEntries(detailMap) > 0
                details{1} = '<h2 id="method_details">Method Details</h2>';
                for m = cls.methods
                    if isKey(detailMap, m.name)
                        details{end+1} = detailMap{m.name}; %#ok<AGROW>
                    end
                end
            end
        end

        function note = inheritNote(obj, cls, pageDir, srcRel)
            % 'Inherited properties and methods are listed on ...' with
            % links to the toolbox superclasses ('' if there are none)
            supers = cellfun(@strtrim, strsplit(cls.super, ','), ...
                'UniformOutput', false);
            supLinks = {};
            for s = supers
                href = obj.Context.Lookup.pathFor(s{1}, pageDir, srcRel);
                if ~isempty(href)
                    supLinks{end+1} = ['<a href="' href '"><code>' ...
                        docgen.Text.escape(s{1}) '</code></a>']; %#ok<AGROW>
                end
            end
            note = '';
            if ~isempty(supLinks)
                plural = '';
                if numel(supLinks) > 1
                    plural = 's';
                end
                note = ['<p class="inherit">Inherited properties and ' ...
                    'methods are listed on the superclass page' plural ' ' ...
                    strjoin(supLinks, ', ') '.</p>'];
            end
        end

        function parts = propertyTable(obj, pg)
            % Properties section of a class page
            E = @docgen.Text.escape;
            parts = {['<h2 id="properties">Properties</h2>' ...
                '<table class="members"><tr><th>Name</th>' ...
                '<th>Description</th></tr>']};
            for p = pg.cls.props
                % help member list first: declaration comments are often
                % code notes
                text = '';
                if isKey(pg.lists, p.name)
                    text = pg.lists{p.name};
                end
                if isempty(text)
                    text = strtrim(p.meta.Description);
                end
                text = docgen.stripName(p.name, text);
                if isempty(text)
                    obj.EmptyDesc(end+1, :) = {pg.srcRel, 'property', p.name};
                end
                parts{end+1} = ['<tr><td><code>' E(p.name) '</code>' ...
                    tagHtml(p.tags) '</td><td>' E(text) '</td></tr>']; %#ok<AGROW>
            end
            parts{end+1} = '</table>';
        end

        function parts = methodTable(obj, pg, detailMap)
            % Methods section of a class page; a method links to its
            % Method Details entry (DETAILMAP), else to its own page
            E = @docgen.Text.escape;
            parts = {['<h2 id="methods">Methods</h2>' ...
                '<table class="members"><tr><th>Name</th>' ...
                '<th>Description</th></tr>']};
            for m = pg.cls.methods
                % a method's own summary line is the fullest text; a first
                % comment that does not name the method is a code comment
                d = strtrim(m.desc);
                nameParts = strsplit(m.name, '.');
                h1Ok = ~isempty(regexpi(d, ['^' regexptranslate('escape', ...
                    nameParts{end}) '\>\s*[-' char(8211) ...
                    char(8212) ':]?'], 'once'));
                text = '';
                if h1Ok
                    text = d;
                end
                if isempty(text) && isKey(pg.lists, m.name)
                    text = pg.lists{m.name};
                end
                if isempty(text)
                    text = d;
                end
                text = docgen.stripName(m.name, text);
                isCtor = strcmp(m.name, pg.cname);
                if isempty(text) && isCtor
                    if pg.abstract
                        text = 'Initialize the object (called by subclass constructors)';
                    else
                        text = ['Create a ' pg.cname ' object'];
                    end
                end
                if isempty(text)
                    obj.EmptyDesc(end+1, :) = {pg.srcRel, 'method', m.name};
                end
                if isCtor && ~pg.abstract && ~pg.staticOnly
                    href = '#creation';
                elseif isKey(detailMap, m.name)
                    href = ['#method_' m.name];
                else
                    href = obj.Context.Lookup.pathFor(m.name, pg.pageDir, pg.srcRel);
                end
                if isempty(href)
                    cellHtml = ['<code>' E(m.name) '</code>'];
                else
                    cellHtml = ['<a href="' href '"><code>' E(m.name) '</code></a>'];
                end
                cellHtml = [cellHtml tagHtml(m.tags)]; %#ok<AGROW>
                if contains(cellHtml, 'abstract')  % the tag already says it
                    text = regexprep(text, '\s*\(abstract\)', '');
                end
                parts{end+1} = ['<tr><td>' cellHtml '</td><td>' E(text) ...
                    '</td></tr>']; %#ok<AGROW>
            end
            parts{end+1} = '</table>';
        end

        function parts = examplesSection(obj, pg)
            % Examples: help-text examples and the example pages that use
            % the class or function
            parts = {};
            related = obj.relatedExamples(pg.stem);
            if isempty(pg.examples) && isempty(related)
                return
            end
            parts{1} = '<h2 id="examples">Examples</h2>';
            if ~isempty(pg.examples)
                d = min(cellfun(@docgen.Text.indent, pg.examples));
                ex = cellfun(@(l) l(d+1:end), pg.examples, 'UniformOutput', false);
                parts{end+1} = ['<pre class="codeinput">' ...
                    docgen.Text.escape(strjoin(ex, newline)) '</pre>'];
            end
            if ~isempty(related)
                parts{end+1} = ['<h3>Related examples</h3><ul class="related">' ...
                    [related{:}] '</ul>'];
            end
        end

        function parts = tipsSection(obj, pg)
            % Tips: the notes after 'See also' in the help text
            parts = {};
            [~, notes] = docgen.Help.seeAlso(pg.hb);
            if ~isempty(notes)
                li = cellfun(@(n) ['<li>' obj.inlineFmt(regexprep(n, ...
                    '^notes?:\s*', '', 'ignorecase'), pg.args, pg.pageDir, ...
                    pg.srcRel) '</li>'], notes, 'UniformOutput', false);
                parts{1} = ['<h2 id="tips">Tips</h2><ul class="tips">' ...
                    [li{:}] '</ul>'];
            end
        end

        function parts = seeAlsoSection(obj, pg)
            % See Also: the 'See also' names and the overview page of the
            % folder
            seeNames = unique(docgen.Help.seeAlso(pg.hb), 'stable');
            links = {};
            for k = 1:numel(seeNames)
                n = seeNames{k};
                href = obj.Context.Lookup.pathFor(n, pg.pageDir, pg.srcRel);
                if ~isempty(href)
                    links{end+1} = ['<a href="' href '"><code>' n '</code></a>']; %#ok<AGROW>
                else
                    % MATLAB function: link works in the MATLAB Help browser
                    links{end+1} = ['<a href="matlab:doc(''' n ''')"><code>' ...
                        n '</code></a>']; %#ok<AGROW>
                end
            end
            [~, base] = fileparts(pg.folder);
            parts = {'<h2 id="seealso">See Also</h2>'};
            if ~isempty(links)
                parts{end+1} = ['<p class="seealso">' strjoin(links, ' | ') '</p>'];
            end
            parts{end+1} = ['<h3>Topics</h3><ul class="related"><li><a href="' ...
                base '_overview.html">' docgen.Text.escape( ...
                obj.Context.overviewTitle(pg.folder)) '</a></li></ul>'];
        end

        function src = relSrc(obj, file)
            src = strrep(file(numel(obj.Context.Toolbox)+2:end), filesep, '/');
        end

        function cls = classInfo(obj, srcRel, qname, lines)
            % cached member information of a class
            if ~isKey(obj.Classes, srcRel)
                cls = docgen.classMembers(qname, lines);
                mc = meta.class.fromName(qname);
                evs = mc.EventList;
                evs = evs(arrayfun(@(e) e.DefiningClass == mc && ~e.Hidden, evs));
                cls.events = struct('name', {evs.Name}, 'desc', {evs.Description});
                obj.Classes{srcRel} = cls;
            end
            cls = obj.Classes{srcRel};
        end

        function buildSubMap(obj, rows)
            % class name -> direct toolbox subclasses (name, src, h1, abstract)
            info = configureDictionary("string", "cell");
            names = {};
            for r = find(strcmp(rows(:, 3), 'class'))'
                src = rows{r, 2};
                lines = docgen.Text.readLines(fullfile(obj.Context.Toolbox, src));
                [~, cname, supers] = docgen.Help.classHeader(lines);
                name = cname;
                if isempty(name)
                    name = docgen.Text.stem(src);
                end
                [~, h1] = docgen.Help.fileSummary(lines);
                cls = obj.classInfo(src, docgen.qualifiedName(src, name), lines);
                ent = struct('name', name, 'src', src, 'h1', h1, ...
                    'abstract', cls.abstract, 'supers', {supers});
                if isKey(info, name)
                    info{name} = [info{name} ent];
                else
                    info{name} = ent;
                    names{end+1} = name; %#ok<AGROW>
                end
            end
            obj.SubMap = configureDictionary("string", "cell");
            for i = 1:numel(names)
                name = names{i};
                for ent = info{name}
                    for s = ent.supers
                        key = strtrim(regexprep(s{1}, '.*\.', ''));
                        if isKey(info, key) && ~strcmp(key, name)
                            if ~isKey(obj.SubMap, key)
                                obj.SubMap{key} = struct('name', {}, 'ents', {});
                            end
                            bucket = obj.SubMap{key};
                            if ~any(strcmp({bucket.name}, name))
                                bucket(end+1) = struct('name', name, ...
                                    'ents', info{name}); %#ok<AGROW>
                                obj.SubMap{key} = bucket;
                            end
                        end
                    end
                end
            end
        end

        function parts = subclassTable(obj, stem, srcRel, pageDir)
            % Subclasses section of an abstract class
            parts = {};
            if ~isKey(obj.SubMap, stem)
                bucket = struct('name', {}, 'ents', {});
            else
                bucket = obj.SubMap{stem};
            end
            rows = struct('name', {}, 'ent', {});
            seen = {};
            for b = bucket
                for ent = b.ents
                    if strcmp(ent.src, srcRel) || any(strcmp(seen, b.name))
                        continue
                    end
                    seen{end+1} = b.name; %#ok<AGROW>
                    rows(end+1) = struct('name', b.name, 'ent', ent); %#ok<AGROW>
                end
            end
            if isempty(rows)
                obj.SyntaxLog(end+1, :) = {srcRel, ['abstract class has no toolbox ' ...
                    'subclasses; Subclasses section omitted']};
                return
            end
            [~, order] = sort(lower({rows.name}));
            E = @docgen.Text.escape;
            parts{1} = ['<h2 id="subclasses">Subclasses</h2>' ...
                '<table class="members"><tr><th>Name</th><th>Description</th></tr>'];
            for r = rows(order)
                shown = docgen.qualifiedName(r.ent.src, r.name);
                href = obj.Context.Lookup.pathFor(shown, pageDir);
                if isempty(href)
                    cellHtml = ['<code>' E(shown) '</code>'];
                else
                    cellHtml = ['<a href="' href '"><code>' E(shown) '</code></a>'];
                end
                text = docgen.stripName(shown, r.ent.h1);
                if isempty(text)
                    text = '(no description)';
                end
                tag = '';
                if r.ent.abstract
                    tag = ' <span class="tag">abstract</span>';
                end
                parts{end+1} = ['<tr><td>' cellHtml tag '</td><td>' E(text) ...
                    '</td></tr>']; %#ok<AGROW>
            end
            parts{end+1} = '</table>';
        end

        function calls = inheritedCalls(obj, lines, display, seen)
            % call lines of the nearest toolbox superclass that documents
            % its constructor, rewritten for class DISPLAY
            calls = {};
            [~, ~, supers] = docgen.Help.classHeader(lines);
            for s = supers
                parts = strsplit(s{1}, '.');
                key = strtrim(parts{end});
                ents = obj.Context.Lookup.entries(key);
                ents = ents(endsWith({ents.src}, '.m'));
                if numel(ents) ~= 1 || any(strcmp(seen, key))
                    continue
                end
                seen{end+1} = key; %#ok<AGROW>
                [~, ~, shb, slines] = docgen.Help.block( ...
                    fullfile(obj.Context.Toolbox, ents.src));
                calls = docgen.Help.callLines(shb, key, display);
                if isempty(calls)
                    ctor = docgen.Help.findFunction(slines, key);
                    if ~isempty(ctor)
                        calls = docgen.Help.callLines(ctor.help, key, display);
                        if isempty(calls) && ~isempty(ctor.inputs)
                            outs = strsplit(ctor.outputs, ',');
                            out = strtrim(outs{1});
                            if isempty(out)
                                out = 'obj';
                            end
                            calls = {docgen.Help.normCall(out, ['(' ctor.inputs ')'], display)};
                        end
                    else
                        calls = obj.inheritedCalls(slines, display, seen);
                    end
                end
                if ~isempty(calls)
                    return
                end
            end
        end

        function [intro, items, ids] = describe(obj, lines, name, display, args, ...
                pageDir, selfSrc, idPrefix)
            % Description in MathWorks style: INTRO is the HTML of the
            % blocks before the first call paragraph, ITEMS the HTML of one
            % description item per call paragraph ('obj = RTk(degree,mesh)
            % creates ...') with the paragraphs that follow it, IDS maps
            % normalised calls to the anchor ids of the items.
            rx = ['^\s*' docgen.Help.callRe(name)];
            introParts = {};
            itemParts = {};
            ids = configureDictionary("string", "cell");
            blocks = textBlocks(lines);
            for b = 1:size(blocks, 1)
                [kind, text] = blocks{b, :};
                tok = {};
                if strcmp(kind, 'p')
                    [tok, e] = regexpi(text, rx, 'tokens', 'end', 'once');
                end
                if ~isempty(tok)
                    out = strtrim(regexprep(tok{1}, '[ =]+$', ''));
                    call = docgen.Help.normCall(out, tok{3}, display);
                    aid = sprintf('%s%d', idPrefix, numel(itemParts) + 1);
                    key = strrep(call, ' ', '');
                    if ~isKey(ids, key)
                        ids{key} = aid;
                    end
                    itemParts{end+1} = {['<div class="description_item" id="' aid ...
                        '"><p><code class="sig">' docgen.Text.escape(call) '</code>' ...
                        obj.inlineFmt(text(e+1:end), args, pageDir, selfSrc) '</p>']}; %#ok<AGROW>
                    continue
                end
                if strcmp(kind, 'pre')
                    h = ['<pre class="helptext">' docgen.Text.escape(text) '</pre>'];
                else
                    h = ['<p>' obj.inlineFmt(text, args, pageDir, selfSrc) '</p>'];
                end
                if isempty(itemParts)
                    introParts{end+1} = h; %#ok<AGROW>
                else
                    itemParts{end}{end+1} = h;
                end
            end
            intro = strjoin(introParts, newline);
            items = cellfun(@(it) [strjoin(it, newline) '</div>'], itemParts, ...
                'UniformOutput', false);
        end

        function out = inlineFmt(obj, text, args, pageDir, selfSrc)
            % HTML of one prose paragraph: string literals and argument
            % names become <code>, toolbox class/function names links
            rx = ['"[^"\n]{1,60}"|(?<!\w)''[^''\s][^''\n]{0,58}''(?!\w)' ...
                '|[A-Za-z_]\w*(\{[^}\s]{1,8}\})?'];
            [starts, ends] = regexp(text, rx, 'start', 'end');
            E = @docgen.Text.escape;
            out = '';
            pos = 1;
            for k = 1:numel(starts)
                out = [out E(text(pos:starts(k)-1))]; %#ok<AGROW>
                tok = text(starts(k):ends(k));
                after = text(ends(k)+1:end);
                pos = ends(k) + 1;
                if tok(1) == '"' || tok(1) == ''''
                    out = [out '<code>' E(tok) '</code>']; %#ok<AGROW>
                    continue
                end
                brace = find(tok == '{', 1);
                if ~isempty(brace)  % indexed argument, e.g. OP{1}
                    base = tok(1:brace-1);
                    if any(strcmp(args, base)) && isArgLike(base)
                        out = [out '<code>' E(docgen.Help.lowerCaps(base)) ...
                            E(tok(brace:end)) '</code>']; %#ok<AGROW>
                    else
                        out = [out E(tok)]; %#ok<AGROW>
                    end
                    continue
                end
                if any(strcmp(args, tok)) && isArgLike(tok)
                    % 'A mesh ...' / 'I ...' at a sentence start are words
                    if ~(any(strcmp(tok, {'A', 'I'})) && ...
                            ~isempty(regexp(after, '^\s+[a-z]', 'once')))
                        out = [out '<code>' E(docgen.Help.lowerCaps(tok)) '</code>']; %#ok<AGROW>
                        continue
                    end
                end
                if ~isempty(regexp(tok, '[a-z]', 'once')) && ...
                        ~isempty(regexp(tok, '[A-Z]', 'once')) && ...
                        numel(tok) > 1 && ~startsWith(after, '(')
                    href = obj.Context.Lookup.pathFor(tok, pageDir, selfSrc);
                    if ~isempty(href) && ~startsWith(href, '../examples/')
                        out = [out '<a href="' href '"><code>' tok '</code></a>']; %#ok<AGROW>
                        continue
                    end
                end
                out = [out E(tok)]; %#ok<AGROW>
            end
            out = [out E(text(pos:end))];
        end

        function details = methodDetails(obj, lines, cls, pageDir, srcRel)
            % method name -> HTML of its Method Details entry, for methods
            % whose help has more than the summary line
            T = docgen.Text;
            details = configureDictionary("string", "cell");
            for m = cls.methods
                n = m.name;
                if contains(n, '.')
                    continue  % set./get. property accessors
                end
                info = docgen.Help.findFunction(lines, n);
                if isempty(info) || numel(info.help) < 2
                    continue
                end
                mhb = info.help;
                h1 = docgen.stripName(n, strtrim(docgen.Help.uncomment(mhb{1})));
                desc = docgen.Help.splitHelp(mhb);
                body = desc(2:end);
                % a code comment right below the help ('% parse input') is
                % not indented like the help text: the help ends there
                for i = 1:numel(body)
                    if ~isempty(strtrim(body{i})) && T.indent(body{i}) < 2
                        body = body(1:i-1);
                        break
                    end
                end
                if all(cellfun(@(l) isempty(strtrim(l)), body))
                    continue
                end
                args = docgen.Help.argTokens(mhb, n);
                [intro, items, ids] = obj.describe(body, n, n, args, pageDir, ...
                    srcRel, [n '_syntax_']);
                calls = docgen.Help.callLines(mhb, n);
                h = {['<div class="method_detail" id="method_' n '">'], ...
                    ['<h3><code>' T.escape(n) '</code></h3>'], ...
                    ['<p class="purpose">' T.escape(h1) '</p>']};
                if ~isempty(calls)
                    h{end+1} = syntaxHtml(calls, ids); %#ok<AGROW>
                end
                if ~isempty(intro)
                    h{end+1} = intro; %#ok<AGROW>
                end
                h = [h items, {'</div>'}]; %#ok<AGROW>
                details{n} = strjoin(h, newline);
            end
        end

        function items = relatedExamples(obj, stem)
            % list items of the example pages that use STEM, most specific
            % examples (fewest referenced names) first
            items = {};
            ex = obj.Context.Examples;
            hits = ex(arrayfun(@(e) any(strcmp(e.refs, stem)), ex));
            if isempty(hits)
                return
            end
            nrefs = arrayfun(@(e) numel(e.refs), hits);
            files = lower({hits.file});
            [~, o1] = sort(files);  % secondary key first (stable sort)
            [~, o2] = sort(nrefs(o1));
            hits = hits(o1(o2));
            for e = reshape(hits(1:min(end, obj.MaxRelated)), 1, [])
                if ~obj.Context.hasExamplePage(e.name)
                    continue
                end
                items{end+1} = ['<li><a href="../examples/' e.name '.html">' ...
                    docgen.exampleTitle(e.title, e.name) '</a> <span class="exname">' ...
                    '(<code>' e.name '</code>)</span></li>']; %#ok<AGROW>
            end
            more = numel(hits) - obj.MaxRelated;
            if more > 0
                items{end+1} = sprintf(['<li>%d more in <a href="../examples/' ...
                    'examples_overview.html">Examples</a></li>'], more);
            end
        end

    end
end

% ---------------------------------------------------------------------------

function tf = isArgLike(tok)
% identifiers that are safe to show as code in prose: all-caps names
% (DEGREE), single capitals and mixed-case identifiers (varStr, s4l)
if isscalar(tok)
    tf = docgen.Text.isUpper(tok);
else
    tf = docgen.Text.isUpper(tok) || ~isempty(regexp(tok, '[a-z][A-Z]|\d|_', 'once'));
end
end

function h = tagHtml(tags)
h = '';
if ~isempty(tags)
    h = sprintf(' <span class="tag">%s</span>', tags{:});
end
end

function [h1, n] = summaryLine(hb, stem)
% summary line of the help lines HB without a leading repetition of STEM;
% a summary wrapped onto a second line is joined (N = 2 lines)
h1 = '';
if ~isempty(hb)
    h1 = strtrim(docgen.Help.uncomment(hb{1}));
end
n = 1;
if numel(hb) > 1 && ~endsWith(h1, '.')
    l1 = strtrim(docgen.Help.uncomment(hb{2}));
    l2 = '';
    if numel(hb) > 2
        l2 = strtrim(docgen.Help.uncomment(hb{3}));
    end
    if ~isempty(l1) && isempty(l2) && isempty(regexp(l1, '[=(]', 'once')) ...
            && isempty(regexpi(hb{2}, docgen.Help.HeadRe, 'once'))
        h1 = [h1 ' ' l1];
        n = 2;
    end
end
h1 = docgen.stripName(stem, h1);
end

function c = crumbs(pg)
% breadcrumb after the start page: main folder > subfolder > name
main = pg.pageDir;
c = {main, [main '_overview.html']};
if ~strcmp(pg.folder, main)
    [~, base] = fileparts(pg.folder);
    c(end+1, :) = {base, [base '_overview.html']};
end
c(end+1, :) = {pg.stem, ''};
end

function t = refType(pg)
% page type shown next to the summary line
if pg.isClass && pg.abstract
    t = 'Abstract class';
elseif pg.isClass
    t = 'Class';
elseif strcmp(pg.kind, 'function')
    t = 'Function';
else
    t = 'Script';
end
end

function parts = functionBody(pg)
% MathWorks function page: Syntax, Description
parts = {['<h2 id="syntax">Syntax</h2>' syntaxHtml(pg.calls, pg.ids)]};
if ~isempty(pg.intro) || ~isempty(pg.items)
    parts{end+1} = '<h2 id="description">Description</h2>';
    if ~isempty(pg.intro)
        parts{end+1} = pg.intro;
    end
    parts = [parts pg.items];
elseif isempty(pg.h1)
    parts{end+1} = ['<h2 id="description">Description</h2>' ...
        '<p><em>(no help text)</em></p>'];
end
end

function h = syntaxHtml(calls, ids)
% Syntax block: one line per call, linked to its description item
rows = cell(1, numel(calls));
for k = 1:numel(calls)
    c = calls{k};
    code = docgen.Text.escape(c);
    key = strrep(c, ' ', '');
    if isKey(ids, key)
        rows{k} = ['<div class="syntax_line"><a href="#' ids{key} '"><code>' ...
            code '</code></a></div>'];
    else
        rows{k} = ['<div class="syntax_line"><code>' code '</code></div>'];
    end
end
h = ['<div class="syntax_block">' [rows{:}] '</div>'];
end

function blocks = textBlocks(lines)
% Split help lines (without '%') into blocks {kind, text} with kind 'p'
% or 'pre'. Lines indented deeper than the help text (ASCII sketches,
% option lists, code) are kept verbatim as 'pre'.
T = docgen.Text;
blocks = cell(0, 2);
nonblank = lines(~cellfun(@(l) isempty(strtrim(l)), lines));
if isempty(nonblank)
    return
end
base = T.indent(nonblank{1});  % indent of the help text
blank = false;
prev = base;
for k = 1:numel(lines)
    l = lines{k};
    if isempty(strtrim(l))
        blank = true;
        continue
    end
    cur = T.indent(l);
    deep = cur >= base + 2;
    if deep
        kind = 'pre';
    else
        kind = 'p';
    end
    % an outdented line ('opts.name:') starts a new paragraph
    shift = strcmp(kind, 'p') && (cur < base) ~= (prev < base);
    prev = cur;
    if ~isempty(blocks) && strcmp(blocks{end, 1}, kind) && ~shift && ...
            (~blank || deep)
        if deep
            if blank
                sep = [newline newline];
            else
                sep = newline;
            end
            body = deblank(l);
        else
            sep = ' ';
            if ~isempty(regexp(blocks{end, 2}, '[\w)]-$', 'once'))
                sep = '';  % word hyphenated at the line end: 'H(div)-'
            end
            body = strtrim(l);
        end
        blocks{end, 2} = [blocks{end, 2} sep body];
    elseif deep
        blocks(end+1, :) = {kind, deblank(l)}; %#ok<AGROW>
    else
        blocks(end+1, :) = {kind, strtrim(l)}; %#ok<AGROW>
    end
    blank = false;
end
for b = 1:size(blocks, 1)  % dedent pre blocks
    if strcmp(blocks{b, 1}, 'pre')
        ls = strsplit(blocks{b, 2}, newline, 'CollapseDelimiters', false);
        ind = cellfun(@docgen.Text.indent, ls(~cellfun(@(x) isempty(strtrim(x)), ls)));
        d = 0;
        if ~isempty(ind)
            d = min(ind);
        end
        ls = cellfun(@(x) x(min(d, numel(x))+1:end), ls, 'UniformOutput', false);
        blocks{b, 2} = strjoin(ls, newline);
    end
end
end