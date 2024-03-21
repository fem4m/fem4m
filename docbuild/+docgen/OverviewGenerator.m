classdef OverviewGenerator < handle
    %OverviewGenerator - Generate the overview pages of the documentation
    %   OBJ = docgen.OverviewGenerator(CTX) prepares the overview pages of
    %   the toolbox described by the docgen.Context CTX, to be written to
    %   CTX.DocFolder: one page per main folder and per subfolder, the
    %   start page index.html and getting_started.html. The page texts are
    %   Markdown drafts in
    %   docbuild/overviews (<main>.md, <main>_<sub>.md, index.md,
    %   getting_started.md); a subfolder without a draft gets a page built
    %   from its Contents.m and the summary lines of its files. Known
    %   class, function and example names are linked to their pages.
    %
    %   The drafts may contain these markers on a line of their own:
    %       {{SUBFOLDERS}}        - table of the subfolders of the folder
    %       {{MEMBERS}}           - table of the classes and functions
    %       {{MEMBERS:<sub>}}     - the same for the subfolder <sub>
    %       {{RELATED_EXAMPLES}}  - examples that use the folder's members
    %
    %   OverviewGenerator functions:
    %       run      - Write the overview pages
    %       mdToHtml - HTML of a Markdown draft
    %
    %   OverviewGenerator properties:
    %       Context - Folders and toolbox data (docgen.Context)
    %
    %   See also genOverviews, docgen.Context, docgen.RefGenerator

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = private)
        Context     % folders and toolbox data (docgen.Context)
    end

    properties (Access = private)
        Rows        % page tree as rows: page, src, kind, folder, main
        Members     % source folder -> cached member rows
    end

    methods

        function obj = OverviewGenerator(ctx)
            %OverviewGenerator - Prepare the generation of the overview pages
            obj.Context = ctx;
            obj.Members = configureDictionary("string", "cell");
            obj.Rows = ctx.Pages;
        end

        function pages = run(obj, mains)
            %run - Write the overview pages
            %   PAGES = run(OBJ) writes all overview pages and returns
            %   their files; run(OBJ,MAINS) only the pages of the main
            %   folders in the cell array MAINS (without index.html and
            %   getting_started.html).
            if nargin < 2
                mains = {};
            end
            T = docgen.Text;
            names = unique(obj.Rows(:, 5));
            pages = cell(1, 0);

            % ---- main folder overviews from drafts ----
            for i = 1:numel(names)
                main = names{i};
                if ~isempty(mains) && ~any(strcmp(mains, main))
                    continue
                end
                draft = obj.Context.draftFile(main);
                if ~isfile(draft)
                    fprintf('NO DRAFT for %s\n', main);
                    continue
                end
                body = obj.mdToHtml(T.readLines(draft), main, main);
                pages{end+1} = obj.write(fullfile(main, [main '_overview.html']), ...
                    docgen.page(obj.Context, [main ' overview'], {main, ''}, ...
                    body, Root='../')); %#ok<AGROW>
            end

            % ---- subfolder overviews: draft if present, else generated ----
            for i = 1:numel(names)
                main = names{i};
                if ~isempty(mains) && ~any(strcmp(mains, main))
                    continue
                end
                folders = unique(obj.Rows(strcmp(obj.Rows(:, 5), main), 4));
                for j = 1:numel(folders)
                    folder = folders{j};
                    if strcmp(folder, main)
                        continue
                    end
                    sub = baseName(folder);
                    draft = obj.Context.draftFile(folder);
                    if isfile(draft)
                        md = T.readLines(draft);  % hand-written page
                    else
                        % built on every run from Contents.m and the
                        % sources, so the page follows the help texts
                        md = T.splitLines(['# ' sub newline newline ...
                            obj.folderDesc(folder) newline newline ...
                            '## Contents' newline newline ...
                            '{{MEMBERS:' sub '}}' newline newline ...
                            '## Related examples' newline newline ...
                            '{{RELATED_EXAMPLES}}' newline]);
                    end
                    body = obj.mdToHtml(md, main, folder);
                    pages{end+1} = obj.write(fullfile(main, [sub '_overview.html']), ...
                        docgen.page(obj.Context, [folder ' overview'], {main, ...
                        [main '_overview.html']; sub, ''}, body, Root='../')); %#ok<AGROW>
                end
            end

            % ---- index.html and getting_started.html ----
            if isempty(mains)
                stems = {'index', 'getting_started'};
                titles = {obj.Context.Site, 'Getting Started'};
                crumbs = {cell(0, 2), {'Getting Started', ''}};
                for k = 1:2
                    md = T.readLines(fullfile(obj.Context.Build, 'overviews', ...
                        [stems{k} '.md']));
                    body = obj.mdToHtml(md, '', '');
                    pages{end+1} = obj.write([stems{k} '.html'], ...
                        docgen.page(obj.Context, titles{k}, crumbs{k}, body)); %#ok<AGROW>
                end
            end
        end

        function html = mdToHtml(obj, lines, pageDir, folder)
            %mdToHtml - HTML of a Markdown draft
            %   HTML = mdToHtml(OBJ,LINES,PAGEDIR,FOLDER) converts the lines
            %   of a draft (headings, paragraphs, lists, tables, code
            %   blocks, inline code/links/bold/italic) for a page in folder
            %   PAGEDIR ('' for the doc root). FOLDER is the source folder
            %   the markers refer to; with FOLDER '' markers are plain text.
            T = docgen.Text;
            marker = '^\s*\{\{(\w+(?::\w+)?)\}\}\s*$';
            ordered = '^\s*[0-9]+\.\s+(.*)';
            blockStart = '^(#{1,3} |```|\||- |\s*[0-9]+\.\s)';
            out = cell(1, 0);
            inUl = false;
            inOl = false;
            inPre = false;
            preBuf = cell(1, 0);
            i = 1;
            n = numel(lines);
            while i <= n
                ln = lines{i};
                if startsWith(ln, '```')
                    if inPre
                        out{end+1} = ['<pre>' T.escape(strjoin(preBuf, newline)) ...
                            '</pre>']; %#ok<AGROW>
                        preBuf = cell(1, 0);
                        inPre = false;
                    else
                        closeLists();
                        inPre = true;
                    end
                    i = i + 1;
                    continue
                end
                if inPre
                    preBuf{end+1} = ln; %#ok<AGROW>
                    i = i + 1;
                    continue
                end
                mm = regexp(ln, marker, 'tokens', 'once');
                if ~isempty(folder) && ~isempty(mm)
                    closeLists();
                    out = [out obj.markerHtml(mm{1}, pageDir, folder)]; %#ok<AGROW>
                    i = i + 1;
                    continue
                end
                hm = regexp(ln, '^(#{1,3}) (.*)', 'tokens', 'once');
                if ~isempty(hm)
                    closeLists();
                    lvl = sprintf('%d', numel(hm{1}));
                    out{end+1} = ['<h' lvl '>' obj.relink(T.escape(hm{2}), pageDir) ...
                        '</h' lvl '>']; %#ok<AGROW>
                    i = i + 1;
                    continue
                end
                if startsWith(ln, '|')
                    closeLists();
                    rows = cell(1, 0);
                    while i <= n && startsWith(lines{i}, '|')
                        s = regexprep(strtrim(lines{i}), '^\|+|\|+$', '');
                        rows{end+1} = cellfun(@strtrim, regexp(s, '\|', 'split'), ...
                            'UniformOutput', false); %#ok<AGROW>
                        i = i + 1;
                    end
                    out = [out obj.tableHtml(rows, pageDir)]; %#ok<AGROW>
                    continue
                end
                if startsWith(ln, '- ')
                    if inOl
                        out{end+1} = '</ol>'; %#ok<AGROW>
                        inOl = false;
                    end
                    if ~inUl
                        out{end+1} = '<ul>'; %#ok<AGROW>
                        inUl = true;
                    end
                    out{end+1} = ['<li>' obj.relink(T.escape(ln(3:end)), pageDir) ...
                        '</li>']; %#ok<AGROW>
                    i = i + 1;
                    continue
                end
                om = regexp(ln, ordered, 'tokens', 'once');
                if ~isempty(om)
                    if inUl
                        out{end+1} = '</ul>'; %#ok<AGROW>
                        inUl = false;
                    end
                    if ~inOl
                        out{end+1} = '<ol>'; %#ok<AGROW>
                        inOl = true;
                    end
                    out{end+1} = ['<li>' obj.relink(T.escape(om{1}), pageDir) ...
                        '</li>']; %#ok<AGROW>
                    i = i + 1;
                    continue
                end
                blank = isempty(strtrim(ln));
                if (inUl || inOl) && blank
                    closeLists();
                    i = i + 1;
                    continue
                end
                if ~blank
                    closeLists();
                    % paragraph: consecutive text lines up to the next
                    % blank line or block element, with inline markdown
                    para = {strtrim(ln)};
                    i = i + 1;
                    while i <= n && ~isempty(strtrim(lines{i})) && ...
                            isempty(regexp(lines{i}, blockStart, 'once')) && ...
                            isempty(regexp(lines{i}, marker, 'once'))
                        para{end+1} = strtrim(lines{i}); %#ok<AGROW>
                        i = i + 1;
                    end
                    out{end+1} = ['<p>' obj.relink(T.escape(strjoin(para, ' ')), ...
                        pageDir) '</p>']; %#ok<AGROW>
                    continue
                end
                i = i + 1;
            end
            closeLists();
            html = strjoin(out, newline);

            function closeLists()
                if inUl
                    out{end+1} = '</ul>';
                    inUl = false;
                end
                if inOl
                    out{end+1} = '</ol>';
                    inOl = false;
                end
            end
        end

    end

    methods (Access = private)

        function out = markerHtml(obj, m, pageDir, folder)
            % HTML of a {{...}} marker line (empty for an unknown marker)
            T = docgen.Text;
            main = strtok(folder, '/');
            out = cell(1, 0);
            if strcmp(m, 'SUBFOLDERS')
                folders = unique(obj.Rows(strcmp(obj.Rows(:, 5), main), 4));
                folders = folders(~strcmp(folders, main));
                [~, order] = sort(lower(folders));
                folders = folders(order);
                tbl = {'<table><tr><th>Subfolder</th><th>Description</th></tr>'};
                for k = 1:numel(folders)
                    b = baseName(folders{k});
                    tbl{end+1} = ['<tr><td><a href="' b '_overview.html">' b ...
                        '</a></td><td>' T.escape(obj.folderDesc(folders{k})) ...
                        '</td></tr>']; %#ok<AGROW>
                end
                tbl{end+1} = '</table>';
                out = {strjoin(tbl, newline)};
            elseif strcmp(m, 'RELATED_EXAMPLES')
                names = obj.relatedExamples(folder);
                if isempty(names)
                    out = {'<p>No example uses the classes of this folder yet.</p>'};
                else
                    links = cellfun(@(r) obj.relink(['`' r '`'], pageDir), names, ...
                        'UniformOutput', false);
                    out = {['<p>' strjoin(links, ', ') '</p>']};
                end
            elseif strcmp(m, 'MEMBERS') || startsWith(m, 'MEMBERS:')
                if startsWith(m, 'MEMBERS:')
                    folder = [main '/' m(9:end)];
                end
                mem = obj.folderMembers(folder);
                if ~isempty(mem)
                    tbl = {'<table><tr><th>Name</th><th>Description</th></tr>'};
                    for k = 1:size(mem, 1)
                        tbl{end+1} = ['<tr><td>' obj.Context.Lookup.linkNames(mem{k, 1}, ...
                            pageDir, true) '</td><td>' ...
                            obj.relink(T.escape(mem{k, 3}), pageDir) ...
                            '</td></tr>']; %#ok<AGROW>
                    end
                    tbl{end+1} = '</table>';
                    out = {strjoin(tbl, newline)};
                end
            end
        end

        function out = tableHtml(obj, rows, pageDir)
            % HTML of a Markdown table given as rows of stripped cells
            T = docgen.Text;
            out = {'<table>'};
            for r = 1:numel(rows)
                cells = rows{r};
                if all(cellfun(@(c) ~isempty(regexp(c, '^[-: ]+$', 'once')), cells))
                    continue  % separator row
                end
                if r == 1
                    out{end+1} = ['<tr>' strjoin(strcat('<th>', cells, '</th>'), '') ...
                        '</tr>']; %#ok<AGROW>
                elseif numel(cells) >= 2
                    % split multi-name cells ('BDMk / BDMk3d', 'A, B'),
                    % link each known name, and drop the name(s) from the
                    % description when it repeats them
                    if isempty(cells{1})
                        cell0 = '-';
                    else
                        cell0 = obj.Context.Lookup.linkNames(cells{1}, pageDir, true);
                    end
                    rest = T.escape(cells{2});
                    names = regexp(strtrim(cells{1}), '\s*[,/]\s*', 'split');
                    for k = 1:numel(names)
                        cleaned = docgen.stripName(names{k}, rest);
                        if ~strcmp(cleaned, rest)
                            rest = cleaned;
                            break
                        end
                    end
                    out{end+1} = ['<tr><td>' cell0 '</td><td>' ...
                        obj.relink(rest, pageDir) '</td></tr>']; %#ok<AGROW>
                else
                    out{end+1} = ['<tr><td>' obj.relink(T.escape(strjoin(cells, ' ')), ...
                        pageDir) '</td></tr>']; %#ok<AGROW>
                end
            end
            out{end+1} = '</table>';
        end

        function text = relink(obj, text, pageDir)
            % inline markdown of HTML-escaped text: `code`, [text](url),
            % **bold**, *italic*; known code and bold names are linked
            text = replaceMatches(text, '`([^`]+)`', ...
                @(t) obj.linked(t{1}, 'code', pageDir));
            text = replaceMatches(text, '\[([^\]]+)\]\(([^)]+)\)', ...
                @(t) ['<a href="' t{2} '">' t{1} '</a>']);
            text = replaceMatches(text, '\*\*([^*]+)\*\*', ...
                @(t) obj.linked(t{1}, 'strong', pageDir));
            text = replaceMatches(text, '(?<!\S)\*([^*<]+)\*(?!\S)', ...
                @(t) ['<em>' t{1} '</em>']);
        end

        function s = linked(obj, name, tag, pageDir)
            % <tag>name</tag>, linked to the page of name if it is known
            s = ['<' tag '>' name '</' tag '>'];
            href = obj.Context.Lookup.pathFor(name, pageDir);
            if ~isempty(href)
                s = ['<a href="' href '">' s '</a>'];
            end
        end

        function desc = folderDesc(obj, folder)
            % first non-empty comment line of the folder's Contents.m,
            % without '%' and without a leading folder name
            desc = '';
            file = fullfile(obj.Context.Toolbox, folder, 'Contents.m');
            if ~isfile(file)
                return
            end
            lines = docgen.Text.readLines(file);
            for k = 1:numel(lines)
                if ~startsWith(lines{k}, '%')
                    continue
                end
                t = strtrim(docgen.Help.uncomment(lines{k}));
                if ~isempty(t)
                    desc = docgen.stripName(baseName(folder), t);
                    return
                end
            end
        end

        function mem = folderMembers(obj, folder)
            % rows name, kind, description of the files of a source
            % folder, in docgen.memberOrder
            if isKey(obj.Members, folder)
                mem = obj.Members{folder};
                return
            end
            main = strtok(folder, '/');
            sel = find(strcmp(obj.Rows(:, 5), main) & strcmp(obj.Rows(:, 4), folder));
            mem = cell(numel(sel), 3);
            for k = 1:numel(sel)
                src = obj.Rows{sel(k), 2};
                stem = docgen.Text.stem(src);
                [~, h1] = docgen.Help.fileSummary( ...
                    docgen.Text.readLines(fullfile(obj.Context.Toolbox, src)));
                mem(k, :) = {stem, obj.Rows{sel(k), 3}, docgen.stripName(stem, h1)};
            end
            mem = mem(docgen.memberOrder(mem(:, 1), mem(:, 2)), :);
            obj.Members{folder} = mem;
        end

        function names = relatedExamples(obj, folder)
            % examples (with a page) that use a member of the folder or of
            % one of its subfolders
            main = strtok(folder, '/');
            folders = unique(obj.Rows(strcmp(obj.Rows(:, 5), main), 4));
            folders = folders(strcmp(folders, folder) | startsWith(folders, [folder '/']));
            members = cell(0, 1);
            for k = 1:numel(folders)
                mem = obj.folderMembers(folders{k});
                members = [members; mem(:, 1)]; %#ok<AGROW>
            end
            ex = obj.Context.Examples;
            uses = arrayfun(@(e) any(ismember(e.refs, members)), ex);
            names = intersect({ex(uses).name}, obj.Context.ExamplePages);  % sorted
        end

        function file = write(obj, rel, txt)
            % write a page below the output folder
            file = fullfile(obj.Context.DocFolder, rel);
            folder = fileparts(file);
            if ~isfolder(folder)
                mkdir(folder);
            end
            docgen.Text.writeText(file, txt);
            fprintf('wrote %s\n', file);
        end

    end
end

function b = baseName(folder)
% last component of a '/'-separated folder
b = folder(find(folder == '/', 1, 'last')+1:end);
end

function out = replaceMatches(text, pat, fn)
% replace each match of PAT by FN(TOKENS), FN a function of the cell
% array of tokens of the match
[toks, s, e] = regexp(text, pat, 'tokens', 'start', 'end');
out = '';
last = 0;
for k = 1:numel(s)
    out = [out text(last+1:s(k)-1) fn(toks{k})]; %#ok<AGROW>
    last = e(k);
end
out = [out text(last+1:end)];
end