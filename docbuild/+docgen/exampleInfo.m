function examples = exampleInfo(toolbox)
%exampleInfo - Name, title and used toolbox names of every example
%   EXAMPLES = docgen.exampleInfo(TOOLBOX) scans the examples
%   TOOLBOX/examples/*.m and returns a struct array, sorted by name, with
%   the fields
%       name  - name of the example, e.g. 'exPoissonP1'
%       file  - path relative to TOOLBOX, e.g. 'examples/exPoissonP1.m'
%       title - first Markdown heading in a comment ('# Title'), or ''
%               (see docgen.exampleTitle)
%       refs  - names of the toolbox .m files outside the examples that
%               occur in the code of the example (as whole words), sorted;
%               commented-out code ('% mark = MarkUniform;') counts,
%               comment text ('% FE-space') does not
%
%   See also docgen.pageTree, docgen.exampleTitle, genExamplesOverview

% Copyright (c) 2026 by the fem4m authors

arguments
    toolbox (1,:) char
end

T = docgen.Text;

% names of the toolbox .m files outside the examples
files = dir(fullfile(toolbox, '**', '*.m'));
exampleFolder = fullfile(toolbox, 'examples');
isExample = strcmp({files.folder}, exampleFolder);
symbols = unique(cellfun(@T.stem, {files(~isExample).name}, 'UniformOutput', false));

names = sort(cellfun(@T.stem, {files(isExample).name}, 'UniformOutput', false));
examples = struct('name', names, 'file', strcat('examples/', names, '.m'), ...
    'title', '', 'refs', {cell(1, 0)})';
for k = 1:numel(examples)
    lines = T.readLines(fullfile(toolbox, examples(k).file));
    for i = 1:numel(lines)
        tok = regexp(strtrim(docgen.Help.uncomment(lines{i})), '#\s*(.+)', ...
            'tokens', 'once');
        if startsWith(strtrim(lines{i}), '%') && ~isempty(tok)
            examples(k).title = strtrim(tok{1});
            break
        end
    end
    code = exampleCode(lines);
    used = cellfun(@(s) ~isempty(regexp(code, ['\<' ...
        regexptranslate('escape', s) '\>'], 'once')), symbols);
    examples(k).refs = symbols(used);
end
end

function code = exampleCode(lines)
% the code of an example as one text: the code lines without their
% comments and the comment lines that are commented-out code, i.e. an
% assignment ('% p = problemX();') or a call ('% plot(mesh)'); rich text
% lines ('%[text] ...') are never code
codeRe = ['^\s*(\[[^\]]*\]|[A-Za-z_][\w.]*(\{[^}]*\}|\([^)]*\))?)\s*=[^=]' ...
    '|^\s*[A-Za-z_][\w.]*\s*\(.*\)\s*;?\s*$'];
code = repmat({''}, size(lines));
for i = 1:numel(lines)
    t = strtrim(lines{i});
    if ~startsWith(t, '%')
        code{i} = docgen.Help.stripComment(lines{i});
    elseif ~startsWith(t, '%[')
        c = strtrim(docgen.Help.uncomment(t));
        if ~isempty(regexp(c, codeRe, 'once'))
            code{i} = docgen.Help.stripComment(c);
        end
    end
end
code = strjoin(code, newline);
end