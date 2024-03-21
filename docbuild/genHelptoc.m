function file = genHelptoc(options)
%genHelptoc - Generate the table of contents of the fem4m documentation
%   genHelptoc() writes toolbox/doc/helptoc.xml, the table of contents of
%   the MATLAB Help browser, from the page tree (docgen.pageTree): the start
%   page, Getting Started, one entry per main folder with its classes and
%   functions (classes first, each alphabetical) and its subfolders, and
%   the examples that have a page (docgen.Context/ExamplePages).
%
%   genHelptoc(DocFolder=DIR) uses the documentation folder DIR instead of
%   toolbox/doc.
%
%   genHelptoc(...,Context=CTX) uses the docgen.Context CTX (and its
%   DocFolder) instead of scanning the toolbox again.
%
%   FILE = genHelptoc(...) returns the file written. checkDocs checks
%   that every target exists, once, and that every page is in the TOC.
%
%   See also genRef, genOverviews, genExamplesOverview, checkDocs, buildtool

% Copyright (c) 2026 by the fem4m authors

arguments
    options.DocFolder (1,1) string = ""
    options.Context docgen.Context {mustBeScalarOrEmpty} = docgen.Context.empty
end
ctx = docgen.Context.resolve(options.Context, options.DocFolder);
file = fullfile(ctx.DocFolder, 'helptoc.xml');

% reference pages as rows: name, page, kind, folder
rows = [cellfun(@docgen.Text.stem, ctx.Pages(:, 2), 'UniformOutput', false), ...
    ctx.Pages(:, [1 3 4])];
folders = rows(:, 4);

lines = {'<?xml version="1.0" encoding="utf-8"?>'
    '<toc version="2.0">'
    ['  <tocitem target="index.html">' docgen.Text.escape(ctx.Site)]
    ['    <tocitem target="getting_started.html" ' ...
    'image="HelpIcon.GETTING_STARTED">Getting Started</tocitem>']};

mains = unique(strtok(folders, '/'));
for k = 1:numel(mains)
    m = mains{k};
    lines{end+1, 1} = ['    <tocitem target="' m '/' m '_overview.html" ' ...
        'image="HelpIcon.FUNCTION">' docgen.Text.escape(m)]; %#ok<AGROW>
    % files of the main folder itself first, then the subfolders with
    % their files nested in the subfolder entry
    lines = [lines; memberItems(rows, m, m, '      ')]; %#ok<AGROW>
    subs = unique(folders(startsWith(folders, [m '/'])));
    for s = 1:numel(subs)
        [~, sub] = fileparts(subs{s});  % last folder name
        lines{end+1, 1} = ['      <tocitem target="' m '/' sub '_overview.html">' ...
            sub]; %#ok<AGROW>
        lines = [lines; memberItems(rows, subs{s}, m, '        ')]; %#ok<AGROW>
        lines{end+1, 1} = '      </tocitem>'; %#ok<AGROW>
    end
    lines{end+1, 1} = '    </tocitem>'; %#ok<AGROW>
end

% examples with a page, sorted by name
lines{end+1, 1} = ['    <tocitem target="examples/examples_overview.html" ' ...
    'image="HelpIcon.EXAMPLES">Examples'];
for name = ctx.ExamplePages
    lines{end+1, 1} = ['      <tocitem target="examples/' name{1} '.html">' ...
        name{1} '</tocitem>']; %#ok<AGROW>
end
lines = [lines; {'    </tocitem>'; '  </tocitem>'; '</toc>'}];

docgen.Text.writeText(file, [strjoin(lines, newline) newline]);
fprintf('wrote %s, %d lines\n', file, numel(lines));
end

function items = memberItems(rows, folder, main, pad)
% TOC entries of the pages in a source folder (docgen.memberOrder)
sel = rows(strcmp(rows(:, 4), folder), :);
sel = sel(docgen.memberOrder(sel(:, 1), sel(:, 3)), :);
items = cell(size(sel, 1), 1);
for k = 1:size(sel, 1)
    items{k} = [pad '<tocitem target="' main '/' sel{k, 2} '">' sel{k, 1} ...
        '</tocitem>'];
end
end