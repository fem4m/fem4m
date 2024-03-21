function pages = pageTree(toolbox)
%pageTree - Reference pages of the documentation, from the toolbox files
%   PAGES = docgen.pageTree(TOOLBOX) lists the documented files of the
%   toolbox folder TOOLBOX. Every class and function file gets a reference
%   page, except
%       Contents.m files,
%       files directly in TOOLBOX,
%       the folders doc, examples and apps (apps are not documented).
%   PAGES is an N-by-5 cell array, one row {page, src, kind, folder, main}
%   per page, e.g. {'Mesh.html', 'mesh/Mesh.m', 'class', 'mesh', 'mesh'},
%   sorted by src; kind is 'class' or 'function', main the first folder
%   of src. The examples are listed by docgen.exampleInfo.
%
%   See also docgen.exampleInfo, docgen.Lookup, docgen.RefGenerator, genHelptoc

% Copyright (c) 2026 by the fem4m authors

arguments
    toolbox (1,:) char
end

files = dir(fullfile(toolbox, '**', '*.m'));
src = strrep(strcat({files.folder}, filesep, {files.name}), [toolbox filesep], '');
src = sort(strrep(src, filesep, '/'));
skip = startsWith(src, {'doc/', 'examples/', 'apps/'}) | ...
    endsWith(src, '/Contents.m') | ~contains(src, '/');
src = src(~skip);

pages = cell(0, 5);
for k = 1:numel(src)
    kind = docgen.Help.fileSummary(docgen.Text.readLines(fullfile(toolbox, src{k})));
    if ~any(strcmp(kind, {'class', 'function'}))
        continue  % scripts and comment-only files get no page
    end
    folder = src{k}(1:find(src{k} == '/', 1, 'last')-1);
    pages(end+1, :) = {[docgen.Text.stem(src{k}) '.html'], src{k}, ...
        kind, folder, strtok(folder, '/')}; %#ok<AGROW>
end
end