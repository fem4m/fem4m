function result = checkDocs(options)
%checkDocs - Check the links and the completeness of the documentation
%   checkDocs() checks all HTML pages of toolbox/doc for
%   - links to pages or files that do not exist, or to a #fragment that
%     is no id or name on the page (broken internal links),
%   - images that do not exist,
%   - pages without an entry in helptoc.xml,
%   - TOC entries whose target does not exist or that repeat a target,
%   - reference pages (docgen.pageTree) whose file is missing,
%   and prints the counts and every problem found.
%
%   checkDocs(DocFolder=DIR) checks the documentation folder DIR instead.
%
%   checkDocs(...,Context=CTX) uses the docgen.Context CTX (and its
%   DocFolder) instead of scanning the toolbox again.
%
%   RESULT = checkDocs(...) returns a struct with the fields broken and
%   missingImages (N-by-2: page, link), notInToc, tocMissing and
%   tocDuplicates (TOC targets) and unreferenced (N-by-2: source, page).
%
%   See also genRef, genOverviews, genHelptoc, buildtool

% Copyright (c) 2026 by the fem4m authors

arguments
    options.DocFolder (1,1) string = ""
    options.Context docgen.Context {mustBeScalarOrEmpty} = docgen.Context.empty
end
ctx = docgen.Context.resolve(options.Context, options.DocFolder);
doc = ctx.DocFolder;

% TOC targets
toc = regexp(fileread(fullfile(doc, 'helptoc.xml')), ...
    'target="([^"#]+?)(?:#[^"]*)?"', 'tokens');
toc = cellfun(@(t) t{1}, toc, 'UniformOutput', false);
toc = toc(~startsWith(toc, 'http'));
[tocTargets, ~, idx] = unique(toc);
tocDuplicates = tocTargets(accumarray(idx(:), 1) > 1);

% all files and folders below doc, relative to doc with '/' separators:
% links are resolved on these paths, so the check does not depend on the
% file separator of the platform
docAbs = dir(doc);
docAbs = docAbs(1).folder;  % doc as dir writes it in the folder field
entries = dir(fullfile(docAbs, '**'));
isFile = ~[entries.isdir];
files = sort(relPath(strcat({entries(isFile).folder}, filesep, ...
    {entries(isFile).name}), docAbs));
folders = setdiff(relPath(unique({entries(~isFile).folder}), docAbs), {''});
known = [files folders];
pages = files(endsWith(files, '.html'));

% reference pages whose file does not exist
planned = strcat(ctx.Pages(:, 5), '/', ctx.Pages(:, 1));
missing = ~ismember(planned, pages);
unreferenced = [ctx.Pages(missing, 2), planned(missing)];

% internal links (with their #fragment) and images of all pages
texts = cellfun(@(p) fileread(fullfile(docAbs, p), 'Encoding', 'UTF-8'), ...
    pages, 'UniformOutput', false);
anchors = dictionary(string(pages), cellfun(@anchorIds, texts, ...
    'UniformOutput', false));  % page -> ids
broken = cell(0, 2);
missingImages = cell(0, 2);
for k = 1:numel(pages)
    t = texts{k};
    base = fileparts(pages{k});  % folder of the page, '' for the doc root
    hrefs = regexp(t, 'href="((?!http|mailto|matlab:)[^"]*)"', 'tokens');
    for h = 1:numel(hrefs)
        [link, frag] = splitFragment(hrefs{h}{1});
        target = pages{k};  % '#frag' refers to the page itself
        if ~isempty(link)
            [ok, target] = exists(link, base, known, docAbs);
            if ~ok
                broken(end+1, :) = {pages{k}, link}; %#ok<AGROW>
                continue
            end
        end
        if ~isempty(frag) && endsWith(target, '.html') && ...
                ~ismember(frag, anchors{target})
            broken(end+1, :) = {pages{k}, hrefs{h}{1}}; %#ok<AGROW>
        end
    end
    srcs = regexp(t, 'src="((?!http|data:)[^"]+)"', 'tokens');
    for s = 1:numel(srcs)
        if ~exists(srcs{s}{1}, base, known, docAbs)
            missingImages(end+1, :) = {pages{k}, srcs{s}{1}}; %#ok<AGROW>
        end
    end
end

notInToc = pages(~ismember(pages, tocTargets));
tocMissing = tocTargets(~ismember(tocTargets, known));

% ---- report ----
fprintf('HTML pages: %d, TOC targets: %d\n', numel(pages), numel(tocTargets));
pairs = @(a, sep, b) cellfun(@(x, y) [x sep y], a, b, 'UniformOutput', false);
report('broken internal links', pairs(broken(:, 1), ' -> ', broken(:, 2)));
report('missing images', pairs(missingImages(:, 1), ' -> ', missingImages(:, 2)));
report('pages not in TOC', notInToc);
report('TOC targets without file', tocMissing);
report('TOC targets listed twice', tocDuplicates);
report('planned pages without file', pairs(unreferenced(:, 2), ' from ', ...
    unreferenced(:, 1)));

result = struct('broken', {broken}, 'missingImages', {missingImages}, ...
    'notInToc', {notInToc}, 'tocMissing', {tocMissing}, ...
    'tocDuplicates', {tocDuplicates}, 'unreferenced', {unreferenced});
end

function p = relPath(files, folder)
% FILES (cell array) relative to FOLDER, with '/' separators; FOLDER
% itself becomes ''
p = cellfun(@(f) strrep(f(min(numel(folder)+2, numel(f)+1):end), ...
    filesep, '/'), files, 'UniformOutput', false);
end

function [tf, target] = exists(link, base, known, docAbs)
% true if the relative LINK of a page in folder BASE (relative to the doc
% folder) points to a file or folder in KNOWN; a link that leaves the doc
% folder, or starts with '/', is checked on disk. TARGET is the path of
% the link relative to the doc folder ('' if it leaves the doc folder).
target = '';
if startsWith(link, '/')
    tf = isfile(link) || isfolder(link);
    return
end
rel = normPath([base '/' link]);
if startsWith(rel, '..')
    p = fullfile(docAbs, strrep(rel, '/', filesep));
    tf = isfile(p) || isfolder(p);
else
    tf = ismember(rel, known);
    target = rel;
end
end

function [link, frag] = splitFragment(href)
% 'page.html#id' -> 'page.html', 'id'; '#id' -> '', 'id'
i = find(href == '#', 1);
if isempty(i)
    link = href;
    frag = '';
else
    link = href(1:i-1);
    frag = href(i+1:end);
end
end

function ids = anchorIds(html)
% targets of '#' links in a page: id and name attributes
ids = regexp(html, '\s(?:id|name)="([^"]*)"', 'tokens');
ids = unique(cellfun(@(c) c{1}, ids, 'UniformOutput', false));
end

function p = normPath(p)
% collapse '.', '..' and repeated '/' as a browser resolves a relative
% link (lexically, also through folders that do not exist); leading '..'
% are kept for a link that leaves the doc folder
parts = strsplit(p, '/');
keep = cell(1, 0);
for k = 1:numel(parts)
    if isempty(parts{k}) || strcmp(parts{k}, '.')
        continue
    elseif strcmp(parts{k}, '..') && ~isempty(keep) && ~strcmp(keep{end}, '..')
        keep(end) = [];
    else
        keep{end+1} = parts{k}; %#ok<AGROW>
    end
end
p = strjoin(keep, '/');
end

function report(what, items)
% print the number of problems of one kind and each problem
fprintf('%s: %d\n', what, numel(items));
for k = 1:numel(items)
    fprintf('  %s\n', items{k});
end
end