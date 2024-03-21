function result = buildDocs(options)
%buildDocs - Generate and check the documentation of the fem4m toolbox
%   RESULT = buildDocs() writes the documentation to toolbox/doc and checks
%   it: stylesheet, overview pages (genOverviews), reference pages
%   (genRef), examples overview (genExamplesOverview), table of contents
%   (genHelptoc), then checkDocs, whose result is returned. The example
%   pages are made by exportExamples; the other pages link only to
%   example pages that exist. The pages of an earlier build are deleted
%   first, except the example pages of examples that still exist.
%
%   RESULT = buildDocs(DocFolder=DIR) uses the documentation folder DIR
%   instead, e.g. a temporary folder for a test.
%
%   See also checkDocs, exportExamples, buildtool

% Copyright (c) 2026 by the fem4m authors

arguments
    options.DocFolder (1,1) string = ""
end
ctx = docgen.Context(options.DocFolder);  % the toolbox is scanned once
removeStalePages(ctx);
cssFolder = fullfile(ctx.DocFolder, 'css');
if ~isfolder(cssFolder)
    mkdir(cssFolder);
end
copyfile(fullfile(ctx.Build, 'doc.css'), cssFolder);
genOverviews(Context=ctx);
genRef(Context=ctx);
genExamplesOverview(Context=ctx);
genHelptoc(Context=ctx);
result = checkDocs(Context=ctx);
end

function removeStalePages(ctx)
% Delete the pages of an earlier build, so pages of renamed or deleted
% files do not survive: all pages outside doc/examples (they are written
% again) and the example pages whose example no longer exists. Only a
% folder with helptoc.xml, i.e. an earlier build, is cleaned.
doc = ctx.DocFolder;
if ~isfile(fullfile(doc, 'helptoc.xml'))
    return
end
doc = dir(doc);
doc = doc(1).folder;  % absolute, as dir writes it in the folder field
examples = fullfile(doc, 'examples');
old = [dir(fullfile(doc, '*.html')); dir(fullfile(doc, '*', '*.html'))];
old = old(~strcmp({old.folder}, examples));
keep = [strcat({ctx.Examples.name}, '.html') {'examples_overview.html'}];
ex = dir(fullfile(examples, '*.html'));
old = [old; ex(~ismember({ex.name}, keep))];
for k = 1:numel(old)
    delete(fullfile(old(k).folder, old(k).name));
end
% remove the folders left empty, e.g. of a main folder that is gone
for f = setdiff(unique({old.folder}), {doc})
    left = dir(f{1});
    if all(ismember({left.name}, {'.', '..', '.DS_Store'}))
        rmdir(f{1}, 's');
    end
end
end