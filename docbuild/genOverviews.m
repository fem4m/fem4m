function gen = genOverviews(mains, options)
%genOverviews - Generate the overview pages of the fem4m documentation
%   genOverviews() writes the overview page of every main folder and
%   subfolder of the page tree (docgen.pageTree) to
%   toolbox/doc/<main folder>, and toolbox/doc/index.html and
%   getting_started.html. The texts come from the Markdown drafts in
%   docbuild/overviews; the examples that use a folder's classes come
%   from docgen.exampleInfo.
%
%   genOverviews(MAINS) only writes the pages of the main folders in the
%   cell array MAINS, e.g. genOverviews({'mesh'}), and not the start pages.
%
%   genOverviews(...,DocFolder=DIR) uses the documentation folder DIR
%   instead of toolbox/doc (pages written there, example pages read).
%
%   genOverviews(...,Context=CTX) uses the docgen.Context CTX (and its
%   DocFolder) instead of scanning the toolbox again.
%
%   GEN = genOverviews(...) also returns the docgen.OverviewGenerator.
%
%   See also genRef, docgen.OverviewGenerator, buildtool

% Copyright (c) 2026 by the fem4m authors

arguments
    mains cell = {}
    options.DocFolder (1,1) string = ""
    options.Context docgen.Context {mustBeScalarOrEmpty} = docgen.Context.empty
end
gen = docgen.OverviewGenerator(docgen.Context.resolve(options.Context, ...
    options.DocFolder));
gen.run(mains);
end