function gen = genRef(folders, options)
%genRef - Generate the HTML reference pages of the fem4m toolbox
%   genRef() writes the reference page of every class and function of the
%   page tree (docgen.pageTree) to toolbox/doc/<main folder> and prints
%   what may need a look in the help texts (members without a description,
%   ambiguous names, Syntax notes). checkDocs checks the links.
%
%   genRef(FOLDERS) only writes the pages of the source folders in the
%   cell array FOLDERS, e.g. genRef({'forms','elements/mixed'}).
%
%   genRef(...,DocFolder=DIR) uses the documentation folder DIR instead
%   of toolbox/doc: pages go to DIR/<main folder>, and the example and
%   overview pages they refer to are looked up in DIR.
%
%   genRef(...,Context=CTX) uses the docgen.Context CTX (and its
%   DocFolder) instead of scanning the toolbox again.
%
%   GEN = genRef(...) also returns the docgen.RefGenerator.
%
%   See also docgen.RefGenerator, docgen.Context, buildtool

% Copyright (c) 2026 by the fem4m authors

arguments
    folders cell = {}
    options.DocFolder (1,1) string = ""
    options.Context docgen.Context {mustBeScalarOrEmpty} = docgen.Context.empty
end
ctx = docgen.Context.resolve(options.Context, options.DocFolder);

% the class members are read with meta.class: the toolbox must be on the path
oldPath = addpath(genpath(ctx.Toolbox));
restorePath = onCleanup(@() path(oldPath));

gen = docgen.RefGenerator(ctx);
gen.run(folders);
gen.printNotes();
end