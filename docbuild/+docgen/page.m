function html = page(ctx, title, crumbs, body, options)
%page - HTML page of the documentation: head, breadcrumb, body and footer
%   HTML = docgen.page(CTX,TITLE,CRUMBS,BODY) returns a page with the
%   <title> TITLE (plain text) and the HTML BODY (char, or a cell array of
%   lines joined with newlines). The breadcrumb starts with a link to the
%   start page, followed by CRUMBS, an N-by-2 cell array {label, href}
%   (plain text label, href '' for the page itself), e.g.
%   {'mesh', 'mesh_overview.html'; 'Mesh', ''}. The name and version of
%   the toolbox in the breadcrumb and footer come from the docgen.Context
%   CTX.
%
%   HTML = docgen.page(...,Root=R) is for a page in a subfolder: R is the
%   path from the page to the doc folder, e.g. '../' (default: '').
%   HTML = docgen.page(...,BodyClass=C) sets the class of <body>.
%
%   See also docgen.RefGenerator, docgen.OverviewGenerator, genExamplesOverview

% Copyright (c) 2026 by the fem4m authors

arguments
    ctx (1,1) docgen.Context
    title {mustBeTextScalar}
    crumbs cell
    body
    options.Root (1,:) char = ''
    options.BodyClass (1,:) char = ''
end
E = @docgen.Text.escape;
if iscell(body)
    body = strjoin(body, newline);
end
site = ctx.Site;
items = {['<a href="' options.Root 'index.html">' E(site) '</a>']};
for k = 1:size(crumbs, 1)
    [label, href] = crumbs{k, :};
    if isempty(href)
        items{end+1} = ['<span>' E(label) '</span>']; %#ok<AGROW>
    else
        items{end+1} = ['<a href="' href '">' E(label) '</a>']; %#ok<AGROW>
    end
end
footer = [site ' documentation'];
if ~isempty(ctx.Version)
    footer = [site ' ' ctx.Version ' documentation'];
end
bodyTag = '<body>';
if ~isempty(options.BodyClass)
    bodyTag = ['<body class="' options.BodyClass '">'];
end
html = ['<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">' ...
    '<meta name="viewport" content="width=device-width, initial-scale=1">' ...
    '<title>' E(char(title)) '</title>' ...
    '<link rel="stylesheet" href="' options.Root 'css/doc.css"></head>' ...
    bodyTag newline ...
    '<nav class="crumb">' strjoin(items, ' &#8250; ') '</nav>' newline ...
    body newline ...
    '<footer>' E(footer) '</footer></body></html>' newline];
end