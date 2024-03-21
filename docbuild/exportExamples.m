function out = exportExamples(names, options)
%exportExamples - Run the examples and export them to HTML
%   OUT = exportExamples() runs every example of toolbox/examples and
%   exports it with its output to toolbox/doc/examples/<name>.html. Each
%   page gets a readable <title> (the '# Title' heading of the example,
%   see docgen.exampleInfo) and a navigation bar linking back to the
%   documentation. Inline code in the text of an example that is exactly
%   the name of a class, a function or another example with a page is
%   linked to that page (docgen.Lookup). The figures are exported in the
%   light graphics theme, also from a desktop in the dark theme. OUT is a
%   struct array with fields name, ok and msg, one element per example.
%
%   OUT = exportExamples(NAMES) runs only the examples in the cell array
%   NAMES, e.g. exportExamples({'exPoissonP1'}).
%
%   OUT = exportExamples(...,DocFolder=DIR) writes the pages to
%   DIR/examples instead of toolbox/doc/examples.
%
%   The examples are run by export(...,Run=true), which runs each one in
%   a cleared base workspace and restores the workspace afterwards; the
%   example sources are not changed. All open figures are closed before
%   each example.
%
%   See also genExamplesOverview, docgen.exampleInfo, buildtool

% Copyright (c) 2026 by the fem4m authors

arguments
    names cell = {}
    options.DocFolder (1,1) string = ""
end
ctx = docgen.Context(options.DocFolder);
toolbox = ctx.Toolbox;
if isempty(names)  % every example in toolbox/examples
    names = {ctx.Examples.name};
end

oldPath = addpath(genpath(toolbox));
restorePath = onCleanup(@() path(oldPath));
outFolder = fullfile(ctx.DocFolder, 'examples');
if ~isfolder(outFolder)
    mkdir(outFolder);
end
% dark lines such as mesh edges vanish on the dark background of the
% dark theme; the temporary value is not stored in the settings
s = settings;
graphicsTheme = s.matlab.appearance.figure.GraphicsTheme;
oldTheme = [];
if hasTemporaryValue(graphicsTheme)
    oldTheme = graphicsTheme.TemporaryValue;
end
graphicsTheme.TemporaryValue = "light";
restoreTheme = onCleanup(@() resetTheme(graphicsTheme, oldTheme));
% the examples that have a page or get one now
lookup = docgen.Lookup(ctx.Pages, union(ctx.ExamplePages, names));
out = struct('name', names, 'ok', false, 'msg', '');
for k = 1:numel(names)
    src = fullfile(toolbox, 'examples', [names{k} '.m']);
    dst = fullfile(outFolder, [names{k} '.html']);
    [out(k).ok, out(k).msg] = exportOne(src, dst);
    if isfile(dst)
        e = ctx.Examples(strcmp({ctx.Examples.name}, names{k}));
        [~, title] = docgen.exampleTitle([e.title], names{k});
        finishPage(dst, title, ctx.Site, lookup, names{k});
    end
    status = 'FAILED';
    if out(k).ok
        status = 'ok';
    end
    fprintf('%s : %s : %s\n', names{k}, status, out(k).msg);
end
end

function [ok, msg] = exportOne(src, dst)
% run the example SRC and export it to DST; an error of the example is
% read from the exported page
% figure(k) would reuse a figure left open by the previous example, and
% that figure would be missing on the page
close all force
try
    % CatchError=false cannot be used: in R2026a export then reports an
    % error for every run, so errors are read from the page
    export(src, dst, Run=true);
catch e
    ok = false;
    msg = ['export: ' e.message];
    return
end
ok = true;
msg = '';
t = fileread(dst, 'Encoding', 'UTF-8');
t = t(strfind(t, '<body'):end);  % the stylesheet names the class, too
s = regexp(t, '<div[^>]*class="[^"]*embeddedOutputsErrorElement', 'once');
if ~isempty(s)
    ok = false;
    % the error message is the first text in the error element
    txt = regexp(t(s:end), '>\s*([^<]*[^<\s])', 'tokens', 'once');
    msg = 'runtime error (see the exported page)';
    if ~isempty(txt)
        msg = ['runtime error: ' txt{1}(1:min(end, 200))];
    end
end
end

function finishPage(file, title, site, lookup, name)
% give the exported page FILE of the example NAME the <title> TITLE
% (plain text, see docgen.exampleTitle), a navigation bar linking back
% to the docs of the toolbox SITE and links for the names in its inline
% code (docgen.Lookup LOOKUP); the bar has its own style, since the page
% does not load doc.css
T = docgen.Text;
navMark = '<!-- fem4m-nav -->';
nav = [navMark '<div style="font: 0.85rem -apple-system, ''Segoe UI'', ' ...
    'Helvetica, Arial, sans-serif; margin: 0.5rem 0 1rem;">' ...
    '<a href="../index.html">' T.escape(site) '</a> &#8250; ' ...
    '<a href="examples_overview.html">Examples</a></div>'];
t = regexprep(fileread(file, 'Encoding', 'UTF-8'), '\r\n?', '\n');
titleTag = ['<title>' T.escape(title) '</title>'];
[s, e] = regexp(t, '<title>.*?</title>', 'once');
if isempty(s)
    t2 = t;
else
    t2 = [t(1:s-1) titleTag t(e+1:end)];
end
if ~contains(t2, '<title>')
    s = strfind(t2, '<head>');
    if ~isempty(s)
        t2 = [t2(1:s(1)+5) titleTag t2(s(1)+6:end)];
    end
end
if ~contains(t2, navMark)  % the mark also shows that the names are linked
    t2 = linkCode(t2, lookup, name);
    % a white page as in doc.css: the page has no background of its own,
    % and a browser in dark mode would show its dark text on a dark one
    t2 = regexprep(t2, '</head>', '<style>html { background: #ffffff; }</style></head>', 'once');
    e = regexp(t2, '<body[^>]*>', 'end', 'once');
    if ~isempty(e)
        t2 = [t2(1:e) nav t2(e+1:end)];
    end
end
if ~strcmp(t2, t)
    T.writeText(file, t2);
end
end

function t = linkCode(t, lookup, name)
% link the inline code of the text of the example NAME that is exactly a
% name known to LOOKUP to its page; the code of the example and its
% outputs are not inline code
pat = '<span class="_richTextNode" style='' font-family: monospace;''>([A-Za-z][\w.]*)</span>';
[spans, tok, parts] = regexp(t, pat, 'match', 'tokens', 'split');
out = parts(1);
for k = 1:numel(spans)
    href = lookup.pathFor(tok{k}{1}, 'examples', ['examples/' name '.m']);
    if isempty(href)
        out{end+1} = spans{k}; %#ok<AGROW>
    else
        out{end+1} = ['<a href="' href '">' spans{k} '</a>']; %#ok<AGROW>
    end
    out{end+1} = parts{k+1}; %#ok<AGROW>
end
t = [out{:}];
end

function resetTheme(graphicsTheme, oldTheme)
% restore the temporary value OLDTHEME of the setting GRAPHICSTHEME ([]:
% none)
if isempty(oldTheme)
    clearTemporaryValue(graphicsTheme);
else
    graphicsTheme.TemporaryValue = oldTheme;
end
end