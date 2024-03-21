function groups = genExamplesOverview(options)
%genExamplesOverview - Generate the examples overview of the documentation
%   genExamplesOverview() writes toolbox/doc/examples/examples_overview.html:
%   the examples grouped by topic, with their titles and the toolbox
%   classes and functions they use (docgen.exampleInfo). Heading, intro
%   text and topics come from docbuild/overviews/examples.md ('## Topic'
%   followed by '- exampleName' lines); examples not listed there appear
%   under 'Other'. Only examples with a page (see exportExamples) are
%   linked.
%
%   genExamplesOverview(DocFolder=DIR) writes the page to DIR/examples
%   instead of toolbox/doc/examples.
%
%   genExamplesOverview(...,Context=CTX) uses the docgen.Context CTX (and
%   its DocFolder) instead of scanning the toolbox again.
%
%   GROUPS = genExamplesOverview(...) returns the topics and their
%   examples as an N-by-2 cell array (see docgen.exampleTopics).
%
%   See also exportExamples, genOverviews, docgen.Context, buildtool

% Copyright (c) 2026 by the fem4m authors

arguments
    options.DocFolder (1,1) string = ""
    options.Context docgen.Context {mustBeScalarOrEmpty} = docgen.Context.empty
end
T = docgen.Text;
ctx = docgen.Context.resolve(options.Context, options.DocFolder);

% ---- examples grouped by topic (docbuild/overviews/examples.md) ----
topics = docgen.exampleTopics(ctx);
if ~isempty(topics.unknown)
    fprintf('examples.md: no example %s\n', strjoin(topics.unknown, ', '));
end
if ~isempty(topics.untopiced)
    fprintf('examples without a topic in examples.md (listed under Other): %s\n', ...
        strjoin(topics.untopiced, ', '));
end
groups = topics.groups;
heading = topics.heading;
intro = topics.intro;

% ---- examples_overview.html ----
% the intro is Markdown like the overview drafts
intro = docgen.OverviewGenerator(ctx).mdToHtml(intro, 'examples', '');
out = {['<h1>' T.escape(heading) '</h1>']};
if ~isempty(intro)
    out{end+1} = intro;
end
for g = 1:size(groups, 1)
    out{end+1} = ['<h2>' T.escape(groups{g, 1}) '</h2><table><tr><th>Example</th>' ...
        '<th>Title</th><th>Uses</th></tr>']; %#ok<AGROW>
    for name = groups{g, 2}
        e = ctx.Examples(strcmp({ctx.Examples.name}, name{1}));
        if ctx.hasExamplePage(e.name)
            link = ['<a href="' e.name '.html"><code>' e.name '</code></a>'];
        else
            fprintf('no page for %s\n', e.name);
            link = ['<code>' e.name '</code>'];
        end
        uses = ctx.Lookup.linkNames(strjoin(e.refs, ', '), 'examples', true);
        out{end+1} = ['<tr><td>' link '</td><td>' docgen.exampleTitle(e.title, e.name) ...
            '</td><td>' uses '</td></tr>']; %#ok<AGROW>
    end
    out{end+1} = '</table>'; %#ok<AGROW>
end
file = fullfile(ctx.DocFolder, 'examples', 'examples_overview.html');
if ~isfolder(fileparts(file))
    mkdir(fileparts(file));
end
T.writeText(file, docgen.page(ctx, heading, {heading, ''}, out, Root='../'));
fprintf('wrote %s\n', file);
fprintf('examples sections: %s\n', strjoin(groups(:, 1), ', '));
end