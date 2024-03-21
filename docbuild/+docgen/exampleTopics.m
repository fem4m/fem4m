function topics = exampleTopics(ctx)
%exampleTopics - Examples grouped by the topics of overviews/examples.md
%   TOPICS = docgen.exampleTopics(CTX) reads docbuild/overviews/examples.md
%   of the docgen.Context CTX: a '# Heading', an intro in Markdown, then
%   '## Topic' lines, each followed by '- exampleName' lines. TOPICS is a
%   struct with the fields
%       heading   - the heading ('Examples' if there is none)
%       intro     - the Markdown lines between the heading and the first
%                   topic (see docgen.OverviewGenerator/mdToHtml)
%       groups    - N-by-2 cell array {topic, example names}; examples not
%                   listed in examples.md come last, under 'Other'
%       unknown   - names listed in examples.md that are no example
%       untopiced - examples not listed in examples.md (sorted)
%
%   See also genExamplesOverview, docgen.exampleInfo

% Copyright (c) 2026 by the fem4m authors

examples = {ctx.Examples.name};
[topics.heading, topics.intro, listed] = readDraft(fullfile(ctx.Build, ...
    'overviews', 'examples.md'));
topics.groups = cell(0, 2);
topics.unknown = cell(1, 0);
seen = cell(1, 0);
for k = 1:size(listed, 1)
    names = listed{k, 2};
    topics.unknown = [topics.unknown names(~ismember(names, examples))];
    names = names(ismember(names, examples));
    seen = [seen names]; %#ok<AGROW>
    if ~isempty(names)
        topics.groups(end+1, :) = {listed{k, 1}, names};
    end
end
rest = examples(~ismember(examples, seen));
[~, order] = sort(lower(rest));
topics.untopiced = reshape(rest(order), 1, []);
if ~isempty(rest)
    topics.groups(end+1, :) = {'Other', topics.untopiced};
end
end

function [heading, intro, topics] = readDraft(file)
% heading, intro lines and topics ({name, examples}) of the draft
heading = 'Examples';
intro = cell(1, 0);
topics = cell(0, 2);
for ln = docgen.Text.readLines(file)
    t = strtrim(ln{1});
    if startsWith(t, '## ')
        topics(end+1, :) = {strtrim(t(4:end)), cell(1, 0)}; %#ok<AGROW>
    elseif startsWith(t, '# ') && isempty(topics)
        heading = strtrim(t(3:end));
    elseif startsWith(t, '- ') && ~isempty(topics)
        topics{end, 2}{end+1} = strrep(strtrim(t(3:end)), '`', '');
    elseif isempty(topics)
        intro{end+1} = ln{1}; %#ok<AGROW>
    end
end
end