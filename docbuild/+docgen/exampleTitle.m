function [html, text] = exampleTitle(title, name)
%exampleTitle - HTML and plain text of an example title
%   [HTML,TEXT] = docgen.exampleTitle(TITLE) renders the inline math of
%   an example title (the '# Title' heading of the example): inside $...$,
%   '_x' or '_{x}' becomes a subscript, '^x' or '^{x}' a superscript,
%   '\_' an underscore and '\times' (also with the doubled backslash of
%   the live script) the multiplication sign. HTML is the escaped
%   title, e.g. 'Adaptive $P\_1$ ...' -> 'Adaptive P<sub>1</sub> ...'.
%   TEXT is the
%   plain text for a <title>, with subscripts inline and superscripts
%   after '^': 'Adaptive P1 ...', '$C^1$ ...' -> 'C^1 ...'.
%
%   [HTML,TEXT] = docgen.exampleTitle(TITLE,NAME) uses the example name
%   NAME instead of an empty TITLE.
%
%   See also docgen.exampleInfo, docgen.RefGenerator, genExamplesOverview

% Copyright (c) 2026 by the fem4m authors

if isempty(title) && nargin > 1
    title = name;
end
html = docgen.Text.escape(char(title));
[math, rest] = regexp(html, '\$([^$]*)\$', 'tokens', 'split');
out = rest{1};
for k = 1:numel(math)
    s = strrep(math{k}{1}, '\_', '_');
    s = regexprep(s, '\\{1,2}times\s*', char(215));
    s = regexprep(s, '_\{?([A-Za-z0-9,]+)\}?', '<sub>$1</sub>');
    s = regexprep(s, '\^\{?([A-Za-z0-9,]+)\}?', '<sup>$1</sup>');
    out = [out s rest{k+1}]; %#ok<AGROW>
end
html = out;

text = regexprep(strrep(html, '<sup>', '^'), '<[^>]+>', '');
entities = {'&lt;', '<'; '&gt;', '>'; '&quot;', '"'; '&#x27;', ''''; '&amp;', '&'};
for k = 1:size(entities, 1)  % &amp; last: '&amp;lt;' stays '&lt;'
    text = strrep(text, entities{k, 1}, entities{k, 2});
end
end