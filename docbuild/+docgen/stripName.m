function text = stripName(name, text)
%stripName - Remove a leading repetition of a name from a description
%   TEXT = docgen.stripName(NAME,TEXT) removes NAME from the start of
%   TEXT when TEXT repeats it:
%   1. TEXT starts with NAME (ignoring case) and a separator '-', '–',
%      '—' or ':' -> the rest after the separator.
%   2. The first word of TEXT is NAME in capitals (e.g. 'PK' for Pk),
%      followed by a capitalised word (MATLAB H1 style 'PK Lagrange
%      element') -> the rest after the first word.
%   For a dotted NAME the last component is tried as well. TEXT is
%   returned unchanged if no rule matches or the rest would be empty.
%
%   See also docgen.RefGenerator

% Copyright (c) 2026 by the fem4m authors

if isempty(text) || isempty(name)
    return
end
cands = {name};
if contains(name, '.')
    parts = strsplit(name, '.');
    cands{end+1} = parts{end};
end
for k = 1:numel(cands)
    c = cands{k};
    e = regexpi(text, ['^' regexptranslate('escape', c) ...
        '\s*[-' char(8211) char(8212) ':]\s*'], 'end', 'once');
    if ~isempty(e) && ~isempty(strtrim(text(e+1:end)))
        text = text(e+1:end);
        return
    end
    if numel(c) < 2
        continue
    end
    tok = regexp(text, '^(\S+)\s+(.*)', 'tokens', 'once');
    if ~isempty(tok) && strcmpi(tok{1}, c) && ...
            docgen.Text.isUpper(tok{1}) && ~isempty(tok{2}) && ...
            docgen.Text.isUpper(tok{2}(1))
        text = tok{2};
        return
    end
end
end