function info = toolboxInfo(toolbox)
%toolboxInfo - Name, summary and version of the toolbox
%   INFO = docgen.toolboxInfo(TOOLBOX) reads TOOLBOX/Contents.m, the one
%   place of the version. Its first two lines are, in the format of ver,
%       % fem4m: A Finite Element Toolbox for MATLAB
%       % Version 1.0 (R2026a) 2026-09-28
%   with ':' or ' - ' after the name. INFO is a struct with the fields
%   name ('fem4m'), summary ('A Finite Element Toolbox for MATLAB') and
%   version ('1.0'); a field is '' if its line is missing.
%
%   See also docgen.Context, buildtool

% Copyright (c) 2026 by the fem4m authors

arguments
    toolbox (1,:) char
end
info = struct('name', '', 'summary', '', 'version', '');
lines = docgen.Text.readLines(fullfile(toolbox, 'Contents.m'));
if ~isempty(lines)
    tok = regexp(lines{1}, '^%\s*([^\s:]+)(?::|\s+-)\s+(.*\S)', 'tokens', 'once');
    if ~isempty(tok)
        [info.name, info.summary] = tok{:};
    end
end
if numel(lines) > 1
    tok = regexp(lines{2}, '^%\s*Version\s+(\S+)', 'tokens', 'once');
    if ~isempty(tok)
        info.version = tok{1};
    end
end
end