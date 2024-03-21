function info = classMembers(qname, lines)
%classMembers - Documented members of a class from its metaclass
%   INFO = docgen.classMembers(QNAME,LINES) reads the class QNAME with
%   meta.class (the class must be on the path) and returns a struct with
%   fields
%       super    - superclass names joined by ', '
%       abstract - true if the class is abstract
%       props    - struct array of the public properties defined by the
%                  class: name, tags (e.g. {'read-only'}), meta (the
%                  meta.property; read its Description only when needed,
%                  MATLAB parses the help text for it, which is slow)
%       methods  - struct array of the public methods defined by the
%                  class: name, tags (e.g. {'static'}), static, desc
%   Hidden members, inherited members, property get/set methods, the
%   implicit default constructor and, for enumerations, the implicit
%   methods (char, eq, ...) and the member constructor are left out. Methods are listed in the
%   order of the source lines LINES; desc is the summary line of the
%   method (the first comment line after its function line).
%
%   See also meta.class, docgen.RefGenerator

% Copyright (c) 2026 by the fem4m authors

mc = meta.class.fromName(qname);
if isempty(mc)
    error('docgen:classMembers:notFound', ...
        'Class %s not found; is the toolbox on the path?', qname);
end
info.super = strjoin({mc.SuperclassList.Name}, ', ');
info.abstract = mc.Abstract;

% ---- properties: declaration order ----
props = mc.PropertyList;
props = props(arrayfun(@(p) p.DefiningClass == mc && ~p.Hidden && ...
    isPublic(p.GetAccess), props));
info.props = struct('name', {}, 'tags', {}, 'meta', {});
for p = reshape(props, 1, [])
    tags = {};
    if p.Abstract
        tags{end+1} = 'abstract'; %#ok<AGROW>
    end
    if p.Constant
        tags{end+1} = 'constant'; %#ok<AGROW>
    elseif any(strcmp(p.SetAccess, {'private', 'protected', 'immutable'}))
        tags{end+1} = 'read-only'; %#ok<AGROW>
    end
    if p.Dependent
        tags{end+1} = 'dependent'; %#ok<AGROW>
    end
    info.props(end+1) = struct('name', p.Name, 'tags', {tags}, 'meta', p);
end

% ---- methods: source order ----
cname = regexprep(qname, '.*\.', '');
meths = mc.MethodList;
meths = meths(arrayfun(@(m) m.DefiningClass == mc && ~m.Hidden && ...
    isPublic(m.Access), meths));
info.methods = struct('name', {}, 'tags', {}, 'static', {}, 'desc', {});
pos = [];
for m = reshape(meths, 1, [])
    [line, desc] = methodSource(lines, m.Name, m.Abstract);
    if line == 0 && strcmp(m.Name, cname)
        continue  % implicit default constructor
    end
    if mc.Enumeration && (line == 0 || strcmp(m.Name, cname))
        continue  % implicit methods of enumerations (char, eq, ...), member constructor
    end
    tags = {};
    if m.Abstract
        tags{end+1} = 'abstract'; %#ok<AGROW>
    end
    if m.Static
        tags{end+1} = 'static'; %#ok<AGROW>
    end
    info.methods(end+1) = struct('name', m.Name, 'tags', {tags}, ...
        'static', m.Static, 'desc', desc);
    pos(end+1) = line; %#ok<AGROW>
end
[~, order] = sort(pos);
info.methods = info.methods(order);
end

function tf = isPublic(access)
tf = (ischar(access) || isstring(access)) && strcmp(access, 'public');
end

function [line, h1] = methodSource(lines, name, isAbstract)
% Line of the definition of method NAME and its summary line: the first
% comment after the function line, else a 'name - ...' comment right
% above it. An abstract method is a signature in a methods (Abstract)
% block; its summary is the comment above (or below) the signature.
line = 0;
h1 = '';
fnRe = ['^\s*function\s+(\[[^\]]*\]\s*=\s*|[\w~]+\s*=\s*)?' ...
    regexptranslate('escape', name) '\s*(\(|\.\.\.|$)'];
for i = docgen.Help.functionLines(lines)
    code = docgen.Help.stripComment(lines{i});
    if ~isempty(regexp(code, fnRe, 'once'))
        line = i;
        break
    end
end
if line > 0
    last = line;  % signature continued with '...'
    while endsWith(deblank(docgen.Help.stripComment(lines{last})), ...
            '...') && last < numel(lines)
        last = last + 1;
    end
    k = last + 1;
    while k <= numel(lines) && startsWith(strtrim(lines{k}), '%')
        t = commentText(lines{k});
        if ~isempty(t)
            h1 = t;
            return
        end
        k = k + 1;
    end
    above = commentAbove(lines, line);
    if ~isempty(above) && ~isempty(regexpi(above, ...
            ['^' regexptranslate('escape', name) '\s*[-:]'], 'once'))
        h1 = above;
    end
    return
end
if ~isAbstract
    return
end
inBlock = false;
sigRe = ['^\s*(\[[^\]]*\]\s*=\s*|\w+\s*=\s*)?' ...
    regexptranslate('escape', name) '\s*(\([^)]*\))?\s*;?\s*$'];
for i = 1:numel(lines)
    code = docgen.Help.stripComment(lines{i});
    if ~isempty(regexp(code, '^\s*methods\s*\(.*Abstract', 'once'))
        inBlock = true;
    elseif inBlock && ~isempty(regexp(code, sigRe, 'once'))
        line = i;
        h1 = commentAbove(lines, i);
        if isempty(h1) && i < numel(lines) && ...
                startsWith(strtrim(lines{i+1}), '%')
            h1 = commentText(lines{i+1});
        end
        return
    end
end
end

function t = commentText(ln)
t = strtrim(docgen.Help.uncomment(strtrim(ln)));
end

function t = commentAbove(lines, i)
% nearest non-empty comment line above line I, through comment lines only
t = '';
k = i - 1;
while k >= 1 && startsWith(strtrim(lines{k}), '%')
    t = commentText(lines{k});
    if ~isempty(t)
        return
    end
    k = k - 1;
end
end