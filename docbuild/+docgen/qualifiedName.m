function q = qualifiedName(src, name)
%qualifiedName - Package-qualified name of a class or function
%   Q = docgen.qualifiedName(SRC,NAME) returns NAME prefixed with the
%   packages of the source file SRC ('/'-separated path relative to the
%   toolbox), e.g. 'pkg.Grid' for SRC 'mesh/+pkg/Grid.m' and NAME 'Grid';
%   NAME itself if SRC is in no package.
%
%   See also docgen.Lookup, docgen.RefGenerator

% Copyright (c) 2026 by the fem4m authors

folders = strsplit(src, '/');
folders = folders(1:end-1);
pkgs = folders(startsWith(folders, '+'));
q = strjoin([cellfun(@(c) c(2:end), pkgs, 'UniformOutput', false) {name}], '.');
end