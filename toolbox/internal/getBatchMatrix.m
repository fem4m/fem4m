function varargout = getBatchMatrix(varargin)
%getBatchMatrix - All combinations of the entries of several vectors
%   [A1,A2,...,An] = getBatchMatrix(V1,V2,...,Vn) returns column vectors
%   of length numel(V1)*numel(V2)*...*numel(Vn) whose rows
%   [A1(k),A2(k),...,An(k)] run through all combinations of entries of
%   V1,...,Vn, with the entries of V1 varying fastest, so that every
%   row describes one experiment of a batch run. Ak has the class of
%   Vk, e.g. a cell array of names. For n >= 2 numeric V1,...,Vn, Ak
%   is Xk(:) for [X1,...,Xn] = ndgrid(V1,...,Vn).
%
%   Example: Combine two finite elements with three problems
%       [feIndx,problemIndx] = getBatchMatrix(1:2,1:3);
%       [feIndx,problemIndx]
%
%   See also getBatchMatrixIndicesToCompute, fem4mModuleTemplate, ndgrid

% Copyright (c) 2026 by the fem4m authors

result = varargin;
for k = 1 : nargin
    result{k} = reshape(result{k},[],1);
end
for k = 1: nargin-1
    result{k+1} = reshape(repmat(result{k+1},1,numel(result{k}))',[],1);
    for l = 1 : k
        result{l} = repmat(result{l},numel(varargin{k+1}),1);
    end
end
varargout = result;
end