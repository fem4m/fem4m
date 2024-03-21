function [a,b,c,d] = intUniqueRows(X)
%intUniqueRows - Stable unique rows as integer indices
%   A = intUniqueRows(X) returns the unique rows of X in the order of
%   their first occurrence, unique(X,'rows','stable').
%
%   [A,B,C] = intUniqueRows(X) also returns the indices B of the first
%   occurrences, A = X(B,:), and C with X = A(C,:), both as uint32
%   column vectors.
%
%   [A,B,C,D] = intUniqueRows(X) also returns the indices D (uint32) of
%   the last occurrences of the rows of A.
%
%   See also unique, ismember

% Copyright (c) 2026 by the fem4m authors

%% new (stable sorting)
[a, b, c] = unique(X,'rows','stable');
b = uint32(b);
c = uint32(c);
if nargout > 3
    [~, perm] = unique(X,'rows','first');
    [~, perm] = sort(perm);
    [~, d] = unique(X,'rows','last');
    d = uint32(d(perm));
end

%% old (sorted)
% [a, b, c] = unique(X,'rows','first');
% if nargout > 3
%     [~, d] = unique(X,'rows','last');
% end

%% matlab code using radixsort
% Z = radixsortrows([X,(1:size(X,1))'],1:size(X,2));
% % count first and last
% b = zeros(size(Z,1),1);
% d = zeros(size(Z,1),1);
% c = zeros(size(Z,1),1);
% indx = 1;
% b(indx) = Z(1,end);
% d(indx) = Z(1,end);
% c(1) = 1;
% for column = 1 : size(Z,1)-1
%    if Z(column,1:end-1) ==  Z(column+1,1:end-1)
%        d(indx) = Z(column+1,end);
%    else
%        indx = indx + 1;
%        b(indx) = Z(column+1,end);
%        d(indx) = Z(column+1,end);
%    end
%    c(column+1) = indx;
% end
% b(indx+1:end) = [];
% d(indx+1:end) = [];
% % permute back
% [b,perm] = sort(b);
% d = d(perm);
% [~,perm] = sort(perm);
% [~,perm2] = sort(Z(:,end));
% c = perm(c(perm2));
% b = uint32(b);
% d = uint32(d);
% c = uint32(c);
% a = uint32(X(b,:));
%
% end
%
% function X = radixsortrows(X,rows)
%     if nargin < 2 || isempty(rows)
%         rows = 1:size(X,2);
%     end
%     maxIndx = max(X(:));
%     % start from the back
%     for k = numel(rows) : -1 : 1
%         row = rows(k);
%         % count
%         sortIndx = accumarray(X(:,row),ones(size(X,1),1),[maxIndx 1]);
%         sortIndx = cumsum([1;sortIndx(1:end-1)]);
%         % collect
%         tmp = zeros(size(X));
%         for column = 1 : size(X,1)
%             tmp(sortIndx(X(column,row)),:) = X(column,:);
%             sortIndx(X(column,row)) = sortIndx(X(column,row)) + 1;
%         end
%         X = tmp;
%     end
% end