function order = memberOrder(names, kinds)
%memberOrder - Order of the classes and functions of a folder
%   ORDER = docgen.memberOrder(NAMES,KINDS) sorts the names in the cell
%   array NAMES, of kind 'class' or 'function' (cell array KINDS): classes
%   first, then functions, each alphabetical ignoring case (ties by
%   name). The table of contents and the member tables of the overview
%   pages list a folder in this order.
%
%   See also genHelptoc, docgen.OverviewGenerator

% Copyright (c) 2026 by the fem4m authors

[~, order] = sort(names);
[~, o] = sort(lower(names(order)));  % sort is stable: ties keep the name order
order = order(o);
[~, o] = sort(~strcmp(kinds(order), 'class'));
order = order(o);
end