function plotConvergence(nrDoF4lvl, error4lvl, name)
%plotConvergence - Plot error curve in loglog scale
%   plotConvergence(NRDOF4LVL,ERROR4LVL) plots ERROR4LVL against
%   NRDOF4LVL with loglog and the line style '-x' in the current axes,
%   sets the title 'Convergence Plot' and shows the legend. With at least
%   two levels the legend entry is the convergence rate, the negative
%   least-squares slope of log(ERROR4LVL) over log(NRDOF4LVL), in
%   parentheses.
%
%   plotConvergence(NRDOF4LVL,ERROR4LVL,NAME) puts the character vector
%   NAME in front of the rate in the legend entry (default '').
%
%   See also drawConvergenceSlope, drawTriangleSlope
%
%   Note: Call hold on to plot several curves in the same axes.

% Copyright (c) 2026 by the fem4m authors

if nargin < 3 || isempty(name)
    name = '';
end

% display convergence rate in legend
if numel(error4lvl)>1
    p = polyfit(log(nrDoF4lvl),log(error4lvl),1);
    name = char([name ' (' num2str(-p(1)) ')']);
end

% loglog plot
loglog(nrDoF4lvl,error4lvl,'-x','DisplayName',name);
title('Convergence Plot');
legend('off');
legend('show');

end