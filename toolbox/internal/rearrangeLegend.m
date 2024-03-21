function rearrangeLegend(neworder,transparent)
%rearrangeLegend - Reorder or fade legend of current axes
%   rearrangeLegend() shows the legend of the current axes.
%
%   rearrangeLegend(NEWORDER) rebuilds the legend with only the graphics
%   objects of the current axes with the indices NEWORDER, counted in the
%   order of creation, in the order of NEWORDER and labeled with their
%   DisplayName ([] keeps the legend as it is).
%
%   rearrangeLegend(...,TRANSPARENT) also removes the legend box and
%   makes the legend background transparent when TRANSPARENT is true
%   (default false).
%
%   See also legend, gca

% Copyright (c) 2026 by the fem4m authors

if nargin < 1
    neworder = [];
end

if nargin < 2 || isempty(transparent)
    transparent = false;
end

hLegend = legend();
hFigure = gca;

% Now re-create the legend
if ~isempty(neworder)
    plots = flipud(get(hFigure, 'children'));
    labels = flipud(get(hFigure.Children, 'DisplayName'));
    legend(plots(neworder), labels(neworder));
end

% set transparent legend
if transparent
    set(hLegend.BoxFace, 'ColorType','truecoloralpha',...
        'ColorData',uint8(255*[1;1;1;0]));
    legend boxoff;
end

end