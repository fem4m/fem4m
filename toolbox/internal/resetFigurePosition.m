function resetFigurePosition(handle,opts)
%resetFigurePosition - Restore default figure position
%   resetFigurePosition() sets the position of the current figure to its
%   default position (property defaultFigurePosition).
%
%   resetFigurePosition(HANDLE) resets the figure HANDLE instead ([] for
%   the current figure).
%
%   resetFigurePosition(...,OPTS) chooses the layout:
%       "inside" - (default) Default position and size
%       other    - Default position, widened by the width of the legend
%                  of the current axes, which moves to the location
%                  'northeastoutside'
%
%   See also figure, legend
%
%   Note: The legend is that of the current axes, also when HANDLE is
%   another figure.

% Copyright (c) 2026 by the fem4m authors

if nargin < 1 || isempty(handle)
    handle = gcf;
end

if nargin < 2 || isempty(opts)
    opts = 'inside';
end

if strcmpi(opts,'inside')
    defaultPosition = get(handle,'defaultFigurePosition');
    set(handle,'Position',defaultPosition);
else
    defaultPosition = get(handle,'defaultFigurePosition');
    leg = legend();
    set(leg,'Units','pixels');
    legPosition = get(leg,'Position');
    position = defaultPosition;
    position(3) = position(3) + legPosition(3);
    set(handle,'Position',position);
    set(leg,'Location','northeastoutside');
end