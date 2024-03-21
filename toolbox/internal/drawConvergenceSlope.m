function drawConvergenceSlope(x_0,y_0,slope,x_begin,x_end,plotLog)
%drawConvergenceSlope - Draw reference slope in convergence plot
%   drawConvergenceSlope(X_0,Y_0,SLOPE,X_BEGIN,X_END) draws with plot in
%   the current axes the line C*x^(-SLOPE) from x = X_BEGIN to x = X_END,
%   with C chosen such that the line passes through the point (X_0,Y_0).
%   On loglog axes it is a straight line with the slope -SLOPE.
%
%   drawConvergenceSlope(...,PLOTLOG) uses C*x^(-SLOPE)*|log(x)| with
%   C*X_0^(-SLOPE)*|log(X_0)| = Y_0 instead when PLOTLOG is true
%   (default false; X_0 must not be 1). Only the two end points are
%   evaluated and joined by a straight segment, which passes through
%   (X_0,Y_0) only when X_0 is X_BEGIN or X_END.
%
%   See also drawTriangleSlope, plotConvergence
%
%   Note: Create the loglog axes and call hold on first.

% Copyright (c) 2026 by the fem4m authors

if nargin < 6 || ~plotLog
    % N^(-rate) approx
    C = y_0*x_0^(slope);
    plot([x_begin,x_end],C.*[x_begin^(-slope),x_end^(-slope)]);
else
    % N^(-1/2) ~ h <=> N ~ h^(-2)
    % O(h^(2*rate)) reduced to O(h^(2*rate)|ln(h)|)
    % N^(-rate)*|log(N)|
    C = y_0*x_0^(slope)/abs(log(x_0));
    plot([x_begin,x_end],C.*[x_begin^(-slope).*abs(log(x_begin)),x_end^(-slope).*abs(log(x_end))]);

end

end