function val = AitkenExtrapolation(s4l)
%AitkenExtrapolation - Aitken extrapolation of a sequence
%   VAL = AitkenExtrapolation(S4L) applies the Aitken delta-squared
%   process to each row of the M-by-N array S4L, one sequence per row (a
%   vector is taken as one sequence), and repeats it on the result, each
%   pass removing two columns, until fewer than three columns remain. VAL
%   (M-by-1) is the last remaining column.
%
%   See also plotConvergence, drawConvergenceSlope
%
%   Note: A zero second difference gives Inf or NaN entries.

% Copyright (c) 2026 by the fem4m authors

if numel(s4l) == max(size(s4l))
    s4l = s4l(:)';
end

while size(s4l,2) >= 3
    N = size(s4l,2);
    s4l = s4l(:,1:N-2)- (s4l(:,2:N-1)-s4l(:,1:N-2)).^2./(s4l(:,3:N)-2*s4l(:,2:N-1)+s4l(:,1:N-2));
end

val = s4l(:,end);