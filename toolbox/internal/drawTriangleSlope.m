function drawTriangleSlope(type,slope,x_coord,y_coord,numbers,c)
%drawTriangleSlope - Draw slope triangle in loglog plot
%   drawTriangleSlope(TYPE,SLOPE,X_COORD,Y_COORD) draws with loglog in
%   the current axes a right triangle whose hypotenuse has the slope
%   -SLOPE (decay rate SLOPE) and whose horizontal leg runs from X_COORD
%   to 10*X_COORD. The triangle gets no legend entry. TYPE must be:
%       "top"    - Triangle above the hypotenuse with the horizontal leg
%                  on top; (X_COORD,Y_COORD) is its upper left corner
%       "bottom" - Triangle below the hypotenuse with the horizontal leg
%                  at the bottom; (X_COORD,Y_COORD) is its lower left
%                  corner
%
%   drawTriangleSlope(TYPE,SLOPE,X_COORD,Y_COORD,NUMBERS) also labels the
%   horizontal leg with 1 and the vertical leg with SLOPE (as a fraction
%   for 1/2, 1/3, 1/4, 2/3, 3/4, 4/3 and 5/4) when NUMBERS is true
%   (default false).
%
%   drawTriangleSlope(...,NUMBERS,C) uses the factor C for the
%   horizontal leg, which then runs from X_COORD to C*X_COORD (default
%   10).
%
%   See also drawConvergenceSlope, plotConvergence
%
%   Note: Call hold on first; otherwise loglog replaces the plot.

% Copyright (c) 2026 by the fem4m authors

if nargin < 5 || isempty(numbers)
    numbers = false;
else
    numbers = logical(numbers);
end

if nargin < 6
    c = 10;
end

FontSize = 14;

if slope == 1/2
    str_slope = '1/2';
elseif slope == 1/3
    str_slope = '1/3';
elseif slope == 1/4
    str_slope = '1/4';
elseif slope == 2/3
    str_slope = '2/3';
elseif slope == 3/4
    str_slope = '3/4';
elseif slope == 4/3
    str_slope = '4/3';
elseif slope == 5/4
    str_slope = '5/4';
else
    str_slope = num2str(slope);
end

if strcmp(type,'bottom')
    orient = 1;
elseif strcmp(type,'top')
    orient = 2;
end

drawDreieck(slope,str_slope,x_coord,y_coord,c,orient,numbers,FontSize);

end


function drawDreieck(slope,str_slope,a1,b1,c,orient,numbers,FontSize)
%drawDreieck - Draw and annotate slope triangle
%   drawDreieck(SLOPE,STR_SLOPE,A1,B1,C,ORIENT,NUMBERS,FontSize) draws the
%   triangle of drawTriangleSlope with the corner (A1,B1), the factor C
%   and ORIENT = 1 ("bottom") or 2 ("top"), hides it from the legend and,
%   when NUMBERS is true, writes 1 and STR_SLOPE at the legs in the font
%   size FontSize.

if orient == 1
    a2 = c * a1; b2 = (c)^(slope) * b1;
    handle = loglog([a1;a2;a1;a1],[b1;b1;b2;b1],'k','LineWidth',1.0);
    x1 = exp((log(a1)+log(a2))/2);
    y1 = b1*0.75;
    x2 = a1*0.75;
    y2 = exp((log(b1)+log(b2))/2);
else
    a2 = c * a1; b2 = (c)^(-slope) * b1;
    handle = loglog([a1;a2;a2;a1],[b1;b1;b2;b1],'k','LineWidth',1.0);
    x1 = exp((log(a1)+log(a2))/2);
    y1 = b1*1.5;
    x2 = a2*1.1;
    y2 = exp((log(b1)+log(b2))/2);
end

hasbehavior(handle,'legend',false);

if numbers
    text(x1,y1,'1','FontSize',FontSize)
    text(x2,y2,str_slope,'FontSize',FontSize)
end

end