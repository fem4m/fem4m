function saveFigureHighDpiPng(filename,resolution,format)
%saveFigureHighDpiPng - Save figure as high-resolution image
%   saveFigureHighDpiPng(FILENAME) saves the current figure to
%   FILENAME.png at 300 dpi; FILENAME is a character vector without
%   extension.
%
%   saveFigureHighDpiPng(FILENAME,RESOLUTION) uses RESOLUTION dpi
%   instead ([] for the default 300).
%
%   saveFigureHighDpiPng(...,FORMAT) saves to FILENAME.png,
%   FILENAME.eps or FILENAME.pdf, where FORMAT (case-insensitive) must
%   be:
%       "png" - (default) Save as PNG
%       "eps" - Save as encapsulated PostScript in color
%       "pdf" - Save as PDF
%   For "eps" and "pdf" the figure units and paper units are set to
%   centimeters and the paper size to the size of the figure, and the
%   file holds a raster image at RESOLUTION dpi (print option -image).
%   Any other FORMAT saves nothing, but still changes these figure
%   properties.
%
%   See also print, exportgraphics

% Copyright (c) 2026 by the fem4m authors

if nargin < 2 || isempty(resolution)
    resolution = 300;
end

if nargin < 3 || isempty(format)
    format = 'png';
end

if strcmpi(format,'png')

    print('-dpng','-image', ['-r',num2str(resolution)], [filename,'.png']);

else
    set(gcf, 'PaperUnits','centimeters');
    set(gcf, 'Units','centimeters');
    pos=get(gcf,'Position');
    set(gcf, 'PaperSize', [pos(3) pos(4)]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition',[0 0 pos(3) pos(4)]);

    if strcmpi(format,'eps')

        print('-depsc','-image', ['-r',num2str(resolution)], [filename,'.eps']);

    elseif strcmpi(format,'pdf')

        print('-dpdf','-image', ['-r',num2str(resolution)], [filename,'.pdf']);

    end
end

end