function mat2data(latexfilename,folder)
%mat2data - Write .dat files of the results referenced in a LaTeX file
%   mat2data(LATEXFILENAME,FOLDER) reads the LaTeX file
%   LATEXFILENAME.tex and writes a data file for every line that
%   contains data/, commented-out lines included, e.g. a reference to
%   data/NAME.dat. NAME has the form MATFILE_GRAPH or
%   MATFILE_GRAPH--BASE: the last '_' separates the name MATFILE of a
%   results file of fem4mResult.save in the subfolder FOLDER of
%   fem4mPath.resultFolder() from the level field GRAPH. The data file
%   has two columns: the level field BASE (default nrDoF) and the level
%   field GRAPH or, for GRAPH of the form A-B, the ratio of the level
%   fields A and B, e.g. an efficiency index. It is written as
%   tab-separated ASCII with double precision to
%   data/MATFILE_GRAPH--BASE.dat; the folder data in the current folder
%   is created if needed.
%
%   Example: Data files of the plots in figures.tex
%       mat2data('figures','run1')
%
%   See also fem4mResult, fem4mPath, readlines
%
%   Note: Only the first reference to data/ on a line is read, and NAME
%   ends at the first '.' after data/. A reference without --BASE is
%   written to data/MATFILE_GRAPH--nrDoF.dat, not to the referenced
%   file.

% Copyright (c) 2026 by the fem4m authors

lines = readlines([latexfilename,'.tex']);
ind = contains(lines,'data/');
lines = lines(ind);

resultsDir = fem4mPath.resultFolder();

if ~exist(fullfile(pwd,'data'),'dir')
    mkdir('data');
end

for j = 1 : numel(lines)

    curLine = lines{j};
    curLine = extractAfter(curLine,['data',filesep]);
    curLine = extractBefore(curLine,'.');
    curLine = strsplit(curLine,'--');
    if numel(curLine)>1
        base = curLine{2};
    else
        base = 'nrDoF';
    end
    curLine = curLine{1};
    indxEnd = strfind(curLine,'_')-1;
    matfilename = curLine(1:indxEnd(end));
    indxBegin = strfind(curLine,'_')+1;
    currentGraph = curLine(indxBegin(end):end);

    m = load(fullfile(resultsDir,folder,[matfilename,'.mat']));
    if contains(currentGraph,'-')
        indx = strfind(currentGraph,'-');
        currentEffNominator = currentGraph(1:indx-1);
        currentEffDenominator = currentGraph(indx+1:end);
        estimatedError = [m.p.level(1:end).(currentEffNominator)];
        eff = [m.p.level(1:end).(currentEffDenominator)];
        eff = (estimatedError./eff)';
        dof = [m.p.level(1:end).(base)];
        data = [dof(:),eff(:)];
    else
        error = [m.p.level(1:end).(currentGraph)];
        dof = [m.p.level(1:end).(base)];
        data = [dof(:),error(:)];
    end

    filenameData = [curLine,'--',base,'.dat'];
    save(fullfile('data',filenameData),'data','-ascii','-double','-tabs');

end

end