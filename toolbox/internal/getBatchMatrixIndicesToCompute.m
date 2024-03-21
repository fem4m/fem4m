function indx = getBatchMatrixIndicesToCompute(resultsDir,N,funcGetFileName)
%getBatchMatrixIndicesToCompute - Experiments without a results file
%   INDX = getBatchMatrixIndicesToCompute(resultsDir,N,funcGetFileName)
%   returns the indices k = 1,...,N, as a column vector, for which the
%   folder resultsDir contains no MAT-file funcGetFileName(k).mat.
%   resultsDir must be an existing folder, e.g.
%   fem4mPath.resultFolder(SUBFOLDER). The function handle
%   funcGetFileName returns the file name of experiment k without
%   extension as a character vector, and the names are compared
%   ignoring case. A batch run then computes only these experiments.
%
%   Example: Experiments of a batch run without saved results
%       r = @(k)fem4mResult(problem='Square',fe='Pk',refine='uniform', ...
%           estimator='residual',solver='direct',param1=k,param2=0, ...
%           maxNrDoF=1e4);
%       name = @(k)fileName(r(k));
%       indx = getBatchMatrixIndicesToCompute(tempdir,3,name)
%
%   See also getBatchMatrix, fem4mResult, fem4mModuleTemplate

% Copyright (c) 2026 by the fem4m authors

files = what(resultsDir);
files = files.mat;
indx = zeros(N,1);
for k = 1 : N
    filename = [funcGetFileName(k),'.mat'];
    filename = repmat(filename,numel(files),1);
    indx(k) = isempty(find(strcmpi(files,filename), 1));
end
indx = find(indx);
end