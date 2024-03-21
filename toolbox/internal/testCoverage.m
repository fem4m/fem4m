function results = testCoverage()
%testCoverage - Run the tests with a code coverage report
%   RESULTS = testCoverage() runs all tests in the tests folder of the
%   fem4m repository, subfolders included, with text output and writes
%   an HTML code coverage report of the toolbox folders elements,
%   forms, geometries, internal, mesh, problems and solve to the folder
%   coverageReport in the current folder. The plotting, figure export,
%   documentation and repository tools of the internal folder are left
%   out, as are the apps and the examples; the tests do not cover them.
%   RESULTS is the matlab.unittest.TestResult array of the run.
%
%   Example: Run the tests and write the coverage report
%       results = testCoverage();
%
%   See also fem4mPath, testsuite, runtests,
%   matlab.unittest.plugins.CodeCoveragePlugin

% Copyright (c) 2026 by the fem4m authors

testdir = fem4mPath.testFolder;
toolboxdir = fullfile(fem4mPath.root,'toolbox');
codeFolders = ["elements" "forms" "geometries" "internal" "mesh" ...
    "problems" "solve"];
% plotting, figure export, documentation and repository tools
excluded = ["drawConvergenceSlope" "drawTriangleSlope" "plotConvergence" ...
    "plotFunc4RefinedMesh" "plotRefFunc4RefinedMesh" "rearrangeLegend" ...
    "resetFigurePosition" "saveFigureHighDpiPng" "spyEps" "fem4mLogo" ...
    "codeIndentation" "mat2data" "fem4mPath" "testCoverage"] + ".m";
excluded = fullfile(toolboxdir,"internal",excluded);

files = [];
for folder = codeFolders
    listing = dir(fullfile(toolboxdir,folder,"**","*.m"));
    files = [files;string(fullfile({listing.folder},{listing.name}))']; %#ok<AGROW>
end
files = setdiff(files,excluded);

import matlab.unittest.TestRunner
import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoverageReport

suite  = testsuite(testdir, 'IncludeSubfolders', true);
runner = TestRunner.withTextOutput;
runner.addPlugin(CodeCoveragePlugin.forFile(files, ...
    'Producing', CoverageReport('coverageReport')));
results = runner.run(suite);

end