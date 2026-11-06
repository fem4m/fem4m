classdef fem4matPath
    %fem4matPath - Locate and prepare the fem4mat repository
    %   fem4matPath provides static methods that return the fem4mat root
    %   folder, add all of its folders except those starting with '.'
    %   to the MATLAB path, and return the results and tests folders,
    %   creating them if needed. The root folder is the repository
    %   folder that contains toolbox/internal/fem4matPath.m, whatever its
    %   name.
    %
    %   fem4matPath functions:
    %       root         - Path of the fem4mat folder
    %       addpath      - Add all fem4mat folders to the path
    %       resultFolder - Return and create results folder
    %       testFolder   - Return and create tests folder
    %
    %   See also fem4matResult, testCoverage, codeIndentation, genpath

    % Copyright (c) 2026 by the fem4mat authors

    methods (Static)

        function path = root()
            %ROOT - Path of the fem4mat folder
            %   PATH = fem4matPath.ROOT() returns the path of the fem4mat
            %   repository folder, the folder that contains
            %   toolbox/internal/fem4matPath.m.
            path = mfilename('fullpath');
            % this file lives in <root>/toolbox/internal
            path = fileparts(fileparts(fileparts(path)));
        end

        function addpath()
            %ADDPATH - Add all fem4mat folders to the path
            %   fem4matPath.ADDPATH() adds every fem4mat folder to the
            %   MATLAB path except the folders whose name starts with '.',
            %   e.g. .git or .buildtool, and their subfolders. As in
            %   genpath, class, package, private and resources folders are
            %   not added; MATLAB finds class, package and private folders
            %   through their parent folders.
            path = fem4matPath.root();
            allpaths = split(string(genpath(path)),pathsep);
            allpaths = allpaths(strlength(allpaths) > 0);
            % only the part below the root counts, the root itself may
            % lie in a folder starting with '.'
            hidden = contains(extractAfter(allpaths,strlength(path)),filesep+".");
            addpath(strjoin(allpaths(~hidden),pathsep));
        end

        function resultsDir = resultFolder(subfolder)
            %resultFolder - Return and create results folder
            %   RESULTSDIR = fem4matPath.resultFolder() returns the
            %   results folder of the fem4mat folder, creating it if
            %   needed.
            %
            %   RESULTSDIR = fem4matPath.resultFolder(SUBFOLDER) returns
            %   the subfolder SUBFOLDER of the results folder, creating
            %   it if needed.
            path = fem4matPath.root();
            if nargin > 0
                resultsDir = fullfile(path,'results',subfolder);
            else
                resultsDir = fullfile(path,'results');
            end
            if ~exist(resultsDir,'dir')
                mkdir(resultsDir);
            end
        end

        function testDir = testFolder(subfolder)
            %testFolder - Return and create tests folder
            %   TESTDIR = fem4matPath.testFolder() returns the tests
            %   folder of the fem4mat folder, creating it if needed.
            %
            %   TESTDIR = fem4matPath.testFolder(SUBFOLDER) returns the
            %   subfolder SUBFOLDER of the tests folder, creating it if
            %   needed.
            path = fem4matPath.root();
            if nargin > 0
                testDir = fullfile(path,'tests',subfolder);
            else
                testDir = fullfile(path,'tests');
            end
            if ~exist(testDir,'dir')
                mkdir(testDir);
            end
        end

    end
end