function plan = buildfile
%buildfile - Build plan of the fem4m toolbox
%   Run the tasks with buildtool from the repository folder:
%
%       buildtool                     % pages and test (default)
%       buildtool docs                % complete documentation: examples, pages, search index
%       buildtool pages               % only the reference, overview and TOC pages
%       buildtool examples            % run the examples, export them to HTML
%       buildtool searchindex         % search index of the help browser
%       buildtool test                % run the tests in tests/
%       buildtool package             % docs + test, then release/fem4m-<version>.mltbx
%       buildtool clean               % delete toolbox/doc and release/ (the task outputs)
%       buildtool -tasks              % list all tasks
%
%   The documentation is generated and checked by docbuild/buildDocs.m
%   (genOverviews, genRef, genExamplesOverview, genHelptoc, checkDocs); the
%   example pages by docbuild/exportExamples.m. toolbox/doc is build
%   output and not under version control: run "buildtool docs" for the
%   complete documentation.
%
%   The examples, pages, searchindex and package tasks are incremental: buildtool skips
%   a task whose inputs and outputs have not changed since its last
%   successful run (it compares file contents and keeps its records in
%   .buildtool/). Changing any toolbox .m file runs all examples again.
%
%   The version is set in one place, line 2 of toolbox/Contents.m
%   ('% Version 1.0 (R2026a) 2026-09-28'); the package and the footer of
%   the documentation pages take it from there.
%   To force a task, delete its outputs first, e.g.
%   buildtool clean("examples") examples.
%
%   See also buildtool, matlab.buildtool.Plan

% Copyright (c) 2026 by the fem4m authors

plan = buildplan(localfunctions);

plan("test") = matlab.buildtool.tasks.TestTask("tests");
plan("clean") = matlab.buildtool.tasks.CleanTask;

% examples first: the reference and overview pages link only to example
% pages that exist, and the search index covers all pages. buildtool runs
% independent dependencies in the order of the task functions in this
% file, not in the listed order, so examplesTask is defined before
% pagesTask.
plan("searchindex").Dependencies = ["examples" "pages"];
plan("docs") = matlab.buildtool.Task( ...
    Description="Build the complete documentation: examples, pages and search index", ...
    Dependencies="searchindex");
plan("package").Dependencies = ["docs" "test"];

% ---- inputs and outputs of the incremental tasks ----
code = "toolbox/**/*.m";  % the toolbox and its examples (doc has no .m files)
% the same without the files directly in toolbox: Contents.m holds the
% version, which is not part of the example pages
subfolderCode = "toolbox/*/**/*.m";
doc = "toolbox/doc";
examples = dir(fullfile(plan.RootFolder, "toolbox", "examples", "*.m"));
[~, stems] = fileparts(string({examples.name}));
examplePages = doc + "/examples/" + stems + ".html";
mains = pageFolders(plan.RootFolder);

% exportExamples also sets the page titles (exampleInfo, exampleTitle)
% and links the names in the text (Lookup)
plan("examples").Inputs = [subfolderCode "docbuild/exportExamples.m" ...
    "docbuild/+docgen/exampleInfo.m" "docbuild/+docgen/exampleTitle.m" ...
    "docbuild/+docgen/Lookup.m"];
plan("examples").Outputs = examplePages;
% the example pages are inputs: the other pages link only to those that
% exist. As a pattern, missing pages (e.g. after buildtool clean) are
% allowed, and plan("examples").Outputs is not used since an input that
% refers to the outputs of a task makes pages depend on examples.
plan("pages").Inputs = [code "docbuild/**/*.m" "docbuild/doc.css" ...
    "docbuild/overviews/*.md" doc + "/examples/*.html"];
plan("pages").Outputs = [doc + ["/index.html" "/getting_started.html" ...
    "/helptoc.xml" "/css" "/examples/examples_overview.html"], doc + "/" + mains];
plan("searchindex").Inputs = [doc + "/**/*.html" "toolbox/info.xml"];
plan("searchindex").Outputs = doc + "/helpsearch-v4_en";
% the packaged files with the README and its images; as an output the
% package is deleted by buildtool clean
plan("package").Inputs = ["toolbox" "README.md" ".github/images"];
plan("package").Outputs = "release/*.mltbx";

plan.DefaultTasks = ["pages" "test"];
end

function examplesTask(context)
% Run all examples and export them with their output to HTML
restorePath = addDocbuild(context); %#ok<NASGU>
results = exportExamples();
failed = results(~[results.ok]);
if ~isempty(failed)
    error("fem4m:build:examples", "Examples failed: %s", ...
        strjoin(string({failed.name}), ", "));
end
end

function pagesTask(context)
% Generate the HTML reference, overview and TOC pages and check the links
restorePath = addDocbuild(context); %#ok<NASGU>
result = buildDocs();
% checkDocs only reports; fail the build on broken links or images and
% on pages that are not in the table of contents
problems = [result.broken; result.missingImages]';  % 2-by-N: page, link
if ~isempty(problems)
    error("fem4m:build:pages", "%d broken links or missing images " + ...
        "(page -> link):%s", size(problems, 2), ...
        sprintf("\n  %s -> %s", problems{:}));
end
toc = [strcat(result.notInToc, " (page not in the TOC)"), ...
    strcat(result.tocMissing, " (TOC target without file)"), ...
    strcat(result.tocDuplicates, " (TOC target listed twice)")];
if ~isempty(toc)
    error("fem4m:build:pages", "%d problems with helptoc.xml:%s", ...
        numel(toc), sprintf("\n  %s", toc{:}));
end
end

function searchindexTask(context)
% Rebuild the search index of the help browser
toolboxFolder = fullfile(context.Plan.RootFolder, "toolbox");
oldPath = path;
restorePath = onCleanup(@() path(oldPath));
addpath(toolboxFolder);  % builddocsearchdb needs info.xml on the path
builddocsearchdb(fullfile(toolboxFolder, "doc"));
end

function packageTask(context)
% Package the toolbox as release/fem4m-<version>.mltbx
% with the name, summary and version of toolbox/Contents.m
root = context.Plan.RootFolder;
toolboxFolder = fullfile(root, "toolbox");
restorePath = addDocbuild(context); %#ok<NASGU>
info = docgen.toolboxInfo(toolboxFolder);
if info.version == "" || info.summary == ""
    error("fem4m:build:package", "toolbox/Contents.m must start with " + ...
        "'%% <name>: <summary>' and '%% Version <version> ...'");
end

% the README of the repository is the readme of the package (Readme,
% since R2026b): staged with its images in the toolbox folder before
% the files are collected, links relative to the toolbox folder
hasReadme = ismember("Readme", properties("matlab.addons.toolbox.ToolboxOptions"));
readme = fullfile(toolboxFolder, "README.md");
readmeImages = fullfile(toolboxFolder, "doc", "readme");
if hasReadme
    if isfile(readme)
        error("fem4m:build:readme", "toolbox/README.md exists; the package task writes it.");
    end
    removeReadme = onCleanup(@() removeStagedReadme(readme, readmeImages));
    copyfile(fullfile(root, ".github", "images"), readmeImages);
    text = fileread(fullfile(root, "README.md"), "Encoding", "UTF-8");
    text = strrep(text, ".github/images/", "doc/readme/");
    text = strrep(text, "](toolbox/", "](");
    text = strrep(text, "href=""toolbox/", "href=""");
    writelines(text, readme, "Encoding", "UTF-8");
end

% fixed identifier: an installed fem4m is updated by a newer release
opts = matlab.addons.toolbox.ToolboxOptions(toolboxFolder, ...
    "be0cfd1f-8db3-4705-8670-a6df90e81928");
opts.ToolboxName = info.name;
opts.ToolboxVersion = info.version;
opts.AuthorName = "Joscha Gedicke";
opts.Summary = info.summary;
opts.Description = "Finite element spaces (Lagrange, DG, hp, " + ...
    "nonconforming, H(div), H(curl), plate and elasticity elements), " + ...
    "variational forms written as strings, adaptive triangular, " + ...
    "tetrahedral, polygonal and surface meshes, problem definitions " + ...
    "for many PDEs, and multigrid and domain decomposition solvers.";
opts.MinimumMatlabRelease = "R2026a";  % search index built with R2026a
% the logo, shown when installing and in the Add-On Manager
opts.ToolboxImageFile = fullfile(root, ".github", "images", "fem4m-logo.png");
if hasReadme
    opts.Readme = readme;
end

% leave out editor backups, Finder files and the machine-specific doc
% registration that MATLAB writes into doc/
files = opts.ToolboxFiles;
exclude = endsWith(files, [".asv" ".DS_Store"]) | ...
    endsWith(files, filesep + "custom_toolbox.json");
opts.ToolboxFiles = files(~exclude);

% the code folders go on the path, the doc folders do not
paths = opts.ToolboxMatlabPath;
docFolder = fullfile(toolboxFolder, "doc");
opts.ToolboxMatlabPath = paths(paths ~= docFolder & ...
    ~startsWith(paths, docFolder + filesep));

releaseFolder = fullfile(root, "release");
if ~isfolder(releaseFolder)
    mkdir(releaseFolder);
end
opts.OutputFile = fullfile(releaseFolder, info.name + "-" + info.version + ".mltbx");
matlab.addons.toolbox.packageToolbox(opts);
fprintf("Packaged %s\n", opts.OutputFile);
end

function removeStagedReadme(readme, images)
% Remove the README and its images staged for the package
if isfile(readme)
    delete(readme);
end
if isfolder(images)
    rmdir(images, "s");
end
end

function restorePath = addDocbuild(context)
% Put docbuild/ on the path until the returned object is cleared
oldPath = addpath(fullfile(context.Plan.RootFolder, "docbuild"));
restorePath = onCleanup(@() path(oldPath));
end

function mains = pageFolders(root)
% Main folders of the toolbox with reference pages (docgen.pageTree),
% each written to a folder of the same name in toolbox/doc
oldPath = addpath(fullfile(root, "docbuild"));
restorePath = onCleanup(@() path(oldPath));
pages = docgen.pageTree(fullfile(root, "toolbox"));
mains = unique(string(pages(:, 5)))';
end