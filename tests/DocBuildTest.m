classdef DocBuildTest < matlab.unittest.TestCase
    %DocBuildTest - Build the documentation and check its links and contents
    %   DocBuildTest generates the documentation (docbuild/buildDocs.m)
    %   into a temporary folder, so toolbox/doc is not touched, and checks
    %   that there are no broken links or missing images, that every page
    %   is in the table of contents and every class and function has its
    %   page, and that every example has a title and a topic in
    %   docbuild/overviews/examples.md. The example pages are not built
    %   (that runs all examples); the other pages then do not link to them.
    %
    %   See also buildDocs, checkDocs

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = private)
        Result   % checkDocs result of the build
        Context  % docgen.Context of the build
    end

    methods (TestClassSetup)

        function buildIntoTemporaryFolder(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(root, 'docbuild')));
            folder = testCase.applyFixture( ...
                matlab.unittest.fixtures.TemporaryFolderFixture).Folder;
            result = [];
            evalc('result = buildDocs(DocFolder=folder);');  % without the build log
            testCase.Result = result;
            testCase.Context = docgen.Context(folder);
        end

    end

    methods (Test)

        function testNoBrokenLinks(testCase)
            testCase.verifyEmpty(testCase.Result.broken, ...
                listing('Broken links (page -> link):', testCase.Result.broken));
        end

        function testNoMissingImages(testCase)
            testCase.verifyEmpty(testCase.Result.missingImages, ...
                listing('Missing images (page -> image):', testCase.Result.missingImages));
        end

        function testEveryPageInToc(testCase)
            testCase.verifyEmpty(testCase.Result.notInToc, ...
                listing('Pages without TOC entry:', testCase.Result.notInToc(:)));
        end

        function testTocTargetsExistOnce(testCase)
            testCase.verifyEmpty(testCase.Result.tocMissing, ...
                listing('TOC targets without file:', testCase.Result.tocMissing(:)));
            testCase.verifyEmpty(testCase.Result.tocDuplicates, ...
                listing('TOC targets listed twice:', testCase.Result.tocDuplicates(:)));
        end

        function testEveryClassAndFunctionHasPage(testCase)
            testCase.verifyEmpty(testCase.Result.unreferenced, ...
                listing('Missing pages (source -> page):', testCase.Result.unreferenced));
        end

        function testEveryExampleHasTitle(testCase)
            info = testCase.Context.Examples;
            untitled = {info(cellfun(@isempty, {info.title})).file};
            testCase.verifyEmpty(untitled, ...
                listing('Examples without a ''# Title'' heading:', untitled(:)));
        end

        function testEveryExampleHasTopic(testCase)
            topics = docgen.exampleTopics(testCase.Context);
            testCase.verifyEmpty(topics.untopiced, listing(['Examples ' ...
                'without a topic in docbuild/overviews/examples.md:'], ...
                topics.untopiced(:)));
            testCase.verifyEmpty(topics.unknown, listing(['Names in ' ...
                'docbuild/overviews/examples.md that are not examples:'], ...
                topics.unknown(:)));
        end

    end

end

function s = listing(heading, items)
% diagnostic text: HEADING and one line per row of the cell array ITEMS
rows = cell(size(items, 1), 1);
for k = 1:size(items, 1)
    rows{k} = strjoin(items(k, :), ' -> ');
end
s = strjoin([{heading}; rows], newline);
end