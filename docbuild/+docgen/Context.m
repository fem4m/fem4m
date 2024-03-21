classdef Context < handle
    %Context - Folders and toolbox data shared by the documentation generators
    %   CTX = docgen.Context() describes the documentation build of the
    %   repository that contains docbuild: its folders and the data derived
    %   from the toolbox sources. CTX = docgen.Context(DOCFOLDER) uses the
    %   documentation folder DOCFOLDER instead of toolbox/doc ("" for the
    %   default). buildDocs passes one context to all generators, so the
    %   toolbox is scanned once.
    %
    %   Context functions:
    %       hasExamplePage - True if an example has a page
    %       draftFile      - Markdown draft of the overview page of a folder
    %       overviewTitle  - Title of the overview page of a folder
    %       resolve        - The given context, or a new one for a doc folder
    %
    %   Context properties:
    %       Site         - Name of the toolbox in titles, breadcrumbs and footers
    %       Version      - Version of the toolbox in footers
    %   Name and version come from toolbox/Contents.m (docgen.toolboxInfo).
    %       Root         - Repository folder
    %       Toolbox      - Toolbox folder
    %       DocFolder    - Documentation folder (default: toolbox/doc)
    %       Build        - docbuild folder (drafts and stylesheet)
    %       Pages        - Reference pages (docgen.pageTree)
    %       Examples     - Name, title and used names of the examples (docgen.exampleInfo)
    %       ExamplePages - Names of the examples that have a page
    %       Lookup       - Name -> page lookup (docgen.Lookup)
    %
    %   See also buildDocs, docgen.pageTree, docgen.exampleInfo, docgen.Lookup

    % Copyright (c) 2026 by the fem4m authors

    properties (SetAccess = private)
        Root        % repository folder
        Toolbox     % toolbox folder
        DocFolder   % documentation folder, written and read (default: toolbox/doc)
        Build       % docbuild folder (drafts and stylesheet)
        Site        % name of the toolbox in titles, breadcrumbs and footers
        Version     % version of the toolbox in footers ('' if there is none)
        Pages       % reference pages, N-by-5 cell array (docgen.pageTree)
        Examples    % name, title and used names of the examples (docgen.exampleInfo)
    end

    properties (Dependent)
        % names of the examples whose page exists in DocFolder/examples,
        % sorted; read on first use (the example pages are made by
        % exportExamples, before the other pages)
        ExamplePages
        % name -> page lookup, built on first use (docgen.Lookup)
        Lookup
    end

    properties (Access = private)
        ExamplePagesCache = []
        LookupCache = []
    end

    methods

        function obj = Context(docFolder)
            %Context - Scan the toolbox for a documentation build
            arguments
                docFolder (1,1) string = ""
            end
            obj.Build = fileparts(fileparts(mfilename('fullpath')));
            obj.Root = fileparts(obj.Build);
            obj.Toolbox = fullfile(obj.Root, 'toolbox');
            obj.DocFolder = fullfile(obj.Toolbox, 'doc');
            if docFolder ~= ""
                obj.DocFolder = char(docFolder);
            end
            info = docgen.toolboxInfo(obj.Toolbox);
            if isempty(info.name)
                error('docgen:Context:name', ['%s must start with ' ...
                    '''%% <name>: <summary>'''],fullfile(obj.Toolbox, 'Contents.m'));
            end
            obj.Site = info.name;
            obj.Version = info.version;
            obj.Pages = docgen.pageTree(obj.Toolbox);
            obj.Examples = docgen.exampleInfo(obj.Toolbox);
        end

        function names = get.ExamplePages(obj)
            if isnumeric(obj.ExamplePagesCache)  % not read yet ([])
                names = {obj.Examples.name};
                obj.ExamplePagesCache = names(cellfun(@(n) isfile(fullfile( ...
                    obj.DocFolder, 'examples', [n '.html'])), names));
            end
            names = obj.ExamplePagesCache;
        end

        function tf = hasExamplePage(obj, name)
            %hasExamplePage - True if the example NAME has a page
            tf = any(strcmp(obj.ExamplePages, name));
        end

        function file = draftFile(obj, folder)
            %draftFile - Markdown draft of the overview page of a folder
            %   FILE = draftFile(CTX,FOLDER) returns the draft of the
            %   overview page of the source folder FOLDER (relative to the
            %   toolbox, '/'-separated): docbuild/overviews/<main>.md for a
            %   main folder, <main>_<sub>.md for a subfolder <main>/.../<sub>.
            %   The file need not exist.
            [main, rest] = strtok(folder, '/');
            name = main;
            if ~isempty(rest)
                name = [main '_' folder(find(folder == '/', 1, 'last')+1:end)];
            end
            file = fullfile(obj.Build, 'overviews', [name '.md']);
        end

        function title = overviewTitle(obj, folder)
            %overviewTitle - Title of the overview page of a folder
            %   TITLE = overviewTitle(CTX,FOLDER) returns the '# ' heading
            %   of the draft of the source folder FOLDER (see draftFile) as
            %   plain text, i.e. without inline Markdown, and the last
            %   folder name (the heading of a generated page) if there is
            %   no draft or heading.
            title = folder(find(folder == '/', 1, 'last')+1:end);
            file = obj.draftFile(folder);
            if ~isfile(file)
                return
            end
            inCode = false;
            for ln = docgen.Text.readLines(file)
                if startsWith(ln{1}, '```')
                    inCode = ~inCode;
                elseif ~inCode && startsWith(ln{1}, '# ')
                    title = regexprep(strtrim(ln{1}(3:end)), {'`([^`]+)`', ...
                        '\[([^\]]+)\]\([^)]+\)', '\*\*([^*]+)\*\*', ...
                        '(?<!\S)\*([^*]+)\*(?!\S)'}, '$1');
                    return
                end
            end
        end

        function lk = get.Lookup(obj)
            if isempty(obj.LookupCache)
                obj.LookupCache = docgen.Lookup(obj.Pages, obj.ExamplePages);
            end
            lk = obj.LookupCache;
        end

    end

    methods (Static)

        function ctx = resolve(ctx, docFolder)
            %resolve - The given context, or a new one for a doc folder
            %   CTX = docgen.Context.resolve(CTX,DOCFOLDER) returns CTX if
            %   it is not empty, else docgen.Context(DOCFOLDER). The
            %   generators use it for their Context and DocFolder options.
            if isempty(ctx)
                ctx = docgen.Context(docFolder);
            end
        end

    end
end