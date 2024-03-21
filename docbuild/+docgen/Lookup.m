classdef Lookup < handle
    %Lookup - Map class and function names to their documentation pages
    %   OBJ = docgen.Lookup(PAGES,EXAMPLES) builds the lookup from the
    %   reference pages PAGES (see docgen.pageTree) and the names of the
    %   examples with a page EXAMPLES (cell array). Every documented class
    %   or function is known by its name (and by its package-qualified
    %   name for package members); example pages are known by their name
    %   in the page folder 'examples'. A name found in more than one place
    %   is ambiguous: it is not linked and collected in Ambiguous.
    %
    %   Lookup functions:
    %       entries   - Pages known for a name
    %       resolve   - Page folder and file of a uniquely known name
    %       pathFor   - Relative link to the page of a name
    %       linkNames - Link every known name in a list of names
    %
    %   Lookup properties:
    %       Ambiguous - Names that were looked up but are not unique
    %
    %   See also docgen.Context, docgen.RefGenerator

    % Copyright (c) 2026 by the fem4m authors

    properties
        Ambiguous = cell(1, 0)  % names that were looked up but are not unique
    end

    properties (Access = private)
        Map  % name -> struct array with fields dir, page, src
    end

    methods

        function obj = Lookup(pages, examples)
            %Lookup - Build the name lookup from the pages
            obj.Map = configureDictionary("string", "cell");
            for i = 1:size(pages, 1)
                [page, src, ~, ~, main] = pages{i, :};
                stem = docgen.Text.stem(src);
                qual = docgen.qualifiedName(src, stem);  % 'mesh/+pkg/Grid.m' -> 'pkg.Grid'
                if ~strcmp(qual, stem)
                    obj.add(qual, main, page, src);
                end
                obj.add(stem, main, page, src);
            end
            for i = 1:numel(examples)
                obj.add(examples{i}, 'examples', [examples{i} '.html'], ...
                    ['examples/' examples{i} '.m']);
            end
        end

        function ents = entries(obj, name)
            %entries - Pages known for a name
            %   ENTS = entries(OBJ,NAME) returns a struct array with fields
            %   dir, page and src (empty if NAME is unknown).
            if isKey(obj.Map, name)
                ents = obj.Map{name};
            else
                ents = struct('dir', {}, 'page', {}, 'src', {});
            end
        end

        function [dir, page] = resolve(obj, name, excludeSrc)
            %resolve - Page folder and file of a uniquely known name
            %   [DIR,PAGE] = resolve(OBJ,NAME,EXCLUDESRC) returns the page
            %   of NAME, ignoring the source file EXCLUDESRC (so a page
            %   does not link to itself). DIR is empty if NAME is unknown
            %   or ambiguous.
            dir = '';
            page = '';
            ents = obj.entries(name);
            if nargin > 2 && ~isempty(excludeSrc) && ~isempty(ents)
                ents = ents(~strcmp({ents.src}, excludeSrc));
            end
            if isempty(ents)
                return
            end
            keys = unique(strcat({ents.dir}, '/', {ents.page}));
            if isscalar(keys)
                dir = ents(1).dir;
                page = ents(1).page;
            elseif ~any(strcmp(obj.Ambiguous, name))
                obj.Ambiguous{end+1} = name;
            end
        end

        function href = pathFor(obj, name, pageDir, excludeSrc)
            %pathFor - Relative link to the page of a name
            %   HREF = pathFor(OBJ,NAME,PAGEDIR,EXCLUDESRC) returns the link
            %   to the page of NAME as seen from a page in folder PAGEDIR,
            %   or '' if NAME is not uniquely known.
            if nargin < 4
                excludeSrc = '';
            end
            [dir, page] = obj.resolve(name, excludeSrc);
            if isempty(dir)
                href = '';
            elseif isempty(pageDir)
                href = [dir '/' page];
            elseif strcmp(dir, pageDir)
                href = page;
            else
                href = ['../' dir '/' page];
            end
        end

        function out = linkNames(obj, text, pageDir, code)
            %linkNames - Link every known name in a list of names
            %   OUT = linkNames(OBJ,TEXT,PAGEDIR,CODE) splits TEXT at the
            %   separators '/ , ; & and', links each known name to its page
            %   and keeps the separators. With CODE true the names are set
            %   in <code>.
            sep = '\s*(?:/|,|;|&|\<and\>)\s*';
            [parts, seps] = regexp(text, sep, 'split', 'match');
            out = '';
            for k = 1:numel(parts)
                part = parts{k};
                if ~isempty(part)
                    name = strtrim(part);
                    href = obj.pathFor(name, pageDir);
                    if isempty(href)
                        out = [out docgen.Text.escape(part)]; %#ok<AGROW>
                    else
                        inner = docgen.Text.escape(name);
                        if code
                            inner = ['<code>' inner '</code>']; %#ok<AGROW>
                        end
                        out = [out '<a href="' href '">' inner '</a>']; %#ok<AGROW>
                    end
                end
                if k <= numel(seps)
                    out = [out docgen.Text.escape(seps{k})]; %#ok<AGROW>
                end
            end
        end

    end

    methods (Access = private)

        function add(obj, name, dir, page, src)
            ent = struct('dir', dir, 'page', page, 'src', src);
            if isKey(obj.Map, name)
                ents = obj.Map{name};
                same = strcmp({ents.dir}, dir) & strcmp({ents.page}, page) ...
                    & strcmp({ents.src}, src);
                if ~any(same)
                    obj.Map{name} = [ents ent];
                end
            else
                obj.Map{name} = ent;
            end
        end

    end
end