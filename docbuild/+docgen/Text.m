classdef Text
    %Text - String helpers shared by the documentation generators
    %
    %   Text functions:
    %       readLines  - Lines of a UTF-8 text file
    %       splitLines - Split text at line breaks
    %       writeText  - Write text to a UTF-8 file
    %       indent     - Number of leading whitespace characters
    %       isUpper    - True for words like 'DEGREE' or 'S4L'
    %       escape     - Escape &, <, >, " and ' for HTML
    %       stem       - File name without folder and extension
    %
    %   See also docgen.Help, docgen.Lookup

    % Copyright (c) 2026 by the fem4m authors

    methods (Static)

        function lines = readLines(file)
            %readLines - Lines of a UTF-8 text file
            %   LINES = docgen.Text.readLines(FILE) returns the lines of
            %   FILE as a row cell array of char, without line breaks; a
            %   final line break does not start an empty line.
            lines = docgen.Text.splitLines(fileread(file, 'Encoding', 'UTF-8'));
        end

        function lines = splitLines(txt)
            %splitLines - Split text at line breaks
            %   LINES = docgen.Text.splitLines(TXT) splits TXT at \n, \r\n
            %   and \r; a final line break does not start an empty line.
            if isempty(txt)
                lines = cell(1, 0);
                return
            end
            lines = regexp(txt, '\r\n|[\n\r]', 'split');
            if isempty(lines{end})
                lines(end) = [];
            end
        end

        function writeText(file, txt)
            %writeText - Write text to a UTF-8 file
            %   docgen.Text.writeText(FILE,TXT) writes the char vector TXT
            %   to FILE, replacing its contents.
            fid = fopen(file, 'w', 'n', 'UTF-8');
            if fid < 0
                error('docgen:Text:write', 'Cannot write %s', file);
            end
            closer = onCleanup(@() fclose(fid));
            fwrite(fid, txt, 'char');
        end

        function n = indent(s)
            %indent - Number of leading whitespace characters
            n = numel(regexp(s, '^\s*', 'match', 'once'));
        end

        function tf = isUpper(s)
            %isUpper - True for words like 'DEGREE' or 'S4L'
            %   TF = docgen.Text.isUpper(S) is true if S has at least one
            %   letter and no lower-case letter.
            tf = any(isstrprop(s, 'upper')) && ~any(isstrprop(s, 'lower'));
        end

        function s = escape(s)
            %escape - Escape &, <, >, " and ' for HTML
            s = strrep(s, '&', '&amp;');
            s = strrep(s, '<', '&lt;');
            s = strrep(s, '>', '&gt;');
            s = strrep(s, '"', '&quot;');
            s = strrep(s, '''', '&#x27;');
        end

        function s = stem(file)
            %stem - File name without folder and extension
            [~, s] = fileparts(file);
            s = char(s);
        end

    end
end