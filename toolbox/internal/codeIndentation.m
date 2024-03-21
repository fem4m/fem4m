function codeIndentation(curPath,insertCopyright)
%codeIndentation - Re-indent all M-files below a folder
%   codeIndentation() rewrites every M-file in the current folder and
%   its subfolders: it removes the empty lines at the end of the file,
%   leaves exactly one empty line before and after the first line
%   that contains the fem4m copyright text, writes the last line
%   without a line break and smart indents the whole file in the
%   MATLAB Editor. Subfolders whose name contains '.' and the file
%   codeIndentation.m itself are skipped. The names of the folders and
%   files are displayed.
%
%   codeIndentation(curPath) rewrites the M-files below the folder
%   curPath instead. An empty curPath stands for the current folder.
%
%   codeIndentation(curPath,insertCopyright) with insertCopyright true
%   also inserts the copyright line, between empty lines, after the
%   help block (the lines 2, 3, ... that contain '%') of every file
%   without one; a file without help first gets an empty comment line.
%   The default is false.
%
%   See also fem4mPath, matlab.desktop.editor.openDocument
%
%   Note: The files are overwritten in place. The smart indent needs the
%   MATLAB Editor (matlab.desktop.editor).

% Copyright (c) 2026 by the fem4m authors

if nargin < 1 || isempty(curPath)
    curPath = pwd();
end

if nargin < 2 || isempty(insertCopyright)
    insertCopyright = false;
end

folder = dir(curPath);
for  j = 1 : size(folder,1)
    if folder(j).isdir && isempty(strfind(folder(j).name,'.'))
        nextPath = fullfile(curPath,folder(j).name);
        disp(nextPath)
        codeIndentation(nextPath,insertCopyright);
    end
end

mFile = what(curPath);
mFile = mFile.m;
curFolder = split(curPath,filesep);
curFolder = curFolder{end};
for i = 1 : size(mFile,1)
    % remove trailing empty lines
    fileID = fopen(fullfile(curPath,mFile{i}),'rt');
    lines = textscan( fileID, '%s', 'Delimiter','\n', 'CollectOutput',true);
    lines = lines{1};
    while isempty(lines{end})
        lines(end) = [];
    end
    fclose(fileID);
    %  add copyright
    % don't add copyright for any 'Contents.m' or any examples file
    if ~strcmpi(mFile{i},'Contents.m') && ~strcmpi(curFolder,'examples')
        fileID = fopen(fullfile(curPath,mFile{i}),'wt');
        fprintf(fileID, '%s\n', lines{1});
        k = 2;
        if insertCopyright && ~any(contains(lines,'% Copyright (c) 2026 by the fem4m authors'))
            while(contains(lines{k},'%'))
                fprintf(fileID, '%s\n', lines{k});
                k = k+1;
            end
            if k == 2
                fprintf(fileID,'%% \n');
            end
            fprintf(fileID,'\n');
            fprintf(fileID,'%% Copyright (c) 2026 by the fem4m authors\n');
            fprintf(fileID,'\n');
            while k <= numel(lines) && isempty(lines{k})
                k = k + 1;
            end
            for j = k : numel(lines)-1
                fprintf(fileID, '%s\n', lines{j});
            end
        else
            indx = find(contains(lines,'% Copyright (c) 2026 by the fem4m authors'),1);
            if isempty(indx)
                indx = -1;
            else
                j = indx-1;
                while j > 0 && isempty(lines{j})
                    lines(j) = [];
                    j = j - 1;
                    indx = indx - 1;
                end
                j = indx+1;
                while j <= numel(lines) && isempty(lines{j})
                    lines(j) = [];
                end
            end
            for j = k : numel(lines)-1
                if j == indx+1 && ~isempty(lines{j})
                    fprintf(fileID,'\n');
                end
                fprintf(fileID, '%s\n', lines{j});
                if j == indx-1 && ~isempty(lines{j})
                    fprintf(fileID,'\n');
                end
            end
        end
        fprintf(fileID, '%s', lines{end});
        fclose(fileID);
    end
    % indentation
    document = matlab.desktop.editor.openDocument(fullfile(curPath,mFile{i}));
    smartIndentContents(document);
    save(document);
    close(document);
    disp(mFile{i})
end

end