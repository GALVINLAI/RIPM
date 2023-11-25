% 该脚本用于生成一个README模板，列出当前文件夹及其子文件夹中的文件和文件夹。
% 生成的README模板以Markdown格式保存在当前文件夹中的README.md文件中。

% 获取当前脚本所在文件夹路径
currentFolder = fileparts(mfilename('fullpath'));

% 生成README.md文件路径
readmeFilePath = fullfile(currentFolder, 'README.md');

% 打开README.md文件以写入模板内容
readmeFile = fopen(readmeFilePath, 'w');

% 写入README标题
fprintf(readmeFile, '# README\n\n');

% 获取当前文件夹及其子文件夹中的文件和文件夹
fileList = dir(currentFolder);

% 遍历文件夹中的内容
for i = 1:numel(fileList)
    item = fileList(i);
    
    % 排除当前文件夹(.)和父文件夹(..)
    if ~strcmp(item.name, '.') && ~strcmp(item.name, '..')
        if item.isdir
            % 如果是文件夹，写入文件夹标题
            fprintf(readmeFile, '## %s\n\n', item.name);
            
            % 获取文件夹中的文件列表
            subfolderPath = fullfile(currentFolder, item.name);
            subfileList = dir(subfolderPath);
            
            % 遍历文件夹中的文件
            for j = 1:numel(subfileList)
                subitem = subfileList(j);
                % 排除子文件夹(.)和父文件夹(..)
                if ~strcmp(subitem.name, '.') && ~strcmp(subitem.name, '..')
                    % 写入文件链接
                    fprintf(readmeFile, '- [%s](%s)\n', subitem.name, fullfile(item.name, subitem.name));
                end
            end
        else
            % 如果是文件，直接写入文件链接
            fprintf(readmeFile, '- [%s](%s)\n', item.name, item.name);
        end
    end
end

% 关闭README.md文件
fclose(readmeFile);

% 显示生成成功消息
disp('README模板生成成功！');
