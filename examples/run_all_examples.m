% ONLY FOR FUNCTION FILES!!

clear all; %#ok<CLALL>
close all;
clc;

List = dir('*.m');

issues = {};

for k = 1 : numel(List)
    s = List(k).name;
    s = s(1:end-2); % remove '.m'
    if ~strcmp(s, mfilename()) % don't run yourself
        try
            eval([s, ';']);
        catch
            issues{end+1} = s; %#ok<SAGROW>
        end
    end
end

if ~isempty(issues)
    warning('There were issues running the following scripts:');
    disp(issues);
end

close all;

