
%%为了实验NLRM的制表脚本,私用不公开.

%%
n_repeat = 20; %%%%%%%%%%%%%
tolKKTres = 1e-8; %%%%%%%%%%%%%

%% [手动确认顺序] 提前切换当前文件夹; 更具问题不同要调参
List = dir('RC_fixedrank_NLRM_nrep20*');
RC_filename_list = {List.name}'; 
new_RC_filename_list = RC_filename_list(end-2:end);
new_RC_filename_list(4:9) = RC_filename_list(1:6);

%% 制表TotalTable是结果
Table = [];
subTable = [];
for k = 1: numel(new_RC_filename_list)
    RC_filename = new_RC_filename_list{k};
    AllResult = readmatrix(RC_filename);
    subTable = Statistics(AllResult, n_repeat, tolKKTres);
    Table = [Table, subTable];
end

TotalTable = NaN(15,9);
TotalTable(1:5,:)=Table(:,1:9);
TotalTable(6:10,:)=Table(:,10:18);
TotalTable(11:15,:)=Table(:,19:27);
%TotalTable(16:20,:)=Table(:,28:36);

function subTable = Statistics(AllResult, n_repeat, tolKKTres)

subTable = NaN(5,3);

for method_ind=1:5 % data(:,i) = [residual; time; iternum; NLRMrelres];
    row = (method_ind-1)*4+1;

    all_residual = AllResult(row,1:n_repeat);
    success_ind = find(all_residual < tolKKTres);

    SuccessRate = length(success_ind)/n_repeat;
    TimeMean = mean(AllResult(row+1,success_ind));
    IternumMean = mean(AllResult(row+2,success_ind));
    % NLRMrelresMean = mean(AllResult(row+3,success_ind));

    subTable(method_ind,:) = [SuccessRate, TimeMean, IternumMean];
end

end
