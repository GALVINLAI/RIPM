function Boss_1_fixedrank_NLRM

%-------------------------- Nonnegative low rank matrix approximation (NLRM) ----
%-------------------------- Measurement of Robustness ---------------------

%%
specifier.matlabversion = 1; % 0 if older than 2015 1 otherwise
% method set = {'RALM','REPM(LQH)','REPM(LSE)','RSQP','RIPM'};
% RIPM uses KrylovIterMethod. 
specifier.ind = [1,1,1,1,1];  % [1,1,1,1,1]

% sd is standard deviation for Gaussian Noise.
sdset = [0, 0.001, 0.01]; %%%%%%%%%%%%% [0, 0.001, 0.01]

% size of problem
rdimset = [20, 30, 40]; %%%%%%%%%%%%% [20, 30, 40]

% KKT residual
tolKKTrespowerset = [8]; %%%%%%%%%%%%%% log10 scale, i.e., 1e-* tolerance

% Number of repeat on same set of data
n_repeat = 20;  %%%%%%%%%%%%% 20;

maxtime = 180; %%%%%%%%%%%%%

AllResultBox = zeros(20, n_repeat+1);

for sd = sdset
    for rdim = rdimset
        cdim = 0.8*rdim; % Automatically decided.
        rankval = max(2, round(0.1*rdim)); % Automatically decided.
        % For rdim in [20, 30, 40], then
        %     cdim in [16, 24, 32],
        %  rankval in [ 2,  3,  4].

        for tolKKTres = tolKKTrespowerset

            AllResult = AllResultBox;

            for repeat = 1 : n_repeat
                %_______Set up data______

                %______Set Object matrix_____
                % Generate a random rdim x cdim matrix origanl A of rank r
                L = rand(rdim, rankval);
                R = rand(rankval, cdim);
                origanl_A =  L*R;
                % Gaussian Noise
                mu = 0; % mean
                vr = sd.^2; % vr is variance.
                GaussianNoise = mu + sqrt(vr)*randn(rdim,cdim);
                A = origanl_A + GaussianNoise; % data matrix A becomes full rank.

                %________Experiment_____
                % common options
                options.tolKKTres = 10^(-tolKKTres); % tolKKTres = 8
                options.maxOuterIter = 1000; % for RALM, REPMs.
                options.maxiter = 10000;  % for RSQP, RIPM.
                options.maxtime = maxtime;  % 60
                options.verbosity = 1;  % 1
                options.rank = rankval;
                % for RIPM
                options.KrylovIterMethod = 1;
                % for RALM, REPMs.
                options.startingtolgradnorm = max(1e-3, 10^(-tolKKTres + 3)); % 1e-3
                options.endingtolgradnorm = 10^(-6);
                options.outerverbosity = options.verbosity;
                % for RSQP
                options.mineigval_correction = 1e-5;  % 1e-5

                %________for initial point_____
                setting.initialpoint =  'random';

                %________Setting________
                setting.repeat = repeat;
                setting.row_dim = rdim;
                setting.col_dim = cdim;
                setting.rank = rankval;

                setting.tolKKTres =  options.tolKKTres;
                setting.maxOuterIter = options.maxOuterIter;
                setting.maxiter = options.maxiter;
                setting.maxtime = options.maxtime;
                setting.verbosity = options.verbosity;

                setting.A = A;
                setting.sd = sd;

                setting.ExperimentName = mfilename();
                setting.SimpleName = setting.ExperimentName(8:end);
                setting.filepath = sprintf('nrep%d_sd%g_Row%d_Col%d_Rank%d_KKTtol%.1e',...
                    setting.repeat,setting.sd,setting.row_dim,setting.col_dim,setting.rank,setting.tolKKTres);

                result = client_fixedrank_NLRM(rdim, cdim, rankval, A, options, specifier, setting);

                result = result(:);
                AllResult(:,repeat) = result; % Add one column of results per repeat.
            end

            AllResult = Statistics(AllResult); % Add a column of statistics to the last column for repeats.
            filename = sprintf('RC_%s_%s.csv',setting.SimpleName,setting.filepath);
            writematrix(AllResult,filename)

            table = reshape(AllResult(:,end),[4,5]);
            table = table'; % data(:,i) = [residual; time; iternum; NLRMrelres];
            filename = sprintf('Table_%s_%s.csv',setting.SimpleName,setting.filepath);
            writematrix(table,filename)

        end
    end
end

%% sub function

    function AllResult = Statistics(AllResult)

        for method_ind=1:5 % data(:,i) = [residual; time; iternum; NLRMrelres];
            row = (method_ind-1)*4+1;

            all_residual = AllResult(row,1:n_repeat);
            success_ind = find(all_residual < options.tolKKTres);

            SuccessRate = length(success_ind)/n_repeat;
            TimeMean = mean(AllResult(row+1,success_ind));
            IternumMean = mean(AllResult(row+2,success_ind));
            NLRMrelresMean = mean(AllResult(row+3,success_ind));

            AllResult(row,n_repeat+1)=SuccessRate;
            AllResult(row+1,n_repeat+1)=TimeMean;
            AllResult(row+2,n_repeat+1)=IternumMean;
            AllResult(row+3,n_repeat+1)=NLRMrelresMean;
        end

    end


end


