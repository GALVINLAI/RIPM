function Boss_2_Stiefel_NonnegativeProjection

%-------------------------- Projection onto Nonnegative Subset of Stiefel Manifold ----
%-------------------------- Measurement of Robustness ---------------------

%%
specifier.matlabversion = 1; % 0 if older than 2015 1 otherwise
% method set = {'RALM','REPM(LQH)','REPM(LSE)','RSQP','RIPM'};
% RIPM uses KrylovIterMethod. 
specifier.ind = [1,1,1,1,1];  % [1,1,1,1,1]

% size of problem
rdimset = [40 50 60 70]; %%%%%%%%%%%%% [40 50 60 70]

% KKT residual
tolKKTrespowerset = [6]; %%%%%%%%%%%%%% log10 scale, i.e., 1e-* tolerance

% Number of repeat on same set of data
n_repeat = 20;  %%%%%%%%%%%%%% 20;

maxtime = 600; %%%%%%%%%%%%%

AllResultBox = zeros(20, n_repeat+1);

for n = rdimset
    k = 0.2*n; % Automatically decided.
    % For rdim in [40 50 60 70], then
    %     cdim in [8  10 12 14]

    for tolKKTres = tolKKTrespowerset

        AllResult = AllResultBox;

        for repeat = 1 : n_repeat
            %_______Set up data______

            %______Set Object matrix_____
            B = rand_Nonneg_St(n,k); % random B
            X1 = (B>0).*(1+rand(n,k)); % Generate the same style as B
            Xstar = X1./sqrt(sum(X1.*X1));% normalize every columns; now, Xstar is Nonnegative_Stiefel
            L = rand(k,k);
            L = L + k*eye(k);
            C = Xstar*L';

            %________Experiment_____
            % common options
            options.tolKKTres = 10^(-tolKKTres); % tolKKTres = 8
            options.maxOuterIter = 1000; % for RALM, REPMs.
            options.maxiter = 10000;  % for RSQP, RIPM.
            options.maxtime = maxtime;  %%%%%% 60
            options.verbosity = 1;  % 1
            % for RIPM 
            options.KrylovIterMethod = 1; %%%%%%%%%%%%%
            % for RALM, REPMs.
            options.startingtolgradnorm = max(1e-3, 10^(-tolKKTres + 3));
            options.endingtolgradnorm = 10^(-6);
            options.outerverbosity = options.verbosity;
            % for RSQP
            options.mineigval_correction = 1e-5;  % 1e-5

            %________Setting________
            setting.repeat = repeat;
            setting.row_dim = n;
            setting.col_dim = k;

            setting.tolKKTres =  options.tolKKTres;
            setting.maxOuterIter = options.maxOuterIter;
            setting.maxiter = options.maxiter;
            setting.maxtime = options.maxtime;
            setting.verbosity = options.verbosity;

            setting.Xstar = Xstar;
            setting.C = C;

            setting.ExperimentName = mfilename();
            setting.SimpleName = setting.ExperimentName(8:end);
            setting.filepath = sprintf('./numerical results/Mod_St/nrep%d_Row%d_Col%d_KKTtol%.1e',...
                setting.repeat,setting.row_dim,setting.col_dim,setting.tolKKTres);

            result = client_Stiefel_NonnegativeProjection(n, k, Xstar, C, options, specifier, setting);

            result = result(:);
            AllResult(:,repeat) = result; % Add one column of results per repeat.
        end

        AllResult = Statistics(AllResult); % Add a column of statistics to the last column for repeats.
        filename = sprintf('./numerical results/Mod_St/RC_%s_%s.csv',setting.SimpleName,setting.filepath);
        writematrix(AllResult,filename)

        table = reshape(AllResult(:,end),[4,5]);
        table = table'; % data(:,i) = [residual; time; iternum; NLRMrelres];
        filename = sprintf('./numerical results/Mod_St/Table_%s_%s.csv',setting.SimpleName,setting.filepath);
        writematrix(table,filename)

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
            OtherMean = mean(AllResult(row+3,success_ind));

            AllResult(row,n_repeat+1)=SuccessRate;
            AllResult(row+1,n_repeat+1)=TimeMean;
            AllResult(row+2,n_repeat+1)=IternumMean;
            AllResult(row+3,n_repeat+1)=OtherMean;
        end

    end


end


