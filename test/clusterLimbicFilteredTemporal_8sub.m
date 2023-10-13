clearvars, clc, close all

%% Plot mean waveforms across subjects 
% Dependencies: matmef and vistasoft github repositories. 
% cd to HAPwave repo
addpath(genpath(pwd)) 

% set local path to your BIDS directory:
myPath          = setLocalDataPath(1);
localDataPath   = myPath.input;

% load the metadata
all_subjects    = {'01','02','03','04','05','06','07','08'};        % List of subjects
all_hemi        = {'r','r','r','l','r','l','l','r'};                % List of hemispheres
all_runs        = {'01','01','01','01','01','01','01','01'};        % List of runs

for ss = 1:length(all_subjects)
    bids_sub        = all_subjects{ss};
    bids_ses        = 'ieeg01';
    bids_task       = 'ccep';
    bids_run        = all_runs{ss};

    % Load metadata and stats
    [events_table,channels_table,electrodes_table,sub_out] = pcc_loadAveragesStats(localDataPath,bids_sub,bids_ses,bids_task,bids_run);
    all_out(ss) = sub_out;
end

% %% Now correct P values for number of comparisons in each subject %%
% Fisrt load the limbic codes in both right and left hemis
area_codes      = {[12123 53 54 12108 12109 12110 12106 12107 59 ...
                    11123 17 18 11108 11109 11110 11106 11107 10]};
nr_subs         = length(all_subjects);

% Loop over subjects
for ss = 1:nr_subs
   
    % List sites that belong to the recording ROI (measured_area)
    these_measured_sites    = find(ismember(all_out(ss).channel_areas,area_codes{1}));
    
    % List sites that belong to the stimulated ROI (stim_area)
    these_stim_sites        = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{1})...
                                 | ismember(all_out(ss).average_ccep_areas(:,2),area_codes{1}));
    
    % prepare for correction of multiple comparisons of p-values
    all_out(ss).hasdata     = NaN(size(all_out(ss).crp_out));
    all_out(ss).crp_p       = NaN(size(all_out(ss).crp_out));
    all_out(ss).a_prime     = NaN(size(all_out(ss).crp_out));
    all_out(ss).cod         = NaN(size(all_out(ss).crp_out));
    all_out(ss).crp_p_adj   = NaN(size(all_out(ss).crp_out));
    all_out(ss).h           = NaN(size(all_out(ss).crp_out));
    all_out(ss).avg_trace_tR = zeros(size(all_out(ss).average_ccep));
    
    % loop over measured sites
    for kk = 1:length(these_measured_sites)
        % loop over the stimulated pairs
        for ll = 1:length(these_stim_sites)
            if ~isempty(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).data)   % we have CRPs
                all_out(ss).hasdata(these_measured_sites(kk), these_stim_sites(ll)) = 1;            % is significant
                all_out(ss).crp_p(these_measured_sites(kk), these_stim_sites(ll))   = all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_projs.p_value_tR;   % p-value
                all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)) = mean(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.al_p); % mean 
                all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll))     = median(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.cod);    % CRP coefficent of determination
                sig_timepoints = find(all_out(ss).tt>0.015 & all_out(ss).tt<=all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.tR);                 % response duration
                all_out(ss).avg_trace_tR(these_measured_sites(kk), these_stim_sites(ll),sig_timepoints) = all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.C; % al_p for tR
            else
                all_out(ss).hasdata(these_measured_sites(kk), these_stim_sites(ll)) = 0;            % is not significant
            end
        end
    end
    pvals                   = all_out(ss).crp_p(all_out(ss).hasdata==1);  % get pVals
    qq                      = 0.05; % false discovery rate
    [h, crit_p, adj_ci_cvrg, adj_p]                 = fdr_bh(pvals,qq,'dep','no'); % Benjamini & Yekutieli FDR correction 
    all_out(ss).crp_p_adj(all_out(ss).hasdata==1)   = adj_p; % Adjusted pVals
    all_out(ss).h(all_out(ss).hasdata==1)           = h; % adjusted pVal is significant
end

% %% Load stim and measure sites %%
% Sort limbic codes by hemisphere
area_names      = {'Hipp','Amyg','PCC','ACC','ANT'};   
area_codes_r    = {[12123 53],[54],[12108 12109 12110],[12106 12107],[59]}; % right
area_codes_l    = {[11123 17],[18],[11108 11109 11110],[11106 11107],[10]}; % left

out             = [];                   % prepare areaByarea structure, with all subjects concatinated for each area
subj_resp_total = zeros(nr_subs,1);     % stim-->measured pair for adjusted FDR
t_win_norm      = [0.015 0.500];        % window for vector length normalization and plotting across subjects

for measure_ind = 1:length(area_names)  % loop through measured sites
    for stim_ind = 1:length(area_names) % now go through stimulated sites
        resp_counter = 0;               % counting all responses across subjects for this connection

        for ss = 1:nr_subs % loop over subjects
            tt = all_out(ss).tt;
            % which hemisphere has coverage
            if isequal(all_hemi{ss},'l')
                area_codes = area_codes_l;
            elseif isequal(all_hemi{ss},'r')
                area_codes = area_codes_r;
            end

            % Get recording ROI (measured_area)
            these_measured_sites    = find(ismember(all_out(ss).channel_areas,area_codes{measure_ind}));
           
            % Get stimulated ROI (stim_area)
            these_stim_sites        = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{stim_ind}) | ...
                                           ismember(all_out(ss).average_ccep_areas(:,2),area_codes{stim_ind}));
            
            % loop over measured sites
            for kk = 1:length(these_measured_sites)
                % loop over the stimulated pairs
                for ll = 1:length(these_stim_sites)
                    if ~isempty(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).data) % ~ same stim/recording site
                        
                        % this is a pair, counting for multiple comparison correction per subject
                        subj_resp_total(ss)     = subj_resp_total(ss) + 1; % set counter
                        
                        % first raw responses
                        plot_responses          = squeeze(all_out(ss).average_ccep(these_measured_sites(kk), these_stim_sites(ll), :));
                        % save outputs
                        resp_counter            = resp_counter + 1;

                        % get CCEP responses for plotting
                        % Scaling to unit length (Euclidean lenght): https://en.wikipedia.org/wiki/Feature_scaling
                        % unit length taken in same window as stats
                        response_vector_length  = sum(plot_responses(all_out(ss).tt > t_win_norm(1) &  all_out(ss).tt < t_win_norm(2)) .^ 2) .^ .5;
                        plot_responses_norm     = plot_responses ./ (response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial
                        out(measure_ind,stim_ind).plot_responses_norm(resp_counter, :) = plot_responses_norm;
                        
                        % store subject index
                        out(measure_ind,stim_ind).subj_ind(resp_counter, :) = ss;
                        
                        % save CRP parms and other params
                        out(measure_ind,stim_ind).elec_relDist(resp_counter, :)= all_out(ss).elec_relDist(these_measured_sites(kk));
                        out(measure_ind,stim_ind).p(resp_counter, :)        = all_out(ss).crp_p_adj(these_measured_sites(kk), these_stim_sites(ll));
                        out(measure_ind,stim_ind).cod(resp_counter, :)      = all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll)); 
                        out(measure_ind,stim_ind).a_prime(resp_counter, :)  = all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)); 
                    end
                end
            end
        end
    end 
end

%% Load matrix and convert to previously used format. This is in response to reviewer 1.11, testing easier forms of processing
srate           = all_out(ss).srate;
ttOrig          = (0:10239)/srate - 2;

V_allsubs       = [];
labels_allsubs  = cell(0, 2); % measured, stim
subNum          = [];

% ordered labels for the rows and columns
areas           = {'HC', 'Amg', 'PCC', 'ACC', 'ANT'};

for ii = 1:size(out, 1)
    for jj = 1:size(out, 2)
        p = out(ii, jj).p; % FDR-adjusted p value
        h = p < 0.05; % significance
        
        resp = out(ii, jj).plot_responses_norm(h, :);
        
        % invert reponses farther than 2.3 mm (PCC measurements only)
        dists = out(ii, jj).elec_relDist(h);
        if ii == 3
            resp(dists > 2.3, :) = -resp(dists > 2.3, :);
        end
        
        subNum          = [subNum; out(ii, jj).subj_ind(h)]; % subject ID
        V_allsubs       = [V_allsubs; resp]; % signals
        labels_allsubs  = [labels_allsubs; repmat({areas{ii}, areas{jj}}, sum(h), 1)];
    end
end

% sort everything by subject
[~, ord]        = sort(subNum);
subNum          = subNum(ord);
V_allsubs       = V_allsubs(ord, :);
labels_allsubs  = labels_allsubs(ord, :);

labels_length   = arrayfun(@(s) sum(subNum == s), 1:8);

% total number of cceps
n               = size(V_allsubs, 1);
assert(n == sum(labels_length) && n == length(labels_allsubs));

% number of subjects
nsubs           = length(labels_length);

%% Normalize, Plot all average CCEPs by type
win         = [0.1, 0.5];

% normalize CCEPs to be comparable cross subjects
V_mag       = vecnorm(V_allsubs(:, ttOrig >= win(1) & ttOrig < win(2)), 2, 2); % magnitude on 100-1000 ms segment, for normalization
V_norm      = V_allsubs ./ V_mag;

% remove stim artifact to prevent smearing during filtering
V_norm      = V_norm(:, ttOrig >= 0.015);
ttNoArt     = ttOrig(ttOrig >= 0.015);

% memory = hipp->PCC, emotional = amyg->ACC
pathIdxes           = struct;
pathIdxes.name      = {'amygAcc', 'amygPcc', 'antPcc', 'hippAcc', 'hippPcc'};
pathIdxes.colors    = [1, 0.3, 0; 0, 1, 0.1; 0, 1, 1; 1, 1, 0; 0.2, 0.2, 1];
pathIdxes.idxes     = getIdxes(labels_allsubs);
           
% Plot all CCEPs by type
figure('Position', [100, 100, 1200, 800]); t = tiledlayout(3, 2);
t.TileSpacing       = 'compact';
t.Padding           = 'compact';
for ii = 1:5
    nexttile;
    plot(ttNoArt, V_norm(pathIdxes.idxes{ii}, :)', 'Color', 0.5*[1, 1, 1]);
    hold on
    plot(ttNoArt, mean(V_norm(pathIdxes.idxes{ii}, :)), 'k-', 'LineWidth', 1.5);
    xlim([0, win(2)]);
    title(pathIdxes.name{ii});
end

%% Low-pass or bandpass filter the signal
% half-power frequency for low-pass or bandpass.
% no filtering = use nan. Low-pass, use 1 value, band-pass, use 2 values (e.g., [8, 13]
% lowpass: 4, 8, and 13
% theta: 4-8, alpha 8-13, beta 13-30
fhps = [{nan}, {4}, {8}, {13}, {[4, 8]}, {[8, 13]}, {[13, 30]}];

for ff = 1:length(fhps)
    
    fhp = fhps{ff};
    
    if any(isnan(fhp)) % don't do any filtering
        fprintf('no filtering done\n');
        V_filt = V_norm;
        outdir = fullfile('output', 'nofilter');
        mkdir(outdir);
    elseif length(fhp) == 1 % lowpass
        d = designfilt('lowpassiir', 'FilterOrder', 4, 'DesignMethod', 'Butter', ...
                   'HalfPowerFrequency', fhp, 'SampleRate', srate);
        outdir = fullfile('output', sprintf('lowpass_%d', fhp));
        mkdir(outdir);
        fprintf('lowpass at F=%d Hz\n', fhp);
        V_filt = filtfilt(d, V_norm')';
    elseif length(fhp) == 2 % bandpass           
       d = designfilt('bandpassiir', 'FilterOrder', 4, ...
                    'DesignMethod', 'butter', ...
                    'HalfPowerFrequency1', fhp(1), ... 
                    'HalfPowerFrequency2', fhp(2), ...
                    'SampleRate', srate);
        outdir = fullfile('output', sprintf('bandpass_%d_%d', fhp(1), fhp(2)));
        mkdir(outdir);
        fprintf('bandpass between F=%d and %d Hz\n', fhp(1), fhp(2));
        V_filt = filtfilt(d, V_norm')';
    end

    % clip after filtering to 100 - 1000 ms to match the DWT data. The raw signals were normalized on this interval (previous section)
    V_filt      = V_filt(:, ttNoArt >= win(1) & ttNoArt < win(2));
    tt          = ttNoArt(ttNoArt >= win(1) & ttNoArt < win(2));

    figure('Position', [100, 100, 900, 600]); t = tiledlayout(3, 2);
    t.TileSpacing = 'compact';
    t.Padding   = 'compact';
    for ii = 1:5
        nexttile;
        plot(tt, V_filt(pathIdxes.idxes{ii}, :)', 'Color', 0.5*[1, 1, 1]);
        hold on
        plot(tt, mean(V_filt(pathIdxes.idxes{ii}, :)), 'k-', 'LineWidth', 1.5);
        title(pathIdxes.name{ii});
    end
    saveas(gcf, fullfile(outdir, 'V_filt'), 'png');


    %% Calculate Principal Components of all categories, using low-passed signals, to determine which PCs to use and variance explained

    % SVD on centered data (centered across time points)
    [U, S, V]   = svd((V_filt - mean(V_filt))', 'econ');
    score       = V*S; % weigh right singular vectors by S to get the PCA scores. Rows are observations, columns correspond to PCs. Note that the data can be reconstructed as score*U'
    varExp      = diag(S^2); 
    varExp      = 100*varExp/sum(varExp);

    % save variance explained to text
    fid         = fopen(fullfile(outdir, 'varExp.txt'), 'w');
    fprintf(fid, '%0.02f\n', varExp(1:20));
    fclose(fid);

    % find top 2 PCs to use by tstat
    tstats = zeros(5, 1);
    for ii = 1:5
        dist1 = score(pathIdxes.idxes{1}, ii);
        dist2 = score(pathIdxes.idxes{5}, ii);
        [~, ~, ~, stats] = ttest2(dist1, dist2);
        tstats(ii) = abs(stats.tstat);
    end
    [~, order]  = sort(tstats, 'descend');
    PCs2Use     = sort(order(1:2)); % PCs to plot and use for LDA
    fprintf('Using PC %d and PC %d, total varExp = %0.2f\n', PCs2Use(:), sum(varExp(PCs2Use)));
    
    % Calculate training accuracy
    training    = [score(pathIdxes.idxes{1}, PCs2Use); score(pathIdxes.idxes{5}, PCs2Use)]; % amygACC first, hippPCC second
    group       = [ones(length(pathIdxes.idxes{1}), 1); 5*ones(length(pathIdxes.idxes{5}), 1)]; % CCEP type to use for classification\
    [class, ~, ~, ~, coeff] = classify(training, training, group); % find decision boundary coefficients for model fit on all data (without leaving any out)
    accTrain    = 100 * sum(class == group) / length(class); % training accuracy
    fprintf('Training accuracy: %0.2f%%\n', accTrain);


    %% Leave-one-subject out CV analysis to calculate accuracy

    accSubs     = nan(nsubs, 1); % CV test accuracy for each subject withheld

    for sub = 1:nsubs

        % divide into test and train sets for DWT CCEPs and labels
        V_filtTest      = V_filt(subNum == sub, :);
        V_filtTrain     = V_filt(subNum ~= sub, :);
        labels_test     = labels_allsubs(subNum == sub, :);
        labels_train    = labels_allsubs(subNum ~= sub, :);

        % PCA on the training set
        [U, S, V]       = svd((V_filtTrain - mean(V_filtTrain))', 'econ');
        scoreTrain      = V*S;

        scoreTest       = (V_filtTest - mean(V_filtTest))*U; % scores of testing set (mean centered projection onto loadings)

        % scores and groupings for training data with amygACC (first) and hippPCC (second)
        idxesTrain      = getIdxes(labels_train);
        training        = [scoreTrain(idxesTrain{1}, PCs2Use); scoreTrain(idxesTrain{5}, PCs2Use)];
        groupTrain      = [ones(length(idxesTrain{1}), 1); 5*ones(length(idxesTrain{5}), 1)];

        % scores and groupings for test data
        idxesTest       = getIdxes(labels_test);
        testing         = [scoreTest(idxesTest{1}, PCs2Use); scoreTest(idxesTest{5}, PCs2Use)];
        groupTest       = [ones(length(idxesTest{1}), 1); 5*ones(length(idxesTest{5}), 1)];

        class           = classify(testing, training, groupTrain); % find decision boundary coefficients for model fit on all data (without leaving any out)
        accSubs(sub)    = 100*sum(class == groupTest) / length(class); % training accuracy

    end

    fprintf('Mean leave-one-subject-out CV accuracy: %0.2f%%\n', mean(accSubs));

    % Real CV accuracy overlaid with null accuracy +/- SD
    figure; hold on
    bar(1:nsubs, accSubs, 'FaceColor', [0.8, 0.8, 0.8]);
    xlim([0, nsubs+1]); ylim([0, 110]);
    set(gca, 'ytick', 0:20:100);
    xlabel('Cross-validation trial'); ylabel('Accuracy (%)');
%     saveas(gcf, fullfile(outdir, sprintf('CVaccuracy_PCs%d-%d', PCs2Use(1), PCs2Use(2))), 'png');
%     saveas(gcf, fullfile(outdir, sprintf('CVaccuracy_PCs%d-%d', PCs2Use(1), PCs2Use(2))), 'svg');

%     % save the individual CV accuracies to file
%     fid = fopen(fullfile(outdir, sprintf('CVaccuracy_PCs%d-%d.txt', PCs2Use(1), PCs2Use(2))), 'w');
%     fprintf(fid, 'TrainAcc\t%0.02f\n', accTrain);
%     fprintf(fid, 'Accuracy\t%0.02f\n', accSubs);
%     fclose(fid);

end


%% Load CV accuracy at each low-pass threshold, plot jitter plots

rng('default'); % for the jitterplot

nsubs       = 8;

dirs        = {'.', 'nofilter', 'lowpass_4', 'lowpass_8', 'lowpass_13', 'bandpass_4_8', 'bandpass_8_13', 'bandpass_13_30'}; % '.' is the current method
cvaccs      = zeros(nsubs, length(dirs));
trainaccs   = zeros(1, length(dirs));

for ii = 1:length(dirs)
    accFile         = glob(fullfile('output', dirs{ii}, 'CVaccurac*.txt'));
    T = readtable(accFile{1}, 'FileType', 'text', 'Delimiter', '\t');
    trainaccs(ii)   = T.Var2(1); % training accuracy
    cvaccs(:, ii)   = T.Var2(2:nsubs+1); % CV accuracy
end

% put in spaces between each type of preprocessing
trainaccs   = [trainaccs(1), 0, trainaccs(2), 0, trainaccs(3:5), 0, trainaccs(6:8)];
cvaccsCell  = num2cell(cvaccs, 1);
cvaccsCell  = [cvaccsCell(1), {[]}, cvaccsCell(2), {[]}, cvaccsCell(3:5), {[]}, cvaccsCell(6:8)];
cvaccsMean  = mean(cvaccs);
cvaccsMean  = [cvaccsMean(1), 0, cvaccsMean(2), 0, cvaccsMean(3:5), 0, cvaccsMean(6:8)];

figure('Position', [200, 200, 400, 200]); hold on
bar(cvaccsMean, 'FaceColor', [0.8, 0.8, 0.8]);
hs = jitterplot(cvaccsCell, 'b');
for ii = 1:length(hs)
    hs(ii).Marker       = '.';
    hs(ii).MarkerSize   = 10;
end
ylim([0, 110]);
set(gca, 'ytick', 0:20:100);
set(gca, 'xtick', []);
saveas(gcf, fullfile('output', 'cvaccs_alltypes'), 'png');
saveas(gcf, fullfile('output', 'cvaccs_alltypes'), 'svg');


%% Functions

% {'amygAcc', 'amygPcc', 'antPcc', 'hippAcc', 'hippPcc'} in that order
function idxes = getIdxes(labels)
    idxes = {find(strcmp(labels(:, 2), 'Amg') & strcmp(labels(:, 1), 'ACC')); % Amyg -> ACC, to-cluster
             find(strcmp(labels(:, 2), 'Amg') & strcmp(labels(:, 1), 'PCC'));
             find(strcmp(labels(:, 2), 'ANT') & strcmp(labels(:, 1), 'PCC'));
             find(strcmp(labels(:, 2), 'HC') & strcmp(labels(:, 1), 'ACC'));
             find(strcmp(labels(:, 2), 'HC') & strcmp(labels(:, 1), 'PCC'))}; % hipp -> PCC, to cluster
end

function idxes = getIdxesOld(labels)
% Returns the indices in labels that correspond to each of the 5 CCEP connections

    idxes = {find(labels(:, 2)==2 & ismember(labels(:, 1), [3, 4])); % Amyg -> ACC, to-cluster
             find(labels(:, 2)==2 & ismember(labels(:, 1), [5, 6, 7]));
             find(labels(:, 2)==9 & ismember(labels(:, 1), [5, 6, 7]));
             find(ismember(labels(:, 2), [1, 8]) & ismember(labels(:, 1), [3, 4]));
             find(ismember(labels(:, 2), [1, 8]) & ismember(labels(:, 1), [5, 6, 7]))}; % hipp -> PCC, to cluster
end
