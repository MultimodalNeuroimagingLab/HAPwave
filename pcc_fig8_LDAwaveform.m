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
    all_out(ss)     = sub_out;
end

% Correct limbic P values for number of comparisons in each subject
% Fisrt set the limbic codes for both hemispheres (Destrieux atlas)
area_codes      = {[12123 53 54 12108 12109 12110 12106 12107 49 11123 17 18 11108 11109 11110 11106 11107 10]};
nr_subs         = length(all_subjects);

% Loop over subjects to get CCEP params of limbic CCEPs
for ss = 1:nr_subs
   
    % Lists recording ROI (measured_area)
    these_measured_sites    = find(ismember(all_out(ss).channel_areas,area_codes{1}));
    
    % Lists stimulated ROI (stim_area)
    these_stim_sites        = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{1})...
                                 | ismember(all_out(ss).average_ccep_areas(:,2),area_codes{1}));
    
    % Prepare variables for correction of multiple comparisons of p-values on CRPs
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
            if ~isempty(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).data)
                all_out(ss).hasdata(these_measured_sites(kk), these_stim_sites(ll)) = 1;            % mark if CRP calculated 
                all_out(ss).crp_p(these_measured_sites(kk), these_stim_sites(ll))   = all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_projs.p_value_tR;     % get the p-value
                all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)) = mean(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.al_p);     % get the alpha coefficient weights
                all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll))     = median(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.cod);    % get the coefficent of determination
                sig_timepoints = find(all_out(ss).tt>0.015 & all_out(ss).tt<=all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.tR);                     % get the response duration
                all_out(ss).avg_trace_tR(these_measured_sites(kk), these_stim_sites(ll),sig_timepoints) = all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.C; % get canonical shape
            else
                all_out(ss).hasdata(these_measured_sites(kk), these_stim_sites(ll)) = 0;            % no CRPs
            end
        end
    end
    % Correction of multiple comparisons for this subject
    pvals                                           = all_out(ss).crp_p(all_out(ss).hasdata==1);
    qq                                              = 0.05;         % false discovery rate
    [h, crit_p, adj_ci_cvrg, adj_p]                 = fdr_bh(pvals,qq,'dep','no'); % Benjamini & Yekutieli FDR correction 
    all_out(ss).crp_p_adj(all_out(ss).hasdata==1)   = adj_p;        % Adjusted pVals
    all_out(ss).h(all_out(ss).hasdata==1)           = h;            % adjusted pVal is significant
end

%% Load params for analysis
% Sort limbic codes by hemisphere
area_names      = {'Hipp','Amyg','PCC','ACC','ANT'};   
area_codes_r    = {[12123 53],[54],[12108 12109 12110],[12106 12107],[49]}; % right
area_codes_l    = {[11123 17],[18],[11108 11109 11110],[11106 11107],[10]}; % left

out             = [];                   % prepare areaByarea structure, with all subjects concatinated for each area
subj_resp_total = zeros(nr_subs,1);     % stim-->measured pair for adjusted FDR
t_win_norm      = [0.015 0.500];        % window for vector length normalization and plotting across subjects

for measure_ind = 1:length(area_names)  % loop over measured sites
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
                        
                        % save parms
                        out(measure_ind,stim_ind).elec_relDist(resp_counter, :)= all_out(ss).elec_relDist(these_measured_sites(kk));
                        out(measure_ind,stim_ind).p(resp_counter, :)        = all_out(ss).crp_p_adj(these_measured_sites(kk), these_stim_sites(ll));
                        out(measure_ind,stim_ind).h(resp_counter, :)        = all_out(ss).h(these_measured_sites(kk), these_stim_sites(ll));
                        out(measure_ind,stim_ind).cod(resp_counter, :)      = all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll)); 
                        out(measure_ind,stim_ind).a_prime(resp_counter, :)  = all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)); 
                    end
                end
            end
        end
    end 
end

%% Load matrix and format for analysis
srate           = all_out(ss).srate;            % get sampling rate from metadata
ttOrig          = (0:10239)/srate - 2;          % set start time

V_allsubs       = [];                           % prepare matrix
labels_allsubs  = cell(0, 2);                   % get labels: measured, stim
subNum          = [];                           % prepare subject labels

% ordered labels for the rows and columns
areas           = {'HC', 'Amg', 'PCC', 'ACC', 'ANT'};

for ii = 1:size(out, 1) % stim sites
    for jj  = 1:size(out, 2) % measured sites
        p               = out(ii, jj).p;    % FDR-adjusted p value
        h               = p < 0.05;
%         h               = logical(out(ii,jj).h);     % p < 0.05 
        
        resp            = out(ii, jj).plot_responses_norm(h, :);
        
        % invert reponses farther than 2.3 mm
        dists           = out(ii, jj).elec_relDist(h);
        if ii == 3
            resp(dists > 2.3, :) = -resp(dists > 2.3, :);
        end
        
        subNum          = [subNum; out(ii, jj).subj_ind(h)];    % subject ID
        V_allsubs       = [V_allsubs; resp];                    % signals
        labels_allsubs  = [labels_allsubs; repmat({areas{ii}, areas{jj}}, sum(h), 1)];
    end
end

% sort everything by subject
[~, ord]        = sort(subNum);
subNum          = subNum(ord);
V_allsubs       = V_allsubs(ord, :);

labels_allsubs  = labels_allsubs(ord, :);
labels_length   = arrayfun(@(s) sum(subNum == s), 1:8);

% Number of CCEPs
n               = size(V_allsubs, 1);
assert(n == sum(labels_length) && n == length(labels_allsubs));

% number of subjects
nsubs           = length(labels_length);

%% Some more data formatting and normalization
% inversion of CCEPs
segWin      = [0.1, 0.5];           % set time window

% take segment from 100 ms to 500 ms
V_seg       = V_allsubs(:, ttOrig >= segWin(1) & ttOrig < segWin(2));
tt          = ttOrig(ttOrig >= segWin(1) & ttOrig < segWin(2));

% Normalize CCEPs
V_seg       = V_seg./vecnorm(V_seg, 2, 2);

% Determine category of each CCEP
% memory = hipp->PCC, emotional = amyg->ACC
pathIdxes           = struct;
pathIdxes.name      = {'amygAcc', 'amygPcc', 'antPcc', 'hippAcc', 'hippPcc'};
pathIdxes.colors    = [1, 0.3, 0; 0, 1, 0.1; 0, 1, 1; 1, 1, 0; 0.2, 0.2, 1];
pathIdxes.idxes     = getIdxes(labels_allsubs);

% Plot all CCEPs by type
figure('Position', [100, 100, 900, 600]); t = tiledlayout(3, 2);
t.TileSpacing       = 'compact';
t.Padding           = 'compact';
for ii = 1:5
    nexttile;
    plot(tt, V_seg(pathIdxes.idxes{ii}, :)', 'Color', 0.5*[1, 1, 1]);
    hold on
    plot(tt, mean(V_seg(pathIdxes.idxes{ii}, :)), 'k-', 'LineWidth', 1.5);
    xlim([0, segWin(2)]);
    title(pathIdxes.name{ii});
end

%% Discrete Wavelet Transform on all CCEPs (DWT)
% References:
% https://ieeexplore.ieee.org/document/1625611
% http://www.mayagupta.org/publications/waveletPCA.pdf
% original source: Reference: I. Daubechies, Ten lectures on wavelets, CBMS, SIAM, 61, 1994, 198-202 and 254-256.

wavType     = 'sym4';
levMax      = wmaxlev(length(tt), wavType); % maximum level L of decomposition

V_wav       = [];
V_recon     = [];
for ii = 1:n
    [wav, l]        = wavedec(V_seg(ii, :), levMax, wavType);
    
    % threshold to denoise
    thresh          = prctile(abs(wav), 95);
    wav(abs(wav) < thresh) = 0;
    
    V_wav(ii, :)    = wav;
    V_recon(ii, :)  = waverec(wav, l, wavType); % reconstruct from clean wavelets
end

% Plot reconstructed signals
figure('Position', [100, 100, 900, 600]); t = tiledlayout(3, 2);
t.TileSpacing       = 'compact';
t.Padding           = 'compact';
for ii = 1:5
    nexttile;
    plot(V_recon(pathIdxes.idxes{ii}, :)', 'Color', 0.5*[1, 1, 1]);
    hold on
    plot(mean(V_recon(pathIdxes.idxes{ii}, :)), 'k-', 'LineWidth', 1.5);
    title(pathIdxes.name{ii});    
end

%% Plot example of DWT + threshold step for panel A)

idx         = 51; % choose a CCEP trial 

figure('Position', [200, 600, 400, 150]);
plot(tt, V_seg(pathIdxes.idxes{1}(idx), :)', 'k-', 'LineWidth', 1); ylim([-0.1, 0.1])
yline(0, 'Color', [0.5, 0.5, 0.5]);

[wav, l]    = wavedec(V_seg(idx, :), levMax, wavType); % wavelet transform

figure('Position', [200, 400, 400, 150]);
plot(wav, 'k-', 'LineWidth', 1); xlim([-inf, inf]); ylim([-0.02, 0.02]);
yline(0, 'Color', [0.5, 0.5, 0.5]);

thresh      = prctile(abs(wav), 95); wav(abs(wav) < thresh) = 0; % Keep only top 5% of coefficients
figure('Position', [200, 200, 400, 150]);
plot(wav, 'k-', 'LineWidth', 1); xlim([-inf, inf]); ylim([-0.02, 0.02]);
yline(0, 'Color', [0.5, 0.5, 0.5]);

%% PCA on all CCEPs (no subject withheld) - determine which PCs to use, plot all CCEPs on PCA space
% SVD on centered data (centered across time points)
[U, S, V]       = svd((V_wav - mean(V_wav))', 'econ');
score           = V*S; % weigh right singular vectors by S to get the PCA scores. Rows are observations, columns correspond to PCs. Note that the data can be reconstructed as score*U'
varExp          = diag(S^2); varExp = 100*varExp/sum(varExp);

% variance explained
figure('Position', [200, 200, 200, 600]); plot(varExp(1:20), 'k-o', 'MarkerFaceColor', 'k');
set(gca, 'ytick', 0:10:50);
title('Variance explained by PCs');
xlim([0.5, 5.5]); ylim([0, 50]);
% % Save figure
% saveas(gcf, fullfile('output', 'varExp'), 'png');
% saveas(gcf, fullfile('output', 'varExp'), 'svg');

% kernel densities of distributions of each PC
xbin            = linspace(-2, 2, 30);
figure('Position', [400, 200, 200, 600]);
for ii = 1:5
    subplot(5, 1, ii);
    pdf1 = hist(score(pathIdxes.idxes{1}, ii), xbin); % AmygACC
    pdf2 = hist(score(pathIdxes.idxes{5}, ii), xbin); % hippPCC
    plot(xbin, pdf1, 'Color', pathIdxes.colors(1, :), 'LineWidth', 1.5); hold on
    plot(xbin, pdf2, 'Color', pathIdxes.colors(5, :), 'LineWidth', 1.5); hold off
    xlim([-2, 2]); ylim([0, 30])
end
% saveas(gcf, fullfile('output', 'PCdistributions'), 'png');
% saveas(gcf, fullfile('output', 'PCdistributions'), 'svg');

% find the 2 PCs between 1-5 that differentiate most between AmygACC and hippPCC conditions
tstats      = zeros(5, 1);
for ii = 1:5
    dist1               = score(pathIdxes.idxes{1}, ii);
    dist2               = score(pathIdxes.idxes{5}, ii);
    [~, ~, ~, stats]    = ttest2(dist1, dist2);
    tstats(ii)          = abs(stats.tstat);
end
[~, order]      = sort(tstats, 'descend');
PCs2Use         = order(1:2); % PCs to plot and use for LDA
PCs2Use         = sort(PCs2Use);
fprintf('Using PC %d and PC %d\n', PCs2Use(:));

figure; hold on
for ii = 1:5
    plot(score(pathIdxes.idxes{ii}, PCs2Use(1)), score(pathIdxes.idxes{ii}, PCs2Use(2)), '.', 'Color', pathIdxes.colors(ii, :), 'MarkerSize', 16, 'DisplayName', pathIdxes.name{ii});
end
xlabel(sprintf('PC %d', PCs2Use(1))); ylabel(sprintf('PC %d', PCs2Use(2)));
axis equal;
legend;

exc = setdiff(1:length(labels_allsubs), cat(1, pathIdxes.idxes{:})); % all non-labelled indices
plot(score(exc, 1), score(exc, 2), '.', 'Color', [0.5, 0.5, 0.5]);

%% Find decision boundary on all CCEP model

training    = [score(pathIdxes.idxes{1}, PCs2Use); score(pathIdxes.idxes{5}, PCs2Use)]; % amygACC first, hippPCC second
group       = [ones(length(pathIdxes.idxes{1}), 1); 5*ones(length(pathIdxes.idxes{5}), 1)]; % CCEP type to use for classification

[class, ~, ~, ~, coeff] = classify(training, training, group);      % find decision boundary coefficients for model fit on all data (without leaving any out)
accTrain                = 100*sum(class == group) / length(class);  % training accuracy
fprintf('Training accuracy: %0.2f%%\n', accTrain);

f           = @(x, y) coeff(1, 2).const + coeff(1, 2).linear(1)*x + coeff(1, 2).linear(2)*y; % 0 = K + x*L(1) + y*L(2);

figure('Position', [200, 200, 600, 600]); hold on

for ii = 1:5
    plot(score(pathIdxes.idxes{ii}, PCs2Use(1)), score(pathIdxes.idxes{ii}, PCs2Use(2)), '.', 'Color', pathIdxes.colors(ii, :), 'MarkerSize', 16, 'DisplayName', pathIdxes.name{ii});
end
plot(score(exc, 1), score(exc, 2), '.', 'Color', [0.5, 0.5, 0.5], 'DisplayName', 'Other CCEPs');
    
h           = fimplicit(f, [-3, 3], '--', 'LineWidth', 1.5, 'Color', [0.3, 0.3, 0.3]); % decision boundary
h.DisplayName= 'Decision boundary between amygACC and hippPCC';
axis equal;
xlim([-2, 2]); ylim([-2, 2]);
set(gca, 'xtick', -3:3); set(gca, 'ytick', -2:2);
xlabel(sprintf('PC %d', PCs2Use(1))); ylabel(sprintf('PC %d', PCs2Use(2)));
hold off

%% Leave-one-subject out CV analysis to calculate accuracy

accSubs     = nan(nsubs, 1); % CV test accuracy for each subject withheld

% store the PCA scores at each CV fold
scoreTrainFolds = cell(nsubs, 1);
scoreTestFolds  = cell(nsubs, 1);

for sub = 1:nsubs

    % divide into test and train sets for DWT CCEPs and labels
    V_wavTest           = V_wav(subNum == sub, :);
    V_wavTrain          = V_wav(subNum ~= sub, :);
    labels_test         = labels_allsubs(subNum == sub, :);
    labels_train        = labels_allsubs(subNum ~= sub, :);

    % PCA on the training set
    [U, S, V]           = svd((V_wavTrain - mean(V_wavTrain))', 'econ');
    scoreTrain          = V*S;

    scoreTest           = (V_wavTest - mean(V_wavTest))*U; % scores of testing set (mean centered projection onto loadings)
    
    % save the scores per fold so they don't need to be recalculated for each perm. Just store the PCs2Use
    scoreTrainFolds{sub} = scoreTrain(:, PCs2Use);
    scoreTestFolds{sub} = scoreTest(:, PCs2Use);

    % scores and groupings for training data with amygACC (first) and hippPCC (second)
    idxesTrain          = getIdxes(labels_train);
    training            = [scoreTrain(idxesTrain{1}, PCs2Use); scoreTrain(idxesTrain{5}, PCs2Use)];
    groupTrain          = [ones(length(idxesTrain{1}), 1); 5*ones(length(idxesTrain{5}), 1)];

    % scores and groupings for test data
    idxesTest           = getIdxes(labels_test);
    testing             = [scoreTest(idxesTest{1}, PCs2Use); scoreTest(idxesTest{5}, PCs2Use)];
    groupTest           = [ones(length(idxesTest{1}), 1); 5*ones(length(idxesTest{5}), 1)];

    %class = classify(testing, training, groupTrain); % find decision boundary coefficients for model fit on all data (without leaving any out)
    [class, ~, ~, ~, coeff] = classify(testing, training, groupTrain); % find decision boundary coefficients for model fit on all data (without leaving any out)
    accSubs(sub)            = 100 * sum(class == groupTest) / length(class); % training accuracy

    fprintf('Accuracy for test sub %d: %0.2f%%\n', sub, accSubs(sub));

end

fprintf('Mean leave-one-subject-out CV accuracy: %0.2f%%\n', mean(accSubs));

%% Calculate null accuracy distribution for each CV fold by permuting CCEP labels within each subject
% Permutation should be done within each subject:
% https://stats.stackexchange.com/questions/536308/permutation-testing-for-machine-learning-permute-entire-set-or-only-training-se
%   More conservative; controls for subject-specific differences in waveform. For example, if one subject only has accAmyg CCEPs, the null accuracy would equal
%   the real accuracy, because there is no way to know if the waveform effects in that subject were due to subject or connection type

rng('default');

nperm       = 100000;
accSubsPerm = nan(nperm, nsubs); % stores CV accuracy for each subject at each permutation
for ii = 1:nperm
    if ~mod(ii, nperm/100), fprintf('.'); end % print dot every 1%
    
    % CCEP labels are randomly permuted within each subject. More fair if each subject has different proportions of each CCEP
    labels_allsubsPerm = cell(size(labels_allsubs));
    for sub = 1:nsubs
        labels_thissub = labels_allsubs(subNum == sub, :);
        labels_allsubsPerm(subNum == sub, :) = labels_thissub(randperm(labels_length(sub)), :);
    end

    for sub = 1:nsubs

        % get training and test subject labels from permuted set
        labels_train            = labels_allsubsPerm(subNum ~= sub, :);
        labels_test             = labels_allsubsPerm(subNum == sub, :);

        % scores and groupings for training data using permuted labels
        idxesTrain              = getIdxes(labels_train);
        training                = [scoreTrainFolds{sub}(idxesTrain{1}, :); scoreTrainFolds{sub}(idxesTrain{5}, :)]; % using scoreTrainFolds from previous section
        groupTrain              = [ones(length(idxesTrain{1}), 1); 5*ones(length(idxesTrain{5}), 1)];

        % scores and groupings for test data using permuted labels
        idxesTest               = getIdxes(labels_test);
        testing                 = [scoreTestFolds{sub}(idxesTest{1}, :); scoreTestFolds{sub}(idxesTest{5}, :)];
        groupTest               = [ones(length(idxesTest{1}), 1); 5*ones(length(idxesTest{5}), 1)];

        class                   = classify(testing, training, groupTrain); % get LDA classification
        accSubsPerm(ii, sub)    = 100*sum(class == groupTest) / length(class); % training accuracy
    end
    
end
fprintf('\n');

%% Plot CV accuracy compared to null accuracies
% p value of each CV fold
p = nan(nsubs, 1);
for sub = 1:nsubs
    p(sub) = sum(accSubsPerm(:, sub) >= accSubs(sub)) / nperm; % right-tailed p-value (H_a: CV accuracy is GREATER than expected by chance)
end
% h = p < 0.05/nsubs; % Bonferroni-corrected
h = fdr_bh(p, 0.05, 'dep');

% Real CV accuracy overlaid with null accuracy +/- SD
figure('Position', [200, 200, 600, 400]); hold on
bar(1:nsubs, accSubs, 'FaceColor', [0.8, 0.8, 0.8]);
errorbar(1:nsubs, mean(accSubsPerm), std(accSubsPerm), 'k.', 'MarkerSize', 18); % mean +/i SD of null accuracy
yline(mean(accSubs), 'Color', 'r', 'LineWidth', 1.5);
text(find(h), 95*ones(sum(h), 1), '*', 'Color', 'k', 'FontSize', 24, 'HorizontalAlignment', 'center');
xlim([0, nsubs+1]); ylim([0, 100]);
set(gca, 'ytick', 0:20:100);
xlabel('Cross-validation trial'); ylabel('Accuracy (%)');

figure('Position', [400, 200, 400, 100]); hold on
histogram(mean(accSubsPerm, 2), 30);
xline(mean(accSubs), 'Color', 'r', 'LineWidth', 1.5);
xlim([0, 100]);

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
