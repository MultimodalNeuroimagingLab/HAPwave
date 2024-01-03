clearvars, clc, close all;

%% Load subjects variables
% Script to plot single subject pannels with cross correlation
% cd to HAPwave repository
addpath(genpath(pwd))

% set local path to your BIDS directory:
myPath          = setLocalDataPath(1);
localDataPath   = myPath.input;

% load the meta data
all_subjects    = {'01','02','03','04','05','06','07','08'}; % 
all_hemi        = {'r','r','r','l','r','l','l','r'};
all_runs        = {'01','01','01','01','01','01','01','01'};

% loop through subjects
for ss              = 1:length(all_subjects)
    bids_sub        = all_subjects{ss};
    bids_ses        = 'ieeg01';
    bids_task       = 'ccep';
    bids_run        = all_runs{ss};

    [events_table,channels_table,electrodes_table,sub_out] = pcc_loadAveragesStats(localDataPath,bids_sub,bids_ses,bids_task,bids_run);
    all_out(ss)     = sub_out;
end

%% Correct P values for number of comparisons in each subject
area_codes      = {[12123 53 54 12108 12109 12110 12106 12107 11123 49 17 18 11108 11109 11110 11106 11107 10]}; % all areas
nr_subs         = length(all_subjects);

for ss = 1:nr_subs % loop over subjects
   
    % Get sites that belong to the measurement ROI (rec_area)
    these_measured_sites    = find(ismember(all_out(ss).channel_areas,area_codes{1}));
    
    % Get sites that belong to the stimulated ROI (stim_area)
    these_stim_sites        = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{1}) | ...
                                   ismember(all_out(ss).average_ccep_areas(:,2),area_codes{1}));
    
    % p-values for correction of multiple comparisons
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
            if ~isempty(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).data) % ~ same stim/recording site
                all_out(ss).hasdata(these_measured_sites(kk), these_stim_sites(ll)) = 1;
                all_out(ss).crp_p(these_measured_sites(kk), these_stim_sites(ll))   = all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_projs.p_value_tR;
                all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)) = mean(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.al_p); % mean alpha prime across trials;
                all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll))     = median(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.cod); % median across trials
                sig_timepoints = find(all_out(ss).tt>0.015 & all_out(ss).tt<=all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.tR);
                all_out(ss).avg_trace_tR(these_measured_sites(kk), these_stim_sites(ll),sig_timepoints) = all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.C;

            else
                all_out(ss).hasdata(these_measured_sites(kk), these_stim_sites(ll)) = 0;
            end
        end
    end
    pvals                                           = all_out(ss).crp_p(all_out(ss).hasdata==1);
    qq                                              = 0.05;
    [h, crit_p, adj_ci_cvrg, adj_p]                 = fdr_bh(pvals,qq,'dep','no');
    all_out(ss).crp_p_adj(all_out(ss).hasdata==1)   = adj_p;
    all_out(ss).h(all_out(ss).hasdata==1)           = h;
end

%% Plot mean of Cross-correlations between ANT and HC stim
% subject 2 and 7 have ANT coverage
% set codes acording to hemis
area_codes_r    = {[12123 53],[54],[12108 12109 12110],[12106 12107],[49]}; % right
area_codes_l    = {[11123 17],[18],[11108 11109 11110],[11106 11107],[10]}; % left
input_color     = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.6350 0.0780 0.1840],[0.8500 0.3250 0.0980]}; % Blue for hippocampus, Green for ANT, red for confidence interval in time
sub_hemi        = {'r','r','r','l','r','l','l','r'};

% ______________________________________________________________________ %
% %% Set subject to test %%
ss              = 7; % Set subject 2 or 7, with ANT electrodes

if isequal(sub_hemi{ss},'l')
    area_codes = area_codes_l;
elseif isequal(sub_hemi{ss},'r')
    area_codes = area_codes_r;
end

% start plotting
figure('Position',[0 0 500 600]),hold on
subplot(2,1,1),hold on

% %% Set connection %%
% Keep in mind area codes:
% Hipp = 1; Amyg = 2; PCC = 3; ACC = 4; ANT = 5   
% Let's start with HC -> PCC connection
measure_ind     = 3;    % measuring in PCC
stim_ind        = 1;    % HC stim

% Get sites in measured & stim ROI
these_measured_sites    = find(ismember(all_out(ss).channel_areas,area_codes{measure_ind}));
these_stim_sites        = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{stim_ind}) | ...
                               ismember(all_out(ss).average_ccep_areas(:,2),area_codes{stim_ind}));

sign_resp       = all_out(ss).crp_p_adj < 0.05;   % plot p<0.05 FDR corrected
sign_resp(all_out(ss).bad_channels,:) = 0;      % bad channels cannot be significant
t_win_cod       = [0.015 0.500];                % window for vector length normalization and plotting across subjects

% get HC->PCC significant waveforms
resp_counter    = 0;    % start counter
resp1_hip2pc    = [];   % prepare matrix with HC->PCC waves

% walk through pairs, test significance and plot if so
 for kk = 1:length(these_measured_sites)
     for ll = 1:length(these_stim_sites)
         if sign_resp(these_measured_sites(kk),these_stim_sites(ll))==1 % significant wave found
             resp_counter   = resp_counter+1;                           % counting
             plot_responses = squeeze(all_out(ss).average_ccep(these_measured_sites(kk),these_stim_sites(ll),:));
             plot_responses(all_out(ss).tt>-0.010 & all_out(ss).tt<0.010) = NaN; % remove stim artifact

             % scaling to unit length(https://en.wikipedia.org/wiki/Feature_scaling)
             % unti length taken in same window as stats
             response_vector_length = sum(plot_responses(all_out(ss).tt>t_win_cod(1) &  all_out(ss).tt<t_win_cod(2)).^2).^.5;
             plot_responses_norm    = plot_responses./(response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial
             
             plot(all_out(ss).tt, ss*.2 + zeros(size(all_out(ss).tt)),'k')
             plot(all_out(ss).tt, ss *.2 + plot_responses_norm,'Color',input_color{4},'LineWidth',1)
             resp1_hip2pc(resp_counter,:) = plot_responses;
         end
     end
 end
xlim([-0.2 .6])%,ylim([-2000 2000])
xlabel('Time (s)')
ylabel('Normalized Voltage')


% %% Set connection %%
% get ANT -> PCC connection
stim_ind        = 5;    % changing to ANT stim
resp_counter    = 0;    % restart counter
resp1_ant2pc    = [];   % prepare matrix with ANT->PCC waves



% Get sites that belong to the measured & stim ROI
these_measured_sites    = find(ismember(all_out(ss).channel_areas,area_codes{measure_ind}));
these_stim_sites        = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{stim_ind}) | ...
                               ismember(all_out(ss).average_ccep_areas(:,2),area_codes{stim_ind}));

sign_resp       = all_out(ss).crp_p_adj<0.05; % plot p<0.05 FDR corrected
sign_resp(all_out(ss).bad_channels,:) = 0; % bad channels can not be significant


% walk through pairs, test significance and plot if so
for kk = 1:length(these_measured_sites)
    for ll = 1:length(these_stim_sites)
        if sign_resp(these_measured_sites(kk),these_stim_sites(ll))==1  % significant wave found
            resp_counter    = resp_counter+1;                           % counting
            plot_responses  = squeeze(all_out(ss).average_ccep(these_measured_sites(kk),these_stim_sites(ll),:));
            plot_responses(all_out(ss).tt>-0.010 & all_out(ss).tt<0.010) = NaN; % remove stim artifact from plot
            
            % Normalize. Scaling to unit length (Euclidean lenght): https://en.wikipedia.org/wiki/Feature_scaling
            % unit length taken in same window as stats
            response_vector_length  = sum(plot_responses(all_out(ss).tt > t_win_cod(1) &  all_out(ss).tt < t_win_cod(2)) .^ 2) .^ .5;
            plot_responses_norm     = plot_responses ./ (response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial
            
            plot(all_out(ss).tt,.15 + ss*.2 + zeros(size(all_out(ss).tt)),'k')
            plot(all_out(ss).tt,.15 + ss*.2 + plot_responses_norm,'Color',input_color{2},'LineWidth',1)
            resp1_ant2pc(resp_counter,:) = plot_responses;
        end
    end
end
title(['Sub-0' num2str(ss)]);

subplot(2,1,2), hold on
c1_all  = [];
nn      = 1;
for kk = 1:size(resp1_ant2pc,1)
    for ll              = 1:size(resp1_hip2pc,1)
        [c1,lags]       = xcorr(resp1_ant2pc((kk), all_out(ss).tt >0.015)',resp1_hip2pc((ll), all_out(ss).tt> 0.015)','coeff');
        c1_all(nn,:)    = c1;                       % getting the cross-correlations
        tt_lags         = lags/all_out(ss).srate;   % getting lags in time
        tt_set          = all_out(ss).tt(all_out(ss).tt>0.015);
        plot(tt_lags,c1,'Color',[.5 .5 .5 .2],'LineWidth',2)
        nn              = nn + 1;
    end
end
plotCurvConf(tt_lags, c1_all); % plots 95% Confidende interval in gray with 50% transparency

mean_c1 = mean(c1_all);
plot(tt_lags,mean_c1,'color','k','LineWidth',1)

% Bootstrapping and resampling time shift.
nr_bootstraps   = 10000;
t_shifts        = NaN(nr_bootstraps,1);
for kk = 1:nr_bootstraps
    this_subset     = randi(size(c1_all,1),size(c1_all,1),1);
    [~,ind_c1]      = max(mean(c1_all(this_subset,:),1));
    t_shifts(kk)    = tt_lags(ind_c1);
end

xline([quantile(t_shifts,.075) quantile(t_shifts,.925)],'-',{median(t_shifts), ' '},'LabelVerticalAlignment','bottom',...
    'LabelOrientation','horizontal', 'LabelHorizontalAlignment','left',...
    'Color',input_color{3},'LineWidth',1) % to rotate label use: 'LabelOrientation','horizontal' | 'vertical'(default),

[max_c1] = max(mean(c1_all,1));
title(['Cross correlation' string(max_c1)])
xlim([-1 1]), ylim([-0.6 1.2])
ylabel('Normalized Cross Correlation')
xlabel('Time shift (s)')

