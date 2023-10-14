clearvars, close all, clc
% startup

%% Load stats across subs
% Dependencies: MNL ieeg basics repository
% cd to HAPwave repository
addpath(genpath(pwd))

% set local path to your BIDS directory:
myPath          = setLocalDataPath(1);
localDataPath   = myPath.input;

% load the meta data
all_subjects    = {'01','02','03','04','05','06','07','08'}; % 
all_hemi        = {'r','r','r','l','r','l','l','r'};
all_runs        = {'01','01','01','01','01','01','01','01'};

% load the meta data
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

for ss = 1:nr_subs % loop over subject
   
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

%% Load stim and measure areas across subjects
% Get area codes of ROI by hemisphere
area_names      = {'Hipp','Amyg','PCC','ACC','ANT'};   
area_codes_r    = {[12123 53],[54],[12108 12109 12110],[12106 12107],[49]}; % right
area_codes_l    = {[11123 17],[18],[11108 11109 11110],[11106 11107],[10]}; % left

nr_subs         = length(all_subjects);
out             = [];                       % prepare areaXarea matrix with all subjects concatinated

subj_resp_total = zeros(nr_subs,1);         % set counter for significant responses at 0

t_win_norm      = [0.015 0.500];            % window for vector length normalization and plotting across subjects

for measure_ind = 1:length(area_codes_r) % loop through the inds 
    for stim_ind = 1:length(area_codes_r)

        resp_counter = 0; % counting all responses across subjects for this connection

        for ss = 1:nr_subs % subject loop
            if isequal(all_hemi{ss},'l')
                area_codes = area_codes_l;
            elseif isequal(all_hemi{ss},'r')
                area_codes = area_codes_r;
            end

            % Get sites that belong to the measurement ROI (rec_area)
            these_measured_sites    = find(ismember(all_out(ss).channel_areas,area_codes{measure_ind}));
            
            % Get sites that belong to the stimulated ROI (stim_area)
            these_stim_sites        = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{stim_ind}) | ...
                                           ismember(all_out(ss).average_ccep_areas(:,2),area_codes{stim_ind}));
            
            % loop over measured sites
            for kk = 1:length(these_measured_sites)
                % loop over the stimulated pairs
                for ll = 1:length(these_stim_sites)
                    if ~isempty(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).data) % ~ same stim/recording site
                        
                        % this is a pair, counting for multiple comparison
                        % correction per subject
                        subj_resp_total(ss) = subj_resp_total(ss) + 1; % set counter
                        
                        % first raw responses
                        plot_responses      = squeeze(all_out(ss).average_ccep(these_measured_sites(kk), these_stim_sites(ll), :));
                        tt                  = all_out(ss).tt;

                        % save outputs
                        resp_counter        = resp_counter + 1;

                        % get CCEP responses for plotting
                        % Scaling to unit length (Euclidean lenght): https://en.wikipedia.org/wiki/Feature_scaling
                        % unit length taken in same window as stats
                        response_vector_length  = sum(plot_responses(all_out(ss).tt > t_win_norm(1) &  all_out(ss).tt < t_win_norm(2)) .^ 2) .^ .5;
                        plot_responses_norm     = plot_responses ./ (response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial
                        out(measure_ind,stim_ind).plot_responses_norm(resp_counter, :) = plot_responses_norm;
                        out(measure_ind,stim_ind).plot_responses(resp_counter, :) = plot_responses;

                        % store subject index
                        out(measure_ind,stim_ind).subj_ind(resp_counter, :) = ss;
                        
                        % save CRP stuff for 
                        out(measure_ind,stim_ind).p(resp_counter, :)        = all_out(ss).crp_p_adj(these_measured_sites(kk), these_stim_sites(ll));
                        out(measure_ind,stim_ind).cod(resp_counter, :)      = all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll)); 
                        out(measure_ind,stim_ind).a_prime(resp_counter, :)  = all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)); 
                    end
                end 
            end
        end 
    end 
end


%% plot average waveforms at per subject (normalized)
sub_color = {[0 0.4470 0.7410],...           % blue
             [0.8500 0.3250 0.0980],...      % orange
             [0.4660 0.6740 0.1880],...      % green
             [0.4940 0.1840 0.5560],...      % purple
             [0.9290 0.6940 0.1250],...      % mustard
             [0.3010 0.7450 0.9330],...      % celeste
             [0.6350 0.0780 0.1840],...      % wine
             [1 0 1]};                       % magenta
    
% Set connection to plot with stim and meas ROI:
% HC= 1; Amg= 2; PCC= 3; ACC= 4
stim_ind        = 1;    % Stimulated ROI  
measure_ind     = 3;    % Measured ROI
% tt              = all_out(1).tt;

figure('Position',[0 0 400 800]), hold on

for ss = 1:8 
    % only this subject
    ss_resps    = out(measure_ind,stim_ind).subj_ind==ss;

    % only significant responses (FDR corrected)
    sign_resp   = out(measure_ind,stim_ind).p<0.05; % adjusted for multiple comparisons

    this_set    = out(measure_ind,stim_ind).plot_responses_norm(sign_resp & ss_resps,:)';        

    % Look for coverage in the stimulated and/or measurement ROI
    if isempty(this_set)
        % no coverage sites of interest 
        % draw dotted line             / Meaning: No electrodes implanted in ROI for this sub
        yline(ss * .2,':', 'Color', [.9 .3 .3],'LineWidth', .5);
        
    else
        % Coverage in both measurement and stimulated ROI
        % black line if no sign responses     / Meaning: CCEPs were not significant under SEPS of interest ):
        plot([-.2 .8], [ss * .2, ss * .2], 'Color', 'k', 'LineWidth', .7);
        
        this_set(tt > -0.010 & tt < 0.010,:) = NaN;
        plot(tt, ss * .2 + this_set, 'Color',sub_color{ss}, 'LineWidth', .5);
        xlim([-0.2 .6]), ylim([0 1.7]);
        ylabel('Voltage (L2-normalized)')
        xlabel('Time (s)')
       
        ss_sign     = width(this_set);
        txt         = 100 * ss_sign/sum(ss_resps)    
        text(-.1, ss * .209, num2str(txt))
        title([area_names{stim_ind} ' -> ' area_names{measure_ind}])
    end
end


