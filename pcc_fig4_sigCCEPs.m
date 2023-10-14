clearvars, close all, clc

%% Plot waveforms for each subject in uV (figure 4)
% Dependencies: mnl_ieeg_basics, matmef, and vistasoft github repositories. 
% cd to HAPwave repo
addpath(genpath(pwd)) 

% set local path to your BIDS directory:
myPath          = setLocalDataPath(1);
localDataPath   = myPath.input;

% load the metadata
all_subjects    = {'01','02','03','04','05','06','07','08'};        % List of subjects
all_hemi        = {'r','r','r','l','r','l','l','r'};                % List of hemispheres
all_runs        = {'01','01','01','01','01','01','01','01'};        % List of runs

for ss              = 1:length(all_subjects)
    bids_sub        = all_subjects{ss};
    bids_ses        = 'ieeg01';
    bids_task       = 'ccep';
    bids_run        = all_runs{ss};

    % Load metadata and stats
    [events_table,channels_table,electrodes_table,sub_out] = ...
        pcc_loadAveragesStats(localDataPath,bids_sub,bids_ses,bids_task,bids_run);
    all_out(ss) = sub_out;
end


%% Correct P values for number of comparisons in each subject
% Fisrt load the limbic codes in both right and left hemis
area_codes  = {[12123 53 54 12108 12109 12110 12106 12107 59 ...
                11123 17 18 11108 11109 11110 11106 11107 10]};
nr_subs     = length(all_subjects);

for ss = 1:nr_subs  % loop over subjects
   
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
            if ~isempty(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).data)   % is a limbic connection
                all_out(ss).hasdata(these_measured_sites(kk), these_stim_sites(ll)) = 1;            % is significant CRP-based
                all_out(ss).crp_p(these_measured_sites(kk), these_stim_sites(ll))   = all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_projs.p_value_tR;  % p-value
                all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)) = mean(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.al_p);  % mean 
                all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll))     = median(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.cod); % CRP coefficent of determination
                sig_timepoints = find(all_out(ss).tt>0.015 & all_out(ss).tt<=all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.tR);                      % response duration
                all_out(ss).avg_trace_tR(these_measured_sites(kk), these_stim_sites(ll),sig_timepoints) = all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.C; % al_p for tR
            else
                all_out(ss).hasdata(these_measured_sites(kk), these_stim_sites(ll))     = 0;         % is not significant
            end
        end
    end
    pvals                                           = all_out(ss).crp_p(all_out(ss).hasdata==1);  % get pVals
    qq                                              = 0.05;                         % false discovery rate
    [h, crit_p, adj_ci_cvrg, adj_p]                 = fdr_bh(pvals,qq,'dep','no');  % Benjamini & Yekutieli FDR correction 
    all_out(ss).crp_p_adj(all_out(ss).hasdata==1)   = adj_p;                        % Adjusted pVals
    all_out(ss).h(all_out(ss).hasdata==1)           = h;                            % adjusted pVal is significant
end


%% Load stim and measure areas across subjects
% Sort limbic codes by hemisphere
area_names      = {'Hipp','Amyg','PCC','ACC'};   
area_codes_r    = {[12123 53],[54],[12108 12109 12110],[12106 12107]};  % right
area_codes_l    = {[11123 17],[18],[11108 11109 11110],[11106 11107]};  % left
out             = [];                   % prepare areaXarea structure with all subjects concatinated
subj_resp_total = zeros(nr_subs,1);     % set counter for significant responses at 0
t_win_norm      = [0.015 0.500];        % time window for pre-processing

for measure_ind = 1:length(area_names)  % loop through measured sites
    for stim_ind = 1:length(area_names) % now go through stimulated sites
        resp_counter = 0;               % counting all responses across subjects for this connection

        for ss = 1:nr_subs % loop over subjects
            tt = all_out(ss).tt;
            % which hemisphere has coverage
            if isequal(all_hemi{ss},'l')
                area_codes  = area_codes_l;
            elseif isequal(all_hemi{ss},'r')
                area_codes  = area_codes_r;
            end

            % Get recording ROI (measured_area)
            these_measured_sites    = find(ismember(all_out(ss).channel_areas,area_codes{measure_ind}));
           
            % Get stimulated ROI (stim_area)
            these_stim_sites        = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{stim_ind})...
                                         | ismember(all_out(ss).average_ccep_areas(:,2),area_codes{stim_ind}));
            
            % loop over measured sites
            for kk = 1:length(these_measured_sites)
                % loop over the stimulated pairs
                for ll = 1:length(these_stim_sites)
                    if ~isempty(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).data) % ~ same stim/recording site
                        
                        % this is a pair, counting for multiple comparison correction per subject
                        subj_resp_total(ss)     = subj_resp_total(ss) + 1;          % set counter
                        
                        % first raw responses
                        plot_responses          = squeeze(all_out(ss).average_ccep(these_measured_sites(kk), these_stim_sites(ll), :));

                        % save outputs
                        resp_counter            = resp_counter + 1;                 % count significant wave
                        out(measure_ind,stim_ind).plot_responses(resp_counter, :)   = plot_responses;   % save significant  waveform
                        out(measure_ind,stim_ind).subj_ind(resp_counter, :)         = ss;               % get subject's index
                        
                        % save CRP stuff for 
                        out(measure_ind,stim_ind).p(resp_counter, :)        = all_out(ss).crp_p_adj(these_measured_sites(kk), these_stim_sites(ll));    % save adjusted p
                        out(measure_ind,stim_ind).cod(resp_counter, :)      = all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll));          % save COD
                        out(measure_ind,stim_ind).a_prime(resp_counter, :)  = all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll));      % save alpha prime
                    end
                end 
            end
        end
    end 
end


%% plot average waveforms at per subject (voltage in uV)
sub_axis = {[-250 700],[-110 300],[-350 300],[-750 2100],[-200 350],[-500 1000],[-300 350],[-250 750]}; % each subject has different y-axis
sub_color = {[0 0.4470 0.7410],...  % blue
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

figure('Position',[0 0 400 800]), hold on

for ss = 1:8 
    subplot(8,1,(ss)), hold on

    % only for this subject
    ss_resps    = out(measure_ind,stim_ind).subj_ind==ss;
    % only significant responses (CRPs FDR adjusted)
    sign_resp   = out(measure_ind,stim_ind).p < 0.05; % adjusted for multiple comparisons
    % find waveforms in this subject 
    this_set    = out(measure_ind,stim_ind).plot_responses(sign_resp & ss_resps,:)';

    % Look for coverage in the stimulated & measured ROI
    if isempty(this_set)
            % no coverage for one or both ROIs: redish dotted line
        yline(ss,':', 'Color', [.9 .3 .3],'LineWidth', .5);    
    else    % Coverage in stim+meas ROI: black line (if black line only, no sign responses for this connection
        plot([-.2 .8], [ss, ss], 'Color', 'k', 'LineWidth', .7);
        this_set(tt > -0.010 & tt < 0.010,:) = NaN;         % remove stim artifact

        % plot waveforms in uV
        plot(tt, ss + this_set, 'Color',sub_color{ss}, 'LineWidth', .5);
        xlim([-0.2 .6]), ylim([sub_axis{ss}]);
        ylabel('Voltage (uV)')

        % get proportion of sign responses
        ss_sign     = width(this_set);
        txt         = 100 * ss_sign /sum(ss_resps)          % percentage of significant responses
        text(-.1, 200, num2str(txt)),
        title([area_names{stim_ind} ' -> ' area_names{measure_ind}])
    end
end
 xlabel('Time (s)')
