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
    bids_sub    = all_subjects{ss};
    bids_ses    = 'ieeg01';
    bids_task   = 'ccep';
    bids_run    = all_runs{ss};

    % Load metadata and stats
    [events_table,channels_table,electrodes_table,sub_out] = pcc_loadAveragesStats(localDataPath,bids_sub,bids_ses,bids_task,bids_run);
    all_out(ss) = sub_out;
end


%% Correct P values for number of comparisons in each subject
% Fisrt load the limbic codes in both right and left hemis
area_codes  = {[12123 53 54 12108 12109 12110 12106 12107 59 ...
                11123 17 18 11108 11109 11110 11106 11107 10]}; 
nr_subs     = length(all_subjects);

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
    pvals                                           = all_out(ss).crp_p(all_out(ss).hasdata==1);  % get pVals
    qq                                              = 0.05; % false discovery rate
    [h, crit_p, adj_ci_cvrg, adj_p]                 = fdr_bh(pvals,qq,'dep','no'); % Benjamini & Yekutieli FDR correction 
    all_out(ss).crp_p_adj(all_out(ss).hasdata==1)   = adj_p; % Adjusted pVals
    all_out(ss).h(all_out(ss).hasdata==1)           = h; % adjusted pVal is significant
end

%% Load stim and measure sites and find N1 responses
% Sort limbic codes by hemisphere
area_names      = {'Hipp','Amyg','PCC','ACC'};   
area_codes_r    = {[12123 53],[54],[12108 12109 12110],[12106 12107]}; % right
area_codes_l    = {[11123 17],[18],[11108 11109 11110],[11106 11107]}; % left

out             = [];                   % prepare an AreaByArea structure, with all subjects concatinated for each area
subj_resp_total = zeros(nr_subs,1);     % prepare for total of adjusted FDR per connection (stim->measured)
t_win_norm      = [0.015 0.500];        % window for vector length normalization and plotting across subjects
tt              = all_out(ss).tt;

for measure_ind = 1:length(area_names)  % loop through measured sites
    for stim_ind = 1:length(area_names) % now go through stimulated sites

        resp_counter = 0; % count all responses across subjects for this connection (stim->measured pair)

        for ss = 1:nr_subs              % loop over subjects
            % which hemisphere has coverage
            if isequal(all_hemi{ss},'l')
                area_codes = area_codes_l;
            elseif isequal(all_hemi{ss},'r')
                area_codes = area_codes_r;
            end

            % Get sites in recording ROI (measured_area)
            these_measured_sites    = find(ismember(all_out(ss).channel_areas,area_codes{measure_ind}));
           
            % Get sites in stimulated ROI (stim_area)
            these_stim_sites        = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{stim_ind}) | ...
                                           ismember(all_out(ss).average_ccep_areas(:,2),area_codes{stim_ind}));

            % loop over measured sites
            for kk = 1:length(these_measured_sites)
                % loop over the stimulated pairs
                for ll = 1:length(these_stim_sites)
                    if ~isempty(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).data) 

                        subj_resp_total(ss)     = subj_resp_total(ss) + 1; % counting pair for multiple comparison correction per subject
                        
                        % first get raw responses
                        plot_responses          = squeeze(all_out(ss).average_ccep(these_measured_sites(kk), these_stim_sites(ll), :));

                        % Is there an N1?
                        params.amplitude_thresh = 3.4;              % standard deviations 3.4;
                        params.n1_peak_range    = 0.050;            % in s
                        params.baseline_tt      = tt>-.5 & tt<-.020;
                        params.peakSign         = 1;
                        % amplitudeThresh is in microvolts, used if standard devisation from baseline is smaller than this value 
                        % threshold for responses is (amplitude_thresh * amplitudeThresh) is baseline STD < amplitudeThresh
                        % threshold for responses is (amplitude_thresh * baseline STD) is baseline STD > amplitudeThresh
                        params.amplitudeThresh  = 10;

                        [posn1_peak_sample,posn1_peak_amplitude,posn1_peak_time] = ccep_detect_n1peak_sEEG(plot_responses,tt,params);
                        params.peakSign         = -1;
                        [negn1_peak_sample,negn1_peak_amplitude,negn1_peak_time] = ccep_detect_n1peak_sEEG(plot_responses,tt,params);
                        
                        % Keep earliest peak
                        if ~isnan(posn1_peak_time) && ~isnan(negn1_peak_time) % if both positive and negative
                            if posn1_peak_time < negn1_peak_time        % positive is earlier
                                n1_peak_amplitude   = posn1_peak_amplitude;
                                n1_peak_time        = posn1_peak_time;
                                n1_peak_sample      = posn1_peak_sample;
                            elseif negn1_peak_time < posn1_peak_time    % negative is earlier
                                n1_peak_amplitude   = negn1_peak_amplitude;
                                n1_peak_time        = negn1_peak_time;
                                n1_peak_sample      = negn1_peak_sample;
                            end
                        elseif ~isnan(posn1_peak_time) && isnan(negn1_peak_time) % positive peak
                            n1_peak_amplitude       = posn1_peak_amplitude;
                            n1_peak_time            = posn1_peak_time;
                            n1_peak_sample          = posn1_peak_sample;
                        elseif ~isnan(negn1_peak_time) && isnan(posn1_peak_time) % negative peak
                            n1_peak_amplitude       = negn1_peak_amplitude;
                            n1_peak_time            = negn1_peak_time;
                            n1_peak_sample          = negn1_peak_sample;
                        else
                            n1_peak_amplitude       = NaN;
                            n1_peak_time            = NaN;
                            n1_peak_sample          = NaN;
                        end
                        % save counts
                        resp_counter                = resp_counter + 1;
                        % get N1 params
                        out(measure_ind,stim_ind).n1ampl(resp_counter, :) = n1_peak_amplitude;
                        out(measure_ind,stim_ind).n1time(resp_counter, :) = n1_peak_time;
                        out(measure_ind,stim_ind).n1sample(resp_counter, :) = n1_peak_sample;
                        
                        % get CCEP responses for plotting
                        % Scaling to unit length (Euclidean lenght): https://en.wikipedia.org/wiki/Feature_scaling
                        % unit length taken in same window as stats
                        response_vector_length      = sum(plot_responses(all_out(ss).tt > t_win_norm(1) &  all_out(ss).tt < t_win_norm(2)) .^ 2) .^ .5;
                        plot_responses_norm         = plot_responses ./ (response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial
                       
                        out(measure_ind,stim_ind).plot_responses_norm(resp_counter, :) = plot_responses_norm;

                        % save subject index
                        out(measure_ind,stim_ind).subj_ind(resp_counter, :) = ss;
                        
                        % save CRP stuff one layer out 
                        out(measure_ind,stim_ind).p(resp_counter, :)        = all_out(ss).crp_p_adj(these_measured_sites(kk), these_stim_sites(ll));
                        out(measure_ind,stim_ind).cod(resp_counter, :)      = all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll)); 
                        out(measure_ind,stim_ind).a_prime(resp_counter, :)  = all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)); 
                    end
                end 
            end  
        end 
    end 
end


%% plot CCEPs and confidence interval
% Plot pannels in figure 2. Each plot containing the average waveform of 
% limbic connections (all stimulated and all measured sites) across subjects
figure('Position',[0 0 1000 800]), hold on

for measure_ind = 1:length(area_codes)                      % loop over measured sites
    for stim_ind = 1:length(area_codes)                     % loop over stimulated sites
        subplot(length(area_codes), length(area_codes), (measure_ind-1) * length(area_codes) + stim_ind),hold on
        
        sign_resp       = out(measure_ind,stim_ind).p < 0.05; % adjusted for multiple comparisons
        this_set        = out(measure_ind,stim_ind).plot_responses_norm(sign_resp==1,:)';        

%         Plot average waveform and confidence interval (all mean responses
%         -from 10 to 12 trials- across subs optional)
        plot(tt(tt>.015 & tt<1), zeros(size(tt(tt>.015 & tt<1))),'k:')  % plot zero line
%         plot(tt,this_set,'color',[.5 .5 .5 .2])                       % OPTIONAL plot all connections
        plot(tt,mean(this_set,2), 'color','k', 'LineWidth',1)           % plot mean of all connections
        plotCurvConf(tt, this_set',[],'0.5');                           % plots 95% Confidende interval in gray with 50% transparency
        xlim([0 1]),ylim([-0.08 0.08])
        title(['stim:' area_names{stim_ind} ' rec:' area_names{measure_ind} ])
    end
    xlim([0 1]),ylim([-0.08 0.08])
end

%% Get proportion of N1s and significant responses (using the CRPs method)
% Plot bars with total of responses, percentage of N1s (of total responses), 
% and percentage of significant responses based on the CRP method
figure('Position',[0 0 300 800]), hold on

for measure_ind = 1:length(area_codes)
    for stim_ind = 1:length(area_codes)
        subplot(length(area_codes), length(area_codes), (measure_ind-1) * length(area_codes) + stim_ind),hold on
        
        % find all CRP-significant
        sign_resp       = out(measure_ind,stim_ind).p<0.05;         % adjusted for multiple comparisons
        % find all N1s 
        n1s             = ~isnan(out(measure_ind,stim_ind).n1ampl); % found and N1 (pos or neg)
        
        sum_cceps       = length(out(measure_ind,stim_ind).p);      % total nr of CCEPs

        bar(100,'FaceColor',[.9 .9 .9])                             % plot bar representing all responses
        text(.7,sum_cceps,{sum_cceps},'color',[.5 .5 .5],...
            'HorizontalAlignment','center','VerticalAlignment','bottom') 
        
        % Plot bar of significant CRP-based
        sum_CRPs        = sum(sign_resp);
        prop_CRPs       = 100 * sum_CRPs / sum_cceps;
        bar(prop_CRPs,'FaceColor',[0 0.4470 0.7410])                % plot bar with sig CRP
        text(1,prop_CRPs,{prop_CRPs},'color',[0 0.4470 0.7410],...
            'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        % now let's do N1s
        sum_N1s         = sum(n1s);
        prop_N1s        = 100 * sum_N1s / sum_cceps;
        bar(prop_N1s,'FaceColor',[0.92 0.69 0.12])                  % plot bar with N1s
        text(1,prop_N1s,{prop_N1s},'color',[0.92 0.69 0.12],...
            'HorizontalAlignment','left','VerticalAlignment','bottom')
        ylim([0 101])
    end
end

%% Now get proportion of N1s of all CRPs
% Plot bars with total of responses, percentage of significant responses
% based on the CRP method, and percentage of N1s detected from the 
% CRP-significant responses. 
figure('Position',[0 0 300 800]), hold on

for measure_ind = 1:length(area_codes)
    for stim_ind = 1:length(area_codes)
        subplot(length(area_codes), length(area_codes), (measure_ind-1) * length(area_codes) + stim_ind),hold on

        % find all significant responses (CRP-based)
        sign_resp       = out(measure_ind,stim_ind).p<0.05;         % adjusted for multiple comparisons CRP-based
       
        % find all N1s 
        n1s             = ~isnan(out(measure_ind,stim_ind).n1ampl) & sign_resp; % found detected n1 AND significant CRP-based
        bar(100,'FaceColor',[.9 .9 .9])                             % plot bar with all responses?
       
        sum_cceps = length(out(measure_ind,stim_ind).p);            % total nr of CCEPs
        text(.7,sum_cceps,{sum_cceps},'color',[.5 .5 .5],...
            'HorizontalAlignment','center','VerticalAlignment','bottom') 
        
        % Plot bar of significant CRP-based
        sum_CRPs        = sum(sign_resp);
        prop_CRPs       = 100*sum_CRPs/sum_cceps;
        bar(prop_CRPs,'FaceColor',[0 0.4470 0.7410])                % plot bar with sig CRP
        text(1,prop_CRPs,{prop_CRPs},'color',[0 0.4470 0.7410],...
            'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        % now let's do N1s
        sum_N1s         = sum(n1s);
        prop_N1s        = 100*sum_N1s/sum_cceps;
        bar(prop_N1s,'FaceColor',[0.92 0.69 0.12])                  % plot bar with N1s
        text(1,prop_N1s,{prop_N1s},'color',[0.92 0.69 0.12],...
            'HorizontalAlignment','left','VerticalAlignment','bottom')
        ylim([0 101])
    end
end
