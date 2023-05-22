clearvars, clc, close all
addpath(genpath(pwd))

% set local path to your BIDS directory:
myPath = setLocalDataPath(1);
localDataPath = myPath.input;

% load the meta data
all_subjects = {'01','02','03','04','05','06','07','08'}; % 
all_hemi = {'r','r','r','l','r','l','l','r'};
all_runs = {'01','01','01','01','01','01','01','01'};

% load the meta data
for ss = 1:length(all_subjects)
    bids_sub = all_subjects{ss};
    bids_ses = 'ieeg01';
    bids_task = 'ccep';
    bids_run = all_runs{ss};

    [events_table,channels_table,electrodes_table,sub_out] = pcc_loadAveragesStats(localDataPath,bids_sub,bids_ses,bids_task,bids_run);
    all_out(ss) = sub_out;
end


%% Correct P values for number of comparisons in each subject

% add thalamic sites to this list
area_codes = {[12123 53 54 12108 12109 12110 12106 12107 11123 17 18 11108 11109 11110 11106 11107]}; % all areas

nr_subs = length(all_subjects);

out = []; % this will be a area X area structure, with all subjects concatinates for each area

subj_resp_total = zeros(nr_subs,1);               % stim-->measured pair for stats FDR correction


resp_counter = 0; % counting all responses across subjects for this connection

for ss = 1:nr_subs % subject loop
   
    % Get sites that belong to the measurement ROI (rec_area)
    these_measured_sites = find(ismember(all_out(ss).channel_areas,area_codes{1}));
    
    % Get sites that belong to the stimulated ROI (stim_area)
    these_stim_sites = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{1}) | ...
        ismember(all_out(ss).average_ccep_areas(:,2),area_codes{1}));
    
    % p-values for correction of multiple comparisons
    all_out(ss).hasdata = NaN(size(all_out(ss).crp_out));
    all_out(ss).crp_p = NaN(size(all_out(ss).crp_out));
    all_out(ss).a_prime = NaN(size(all_out(ss).crp_out));
    all_out(ss).cod = NaN(size(all_out(ss).crp_out));
    all_out(ss).crp_p_adj = NaN(size(all_out(ss).crp_out));
    all_out(ss).h = NaN(size(all_out(ss).crp_out));
    all_out(ss).avg_trace_tR = zeros(size(all_out(ss).average_ccep));
    
    % loop over measured sites
    for kk = 1:length(these_measured_sites)
        % loop over the stimulated pairs
        for ll = 1:length(these_stim_sites)
            if ~isempty(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).data) % ~ same stim/recording site
                all_out(ss).hasdata(these_measured_sites(kk), these_stim_sites(ll)) = 1;
                all_out(ss).crp_p(these_measured_sites(kk), these_stim_sites(ll)) = all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_projs.p_value_tR;
                all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)) = mean(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.al_p); % mean alpha prime across trials;
                all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll)) = median(all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.cod); % median across trials
                sig_timepoints = find(all_out(ss).tt>0.015 & all_out(ss).tt<=all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.tR);
                all_out(ss).avg_trace_tR(these_measured_sites(kk), these_stim_sites(ll),sig_timepoints) = all_out(ss).crp_out(these_measured_sites(kk), these_stim_sites(ll)).crp_parms.C;

            else
                all_out(ss).hasdata(these_measured_sites(kk), these_stim_sites(ll)) = 0;
            end
        end
    end
    pvals = all_out(ss).crp_p(all_out(ss).hasdata==1);
    qq = 0.05;
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,qq,'pdep','no');
    all_out(ss).crp_p_adj(all_out(ss).hasdata==1) = adj_p;
    all_out(ss).h(all_out(ss).hasdata==1) = h;
end

%% Load stim and measure areas across subjects
set_color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
% blue, orange, green, purple, mustard, celeste, wine

% changes to the following:
area_names = {'Hipp','Amyg','PCC','ACC'};   
area_codes_r = {[12123 53],[54],[12108 12109 12110],[12106 12107]}; % right
area_codes_l = {[11123 17],[18],[11108 11109 11110],[11106 11107]}; % left

sub_hemi = {'r','r','r','l','r','l','l','r'};

nr_subs = length(sub_hemi);

out = []; % this will be a area X area structure, with all subjects concatinates for each area

subj_resp_total = zeros(nr_subs,1);               % stim-->measured pair for stats FDR correction

t_win_norm = [0.015 0.500]; % window for vector length normalization and plotting across subjects

for measure_ind = 1:length(area_codes_r) % loop through the inds 
    for stim_ind = 1:length(area_codes_r)

        resp_counter = 0; % counting all responses across subjects for this connection

        for ss = 1:nr_subs % subject loop
            if isequal(sub_hemi{ss},'l')
                area_codes = area_codes_l;
            elseif isequal(sub_hemi{ss},'r')
                area_codes = area_codes_r;
            end

            % Get sites that belong to the measurement ROI (rec_area)
            these_measured_sites = find(ismember(all_out(ss).channel_areas,area_codes{measure_ind}));
            
            % Get sites that belong to the stimulated ROI (stim_area)
            these_stim_sites = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{stim_ind}) | ...
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
                        plot_responses = squeeze(all_out(ss).average_ccep(these_measured_sites(kk), these_stim_sites(ll), :));
                        tt = all_out(ss).tt;

                        % Is there an N1?
                        params.amplitude_thresh = 3.4;% standard deviations 3.4;
                        params.n1_peak_range = 0.050; %in s
                        params.baseline_tt = tt>-.5 & tt<-.020;
                        params.peakSign = 1;
                        % amplitudeThresh is in microvolts, used if standard devisation from baseline is smaller than this value 
                        % threshold for responses is (amplitude_thresh * amplitudeThresh) is baseline STD < amplitudeThresh
                        % threshold for responses is (amplitude_thresh * baseline STD) is baseline STD > amplitudeThresh
                        params.amplitudeThresh = 10;

                        [posn1_peak_sample,posn1_peak_amplitude,posn1_peak_time] = ccep_detect_n1peak_sEEG(plot_responses,tt,params);
                        params.peakSign = -1;
                        [negn1_peak_sample,negn1_peak_amplitude,negn1_peak_time] = ccep_detect_n1peak_sEEG(plot_responses,tt,params);
                        
                        % Keep earliest peak
                        if ~isnan(posn1_peak_time) && ~isnan(negn1_peak_time) % if both positive and negative
                            if posn1_peak_time < negn1_peak_time % positive is earlier
                                n1_peak_amplitude = posn1_peak_amplitude;
                                n1_peak_time = posn1_peak_time;
                                n1_peak_sample = posn1_peak_sample;
                            elseif negn1_peak_time < posn1_peak_time
                                n1_peak_amplitude = negn1_peak_amplitude;
                                n1_peak_time = negn1_peak_time;
                                n1_peak_sample = negn1_peak_sample;
                            end
                        elseif ~isnan(posn1_peak_time) && isnan(negn1_peak_time)
                            n1_peak_amplitude = posn1_peak_amplitude;
                            n1_peak_time = posn1_peak_time;
                            n1_peak_sample = posn1_peak_sample;
                        elseif ~isnan(negn1_peak_time) && isnan(posn1_peak_time)
                            n1_peak_amplitude = negn1_peak_amplitude;
                            n1_peak_time = negn1_peak_time;
                            n1_peak_sample = negn1_peak_sample;
                        else
                            n1_peak_amplitude = NaN;
                            n1_peak_time = NaN;
                            n1_peak_sample = NaN;
                        end
                        % save outputs
                        resp_counter = resp_counter + 1;

                        % get N1 responses
                        out(measure_ind,stim_ind).n1ampl(resp_counter, :) = n1_peak_amplitude;
                        out(measure_ind,stim_ind).n1time(resp_counter, :) = n1_peak_time;
                        out(measure_ind,stim_ind).n1sample(resp_counter, :) = n1_peak_sample;
                        
                        % get CCEP responses for plotting
                        % Scaling to unit length (Euclidean lenght): https://en.wikipedia.org/wiki/Feature_scaling
                        % unit length taken in same window as stats
                        response_vector_length = sum(plot_responses(all_out(ss).tt > t_win_norm(1) &  all_out(ss).tt < t_win_norm(2)) .^ 2) .^ .5;
                        plot_responses_norm = plot_responses ./ (response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial
                        out(measure_ind,stim_ind).plot_responses_norm(resp_counter, :) = plot_responses_norm;
                        
                        % store subject index
                        out(measure_ind,stim_ind).subj_ind(resp_counter, :) = ss;
                        
                        % save CRP stuff for 
                        out(measure_ind,stim_ind).p(resp_counter, :) = all_out(ss).crp_p_adj(these_measured_sites(kk), these_stim_sites(ll));
                        out(measure_ind,stim_ind).cod(resp_counter, :) = all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll)); 
                        out(measure_ind,stim_ind).a_prime(resp_counter, :) = all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)); 
                    end
                end % done looping through stim pairs
            end  % done looping through measured electrode

        end % done subject loop
    end 
end


%% plot CCEPs with N1 in dots and std error 
 
tt = all_out(1).tt;

figure('Position',[0 0 1000 800]), hold on
for measure_ind = 1:length(area_codes)
    for stim_ind = 1:length(area_codes)
        subplot(length(area_codes), length(area_codes), (measure_ind-1) * length(area_codes) + stim_ind),hold on
        
        sign_resp = out(measure_ind,stim_ind).p<0.05; % adjusted for multiple comparisons
       
        this_set = out(measure_ind,stim_ind).plot_responses_norm(sign_resp==1,:)';        

%         Plot all response, mean and confidence interval
%         plot(tt,this_set,'color',[.5 .5 .5 .2])                   % just plot all responses
        plot(tt,mean(this_set,2), 'color','k', 'LineWidth',1)       % plot mean of all responses
        plot(tt(tt>.015 & tt<1),zeros(size(tt(tt>.015 & tt<1))),'k:') % plot zero line
        plotCurvConf(tt, this_set');                                % plots 95% Confidende interval in gray with 50% transparency
        
        % if a significant N1, plot it
        for kk = 1:length(sign_resp) % adjusted for multiple comparisons
            if sign_resp(kk)==1 % significant CRP
%                 if ~isnan(out(measure_ind,stim_ind).n1ampl(kk)) % found and N1 (pos or neg)
%                     plot(tt(out(measure_ind,stim_ind).n1sample(kk)),... 
%                         out(measure_ind,stim_ind).plot_responses_norm(kk,out(measure_ind,stim_ind).n1sample(kk)),...
%                         'r.')
%                 end
            end
        end
              
        xlim([0 1]),ylim([-0.08 0.08])
        title(['stim:' area_names{stim_ind} ' rec:' area_names{measure_ind} ])
    end
    xlim([0 1]),ylim([-0.08 0.08])
end

% xlim([-0.2 .6]), ylim([10 150]);
% xlabel('Time (s)'), ylabel('Normalized Amplitude')

%% now get proportion of N1s of all CRPs

figure('Position',[0 0 300 800]), hold on

for measure_ind = 1:length(area_codes)
    for stim_ind = 1:length(area_codes)
        subplot(length(area_codes), length(area_codes), (measure_ind-1) * length(area_codes) + stim_ind),hold on

        % all significant responses 
        sign_resp = out(measure_ind,stim_ind).p<0.05; % adjusted for multiple comparisons
       
        % all N1s 
        n1s = ~isnan(out(measure_ind,stim_ind).n1ampl); % found and N1 (pos or neg)
        
%         bar(length(out(measure_ind,stim_ind).p),'FaceColor',[.9 .9 .9]) % plot bar with all responses?
        bar(100,'FaceColor',[.9 .9 .9]) % plot bar with all responses?
        
        sum_cceps = length(out(measure_ind,stim_ind).p); % total nr of CCEPs
        text(.7,sum_cceps,{sum_cceps},'color',[.5 .5 .5],...
            'HorizontalAlignment','center','VerticalAlignment','bottom') 
        
        bar(100*sum(sign_resp)/sum_cceps,'FaceColor',[0 0.4470 0.7410])         % plot bar with sig CRP
        sum_CRPs = sum(sign_resp);
        text(1,sum_CRPs,{sum_CRPs},'color',[0 0.4470 0.7410],...
            'HorizontalAlignment','center','VerticalAlignment','bottom')
        
        bar(100*sum(n1s)/sum_cceps,'FaceColor',[0.92 0.69 0.12])               % plot bar with N1s
        sum_N1s = sum(n1s);
        text(1,sum_N1s,{sum_N1s},'color',[0.92 0.69 0.12],...
            'HorizontalAlignment','left','VerticalAlignment','bottom')
        
        ylim([0 101])
    end
end
