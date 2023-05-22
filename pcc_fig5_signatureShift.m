clearvars, close all, clc
% startup

%% Plot figure showing how positive and negative peaks come from rec electrodes in different depths
% Dependencies: MNL ieeg basics repository
% cd to HAPwave repository

addpath(genpath(pwd))

% set local path to your BIDS directory:
myPath = setLocalDataPath(1);
localDataPath = myPath.input;

% load the meta data
all_subjects = {'01','02','03','04','05','06','07'}; % S4 & S5 have temporal delay
all_hemi = {'r','r','r','l','r','l','l'};
all_runs = {'01','01','01','01','01','01','01',};

nr_subs = length(all_subjects);
% t_win_cod = [0.015 .5];

for ss = 1:nr_subs

    bids_sub = all_subjects{ss};
    bids_ses = 'ieeg01';
    bids_task = 'ccep';
    bids_run = all_runs{ss};

    [events_table,channels_table,electrodes_table,sub_out] = pcc_loadAveragesStats(localDataPath,bids_sub,bids_ses,bids_task,bids_run);
    all_out(ss) = sub_out;

end

%% Correct P values for number of comparisons in each subject

area_codes = {[12123 53 54 12108 12109 12110 12106 12107 11123 59 17 18 11108 11109 11110 11106 11107 10]}; % all areas

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
area_names = {'Hipp','Amyg','PCC','ACC','ANT'};   
area_codes_r = {[12123 53],[54],[12108 12109 12110],[12106 12107],[59]}; % right
area_codes_l = {[11123 17],[18],[11108 11109 11110],[11106 11107],[10]}; % left

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

                        % save outputs
                        resp_counter = resp_counter + 1;

                        % get CCEP responses for plotting
                        % Scaling to unit length (Euclidean lenght): https://en.wikipedia.org/wiki/Feature_scaling
                        % unit length taken in same window as stats
                        response_vector_length = sum(plot_responses(all_out(ss).tt > t_win_norm(1) &  all_out(ss).tt < t_win_norm(2)) .^ 2) .^ .5;
                        plot_responses_norm = plot_responses ./ (response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial
                        out(measure_ind,stim_ind).plot_responses_norm(resp_counter, :) = plot_responses_norm;
                        out(measure_ind,stim_ind).plot_responses(resp_counter, :) = plot_responses;

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




%% Plot polarity shift by recording site
% rgb_color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
% blue, orange, green, purple, mustard, celeste, wine


% stim_area = 1; % 1 for hippocampus, 2 for amygdala, 3 for thalamus
% area_names = {'Hipp','Amyg','PCC','ACC',''};   
measure_ind = 3;
stim_ind = 1;

these_measured_sites = area_codes{measure_ind};
these_stim_sites = area_codes{stim_ind};

out_plot_responses_norm = [];
out_subj_ind = [];
resp_counter = 0;

figure('Position',[0 0 500 350]), hold on; %('Position',[0 0 600 200]), hold on

for ss = 4 %1:7
% for ss = bids_sub

    
%     if stim_area==1
%         these_stim_sites = find(all_out(ss).stim_area==1 | all_out(ss).stim_area==8); % hippocampal formation 
%     elseif stim_area==2
%         these_stim_sites = find(all_out(ss).stim_area==2); % amygdala
%     elseif stim_area==3
%         these_stim_sites = find(all_out(ss).stim_area==9); % thalamus
%     end
%   
    plot(all_out(ss).tt, ss + zeros(size(all_out(ss).tt)),'k')
    
    % only significant responses (FDR corrected)
    sign_resp = out(measure_ind,stim_ind).p<0.05; % adjusted for multiple comparisons
    
    % walk through pairs, test significance and plot if so
    for kk = 1:length(these_measured_sites)
        for ll = 1:length(these_stim_sites)
            if sign_resp(these_measured_sites(kk),these_stim_sites(ll))==1
% %                 plot_responses = squeeze(all_out(ss).average_ccep(these_measured_sites(kk),these_stim_sites(ll),:));
%               plot_responses = out(these_measured_sites(kk),these_stim_sites(ll)).plot_responses(sign_resp,:)';
%               plot_responses(all_out(ss).tt>-0.010 & all_out(ss).tt<0.010) = NaN;
%                 
%                 % Normalize. Scaling to unit length (Euclidean lenght): https://en.wikipedia.org/wiki/Feature_scaling
%                 % unit length taken in same window as stats
%                 response_vector_length = sum(plot_responses(all_out(ss).tt > t_win_cod(1) &  all_out(ss).tt < t_win_cod(2)) .^ 2) .^ .5;
%                 plot_responses_norm = plot_responses ./ (response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial

                % we looked at the data and visually split out 2 measure
                % sites that produce different signatures
                if kk==1 % measurement el 1 (RY1/MPC)
                    plot(all_out(ss).tt, ss + plot_responses_norm,'color',[0 0.4470 0.7410 ],'LineWidth',.85)
                    disp(['blue ' all_out(ss).channel_names(these_measured_sites(kk))])
                elseif kk==2 % measurement el 2 (RY2/MPC)
                    plot(all_out(ss).tt, ss + plot_responses_norm,'color',[1 .8 .1],'LineWidth',.85)
                    disp(['yellow ' all_out(ss).channel_names(these_measured_sites(kk))])
                elseif kk==3 % measurement el 3 ((RZ1/PDC)
                    plot(all_out(ss).tt, ss + plot_responses_norm,'color',[.8 .8 .8],'LineWidth',.85)
                    disp(['gray ' all_out(ss).channel_names(these_measured_sites(kk))])
%                 elseif kk==4 % measurement el 3 ((RZ1/PDC)
%                     plot(all_out(ss).tt, ss + plot_responses_norm,'color',[1 .8 .1],'LineWidth',.85)
%                     disp(['yellow ' all_out(ss).channel_names(these_measured_sites(kk))])
                end
                
                % save outputs
                resp_counter = resp_counter+1;
                out_plot_responses_norm(resp_counter,:) = plot_responses_norm;
                out_subj_ind(resp_counter,:) = ss;
            end
        end
    end
end

title('Blue=RY1;Yellow=RY2')
xlim([-0.2 .6])%, ylim([-2000 4000]);
xlabel('time (s)')%, ylabel('amplitude (uV)')

set(gcf,'PaperPositionMode','auto')
% print(fullfile('./local',figName),'-dpng');%,'-r700';
% print(fullfile('./local',figName),'-painters','-depsc')%,'-r500',)

