
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
area_codes = {[12123 53 54 12108 12109 12110 12106 12107 11123 59 17 18 11108 11109 11110 11106 11107 10]}; % all areas

nr_subs = length(all_subjects);

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


%% Posterior cingulate
rgb_color = {[0.3010 0.7450 0.9330],[1 .8 .1],[0.9290 0.6940 0.1250],[0 0.4470 0.7410]};
rgb_label = {'blue', 'yellow', 'orange', 'light blue'};

area_codes_r = {[12123 53],[54],[12108 12109 12110],[12106 12107],[59]}; % right
area_codes_l = {[11123 17],[18],[11108 11109 11110],[11106 11107],[10]}; % left
area_names = {'Hipp','Amyg','PCC','ACC','ANT'};   

sub_hemi = {'r','r','r','l','r','l','l','r'};

figure('Position',[0 0 500 350]), hold on; %('Position',[0 0 600 200]), hold on

out_plot_responses_norm = [];
out_subj_ind = [];
resp_counter = 0;

ss = 5;

if isequal(sub_hemi{ss},'l')
    area_codes = area_codes_l;
elseif isequal(sub_hemi{ss},'r')
    area_codes = area_codes_r;
end

stim_ind = 1;
measure_ind = 3;

these_measured_sites = find(ismember(all_out(ss).channel_areas,area_codes{measure_ind}));
% Get sites that belong to the stimulated ROI (stim_area)
these_stim_sites = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{stim_ind}) | ...
    ismember(all_out(ss).average_ccep_areas(:,2),area_codes{stim_ind}));

plot(all_out(ss).tt, ss + zeros(size(all_out(ss).tt)),'k')

sign_resp = all_out(ss).crp_p_adj<0.05; % plot p<0.05 FDR corrected

t_win_cod = [0.015 0.500]; % window for vector length normalization and plotting across subjects

% walk through pairs, test significance and plot if so
for kk = 1:length(these_measured_sites)
    for ll = 1:length(these_stim_sites)
        if sign_resp(these_measured_sites(kk),these_stim_sites(ll))==1
            plot_responses = squeeze(all_out(ss).average_ccep(these_measured_sites(kk),these_stim_sites(ll),:));
            plot_responses(all_out(ss).tt>-0.010 & all_out(ss).tt<0.010) = NaN;
            
            % Normalize. Scaling to unit length (Euclidean lenght): https://en.wikipedia.org/wiki/Feature_scaling
            % unit length taken in same window as stats
            response_vector_length = sum(plot_responses(all_out(ss).tt > t_win_cod(1) &  all_out(ss).tt < t_win_cod(2)) .^ 2) .^ .5;
            plot_responses_norm = plot_responses ./ (response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial

            % we looked at the data and visually split out 2 measure
            % sites that produce different signatures
            if kk==1        % waves el1p1, superficial (>2.33mm)
                plot(all_out(ss).tt, ss + plot_responses_norm,'color',rgb_color{1},'LineWidth',.85)
                disp([rgb_label{1} all_out(ss).channel_names(these_measured_sites(kk))])
            elseif kk==2    % waves el2p1, deep (<2.33mm)
                plot(all_out(ss).tt, ss + plot_responses_norm,'color',rgb_color{2},'LineWidth',.85)
                disp([rgb_label{2} all_out(ss).channel_names(these_measured_sites(kk))])
            elseif kk==3    % waves el1p2, deep (<2.33mm)
                plot(all_out(ss).tt, ss + plot_responses_norm,'color',rgb_color{3},'LineWidth',.85)
                disp([rgb_label{3} all_out(ss).channel_names(these_measured_sites(kk))])
            elseif kk==4    % waves el1p3, superficial (>2.33mm)
                plot(all_out(ss).tt, ss + plot_responses_norm,'color',rgb_color{4},'LineWidth',.85)
                disp([rgb_label{4} all_out(ss).channel_names(these_measured_sites(kk))])
            end
            
            % save outputs
            resp_counter = resp_counter+1;
            out_plot_responses_norm(resp_counter,:) = plot_responses_norm;
            out_subj_ind(resp_counter,:) = ss;
        end
    end
end


title('PCC depths','Blue, superficial | Yellow, deep')
xlim([-0.2 .6])%, ylim([-2000 4000]);
xlabel('time (s)')%, ylabel('amplitude (uV)')

set(gcf,'PaperPositionMode','auto')
% print(fullfile('./local',figName),'-dpng');%,'-r700';
% print(fullfile('./local',figName),'-painters','-depsc')%,'-r500',)

