clearvars, clc, close all
% startup

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
    crpFile = fullfile(localDataPath,'derivatives','stats',['sub-' bids_sub],...
            ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_crp_noCAR.mat']);
    [events_table,channels_table,electrodes_table,sub_out] = pcc_loadAveragesStats(localDataPath,bids_sub,bids_ses,bids_task,bids_run,crpFile);
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
sub_color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
% blue, orange, green, purple, mustard, celeste, wine


% changes to the following:
% area_names = {'Hipp','Amyg','PCC','ACC'};   % Left ANT out
% area_codes = {[1 8],[2],[5 6 7],[3 4]};
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

                        % later normalize
                        % Scaling to unit length (Euclidean lenght): https://en.wikipedia.org/wiki/Feature_scaling
                        % unit length taken in same window as stats
                        response_vector_length = sum(plot_responses(all_out(ss).tt > t_win_norm(1) &  all_out(ss).tt < t_win_norm(2)) .^ 2) .^ .5;
                        plot_responses_norm = plot_responses ./ (response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial

                        % save outputs
                        resp_counter = resp_counter + 1;

                        % get responses
                        out(measure_ind,stim_ind).plot_responses_norm(resp_counter, :) = plot_responses_norm;
                        out(measure_ind,stim_ind).plot_responses(resp_counter, :) = plot_responses;
                        out(measure_ind,stim_ind).C(resp_counter, :) = squeeze(all_out(ss).avg_trace_tR(these_measured_sites(kk), these_stim_sites(ll), :));

                        out(measure_ind,stim_ind).subj_ind(resp_counter, :) = ss;

                        out(measure_ind,stim_ind).p(resp_counter, :) = all_out(ss).crp_p_adj(these_measured_sites(kk), these_stim_sites(ll));
                        out(measure_ind,stim_ind).cod(resp_counter, :) = all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll)); 
                        out(measure_ind,stim_ind).a_prime(resp_counter, :) = all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)); 
                    end
                end % done looping through stim pairs
            end  % done looping through measured electrode

        end % done subject loop
    end 
end


%% plot average and std error 
 
tt = all_out(1).tt;

figure('Position',[0 0 1000 800]), hold on
for measure_ind = 1:length(area_codes)
    for stim_ind = 1:length(area_codes)
        subplot(length(area_codes), length(area_codes), (measure_ind-1) * length(area_codes) + stim_ind),hold on

        sign_resp = out(measure_ind,stim_ind).p<0.05; % adjusted for multiple comparisons
%         sign_resp = out(measure_ind,stim_ind).cod>0.30; 

        this_set = out(measure_ind,stim_ind).plot_responses_norm(sign_resp==1,:)';        

%         Plot all response, mean and confidence interval
        plot(tt,this_set,'color',[.5 .5 .5 .2])                   % just plot all responses
        plot(tt,mean(this_set,2), 'color','k', 'LineWidth',1)       % plot mean of all responses
        plot(tt(tt>.015 & tt<1),zeros(size(tt(tt>.015 & tt<1))),'k:') % plot zero line
        plotCurvConf(tt, this_set');                                % plots 95% Confidende interval in gray with 50% transparency

%         % Plot significance along timepoints
%         [h, p] = ttest(this_set');                                  % get ttest & p values
%         hAdj = fdr_bh(p, 0.05, 'pdep');                              % adjust ttest with fdr for any dependency structure ('dep')
%         [M,I] = max(hAdj);
        % show standard deviation on top
%         plot(tt(find(hAdj)), -0.065*hAdj(find(hAdj)), 'b.');
%         ylabel(M)

%        % plot first pc
%         [u,s,v] = svd(this_set(tt>.015 & tt<1,:),'econ');
%         plot(tt(tt>.015 & tt<1),u(:,2)) % PC1 seems to capture limbic-H-wave better, but not perfect


        xlim([0 1]),ylim([-0.08 0.08])
        title(['stim:' area_names{stim_ind} ' rec:' area_names{measure_ind} ])
    end
    xlim([0 1]),ylim([-0.08 0.08])
end


% xlim([-0.2 .6]), ylim([10 150]);
% xlabel('Time (s)'), ylabel('Normalized Amplitude')

% Saving .png and .eps
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300','-vector',fullfile(localDataPath,'derivatives','figures',...
    ['allsub_CAR']))
print('-depsc2','-r300','-vector',fullfile(localDataPath,'derivatives','figures',...
    ['allsub_CAR']))



