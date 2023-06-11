clearvars, clc, close all
% startup 

%% Load stats across subs
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

area_codes = {[12123 53 54 12108 12109 12110 12106 12107 11123 49 17 18 11108 11109 11110 11106 11107 10]}; % all areas

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
    [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals,qq,'dep','no');
    all_out(ss).crp_p_adj(all_out(ss).hasdata==1) = adj_p;
    all_out(ss).h(all_out(ss).hasdata==1) = h;
end

%% Load stim and measure areas across subjects
set_color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
% blue, orange, green, purple, mustard, celeste, wine

area_names = {'Hipp','Amyg','PCC','ACC','ANT'};   
area_codes_r = {[12123 53],[54],[12108 12109 12110],[12106 12107],[49]}; % right
area_codes_l = {[11123 17],[18],[11108 11109 11110],[11106 11107],[10]}; % left

sub_hemi = {'r','r','r','l','r','l','l','r'};

nr_subs = length(sub_hemi);

out = []; % this will be a area X area structure, with all subjects concatinates for each area

subj_resp_total = zeros(nr_subs,1);               % stim-->measured pair for stats FDR correction

t_win_norm = [0.015 0.500]; % window for vector length normalization and plotting across subjects

for measure_ind = 1:length(area_names) % loop through the inds 
    for stim_ind = 1:length(area_names)

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
                        
                        % store subject index
                        out(measure_ind,stim_ind).subj_ind(resp_counter, :) = ss;
                        
                        % save CRP stuff for linear mixed-effects model
                        out(measure_ind,stim_ind).p(resp_counter, :) = all_out(ss).crp_p_adj(these_measured_sites(kk), these_stim_sites(ll));
                        out(measure_ind,stim_ind).cod(resp_counter, :) = all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll)); 
                        out(measure_ind,stim_ind).a_prime(resp_counter, :) = all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)); 
                    end
                end % done looping through stim pairs
            end  % done looping through measured electrode

        end % done subject loop
    end 
end

%% General linear fixed effects model

% area_names = {'Hipp','Amyg','PCC','ACC'};   
% out(measure_ind,stim_ind)

pcc_hc_amp = out(3,1).cod;
pcc_hc_stimgroup = cell(size(pcc_hc_amp)); 
pcc_hc_measgroup = cell(size(pcc_hc_amp)); 
[pcc_hc_stimgroup{:}] = deal('HC');
[pcc_hc_measgroup{:}] = deal('PCC');
sub_ind1 = out(3,1).subj_ind;

acc_hc_amp = out(4,1).cod;
acc_hc_stimgroup = cell(size(acc_hc_amp)); % hc stim gets code 1
acc_hc_measgroup = cell(size(acc_hc_amp)); % acc measure gets code 4
[acc_hc_stimgroup{:}] = deal('HC');
[acc_hc_measgroup{:}] = deal('ACC');
sub_ind2 = out(4,1).subj_ind;

pcc_amg_amp = out(3,2).cod;
pcc_amg_stimgroup = cell(size(pcc_amg_amp)); % amg stim gets code 2
pcc_amg_measgroup = cell(size(pcc_amg_amp)); % pcc measure gets code 3
[pcc_amg_stimgroup{:}] = deal('AMG');
[pcc_amg_measgroup{:}] = deal('PCC');
sub_ind3 = out(3,2).subj_ind;

acc_amg_amp = out(4,2).cod;
acc_amg_stimgroup = cell(size(acc_amg_amp)); % amg stim gets code 2
acc_amg_measgroup = cell(size(acc_amg_amp)); % acc measure gets code 4
[acc_amg_stimgroup{:}] = deal('AMG');
[acc_amg_measgroup{:}] = deal('ACC');
sub_ind4 = out(4,2).subj_ind;

sub_ind = [sub_ind1; sub_ind2; sub_ind3; sub_ind4];
y = [pcc_hc_amp(:); acc_hc_amp(:); pcc_amg_amp(:); acc_amg_amp(:)];
stim_site = [pcc_hc_stimgroup(:); acc_hc_stimgroup(:); pcc_amg_stimgroup(:); acc_amg_stimgroup(:)]; % stim group
measure_site = [pcc_hc_measgroup(:); acc_hc_measgroup(:); pcc_amg_measgroup(:); acc_amg_measgroup(:)]; % measurement group

tbl = table(y,categorical(stim_site),categorical(measure_site),categorical(sub_ind),...
    'VariableNames',{'ccep_val','stim_site','measure_site','sub_ind'});

lme = fitlme(tbl,'ccep_val ~ 1 + stim_site*measure_site + (1|sub_ind)') % with interaction

% now test only for PCC
sub_ind = [sub_ind1; sub_ind3];
y = [pcc_hc_amp(:); pcc_amg_amp(:)];
stim_site = [pcc_hc_stimgroup(:); pcc_amg_stimgroup(:)]; % stim group
tbl = table(y,categorical(stim_site),categorical(sub_ind),...
    'VariableNames',{'ccep_val','stim_site','sub_ind'});
lme = fitlme(tbl,'ccep_val ~ 1 + stim_site + (1|sub_ind)') % with interaction

% now test only for ACC
sub_ind = [sub_ind2; sub_ind4];
y = [acc_hc_amp(:); acc_amg_amp(:)];
stim_site = [acc_hc_stimgroup(:); acc_amg_stimgroup(:)]; % stim group
tbl = table(y,categorical(stim_site),categorical(sub_ind),...
    'VariableNames',{'ccep_val','stim_site','sub_ind'});
lme = fitlme(tbl,'ccep_val ~ 1 + stim_site + (1|sub_ind)') % with interaction

%% Distribution plot

sub_ind = [sub_ind1; sub_ind2; sub_ind3; sub_ind4];
y = [pcc_hc_amp(:); acc_hc_amp(:); pcc_amg_amp(:); acc_amg_amp(:)];
stim_site = [pcc_hc_stimgroup(:); acc_hc_stimgroup(:); pcc_amg_stimgroup(:); acc_amg_stimgroup(:)]; % stim group
measure_site = [pcc_hc_measgroup(:); acc_hc_measgroup(:); pcc_amg_measgroup(:); acc_amg_measgroup(:)]; % measurement group

plot_groups = zeros(size(stim_site));
plot_groups(ismember(stim_site,'HC') & ismember(measure_site,'ACC')) = 1;
plot_groups(ismember(stim_site,'AMG') & ismember(measure_site,'ACC')) = 2;
plot_groups(ismember(stim_site,'HC') & ismember(measure_site,'PCC')) = 3;
plot_groups(ismember(stim_site,'AMG') & ismember(measure_site,'PCC')) = 4;

% X = [X; NaN; NaN; NaN; NaN];
% plot_groups = [plot_groups; 1; 2; 3; 4];

% figure,hold on
figure('Position',[0 0 250 400]), hold on;
distributionPlot(y,'groups',plot_groups,'color',[.7 .7 1],'distWidth',1,'showMM',0);%

% plotSpread
plotSpread(y, 'distributionIdx', plot_groups, 'xValues', 1:4, 'spreadWidth', .5,...
    'distributionColors', [0.5, 0.5, 0.5]);

plot([1:4],[median(y(plot_groups==1)) median(y(plot_groups==2)) median(y(plot_groups==3)) median(y(plot_groups==4))],'k.','MarkerSize',30);

xlim([0 5]),ylim([-0.6 1.4]),ylabel('Explained Variance (R^2)')
set(gca,'YTick',[0:.2:1],'XTick',[1:4],'XTickLabel',{'H-PCC','A-PCC','H-ACC','A-ACC'})


