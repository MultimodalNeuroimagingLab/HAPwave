clearvars, clc, close all

%% Plot interaction effect & LME model (fig 4.f)
% Dependencies: mnl_ieeg_basics, matmef, and vistasoft github repositories. 
% cd to HAPwave repo
addpath(genpath(pwd)) 

% set local path to your BIDS directory:
myPath          = setLocalDataPath(1);
localDataPath   = myPath.input;

% load the metadata
all_subjects    = {'01','02','03','04','05','06','07','08'};   % List of subjects
all_hemi        = {'r','r','r','l','r','l','l','r'};               % List of hemispheres
all_runs        = {'01','01','01','01','01','01','01','01'};       % List of runs

for ss              = 1:length(all_subjects)
    bids_sub        = all_subjects{ss};
    bids_ses        = 'ieeg01';
    bids_task       = 'ccep';
    bids_run        = all_runs{ss};

    % Load metadata and stats
    [events_table,channels_table,electrodes_table,sub_out] = pcc_loadAveragesStats(localDataPath,bids_sub,bids_ses,bids_task,bids_run);
    all_out(ss)     = sub_out;
end

%% Correct P values for number of comparisons in each subject
% Fisrt load the limbic codes in both right and left hemis
area_codes  = {[12123 53 54 12108 12109 12110 12106 12107 59 11123 17 18 11108 11109 11110 11106 11107 10]}; % all areas
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
    qq                                              = 0.05;     % false discovery rate
    [h, crit_p, adj_ci_cvrg, adj_p]                 = fdr_bh(pvals,qq,'dep','no'); % Benjamini & Yekutieli FDR correction 
    all_out(ss).crp_p_adj(all_out(ss).hasdata==1)   = adj_p;    % Adjusted pVals
    all_out(ss).h(all_out(ss).hasdata==1)           = h;        % adjusted pVal is significant
end


%% Load stim and measure areas across subjects
% Sort limbic codes by hemisphere
area_names      = {'Hipp','Amyg','PCC','ACC'};   
area_codes_r    = {[12123 53],[54],[12108 12109 12110],[12106 12107]}; % right
area_codes_l    = {[11123 17],[18],[11108 11109 11110],[11106 11107]}; % left
nr_subs         = length(all_subjects); % number of subjects

out             = [];                   % prepare areaXarea matrix to concatinate across subs

subj_resp_total = zeros(nr_subs,1);     % stim-->measured pair for adjusted FDR

t_win_norm      = [0.015 0.500];        % window for vector length normalization and plotting across subjects

for measure_ind = 1:length(area_names)  % loop over measured sites
    for stim_ind = 1:length(area_names) % now go through stimulated sites
        resp_counter = 0;               % counting all responses across subjects for this connection

        for ss = 1:nr_subs              % loop over subjects
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
                        subj_resp_total(ss) = subj_resp_total(ss) + 1; % set counter
                        
                        % first raw responses
                        plot_responses      = squeeze(all_out(ss).average_ccep(these_measured_sites(kk), these_stim_sites(ll), :));

                        % save outputs
                        resp_counter        = resp_counter + 1;

                        % get CCEP responses for plotting
                        % Scaling to unit length (Euclidean lenght): https://en.wikipedia.org/wiki/Feature_scaling
                        % unit length taken in same window as stats
                        response_vector_length  = sum(plot_responses(all_out(ss).tt > t_win_norm(1) &  all_out(ss).tt < t_win_norm(2)) .^ 2) .^ .5;
                        plot_responses_norm     = plot_responses ./ ...
                            (response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial
                        out(measure_ind,stim_ind).plot_responses_norm(resp_counter, :) = plot_responses_norm;
                        
                        % store subject index
                        out(measure_ind,stim_ind).subj_ind(resp_counter, :) = ss;
                        
                        % save CRP parms for linear mixed-effects model
                        out(measure_ind,stim_ind).p(resp_counter, :)        = all_out(ss).crp_p_adj(these_measured_sites(kk), these_stim_sites(ll));
                        out(measure_ind,stim_ind).cod(resp_counter, :)      = all_out(ss).cod(these_measured_sites(kk), these_stim_sites(ll)); 
                        out(measure_ind,stim_ind).a_prime(resp_counter, :)  = all_out(ss).a_prime(these_measured_sites(kk), these_stim_sites(ll)); 
                    end
                end
            end
        end
    end 
end

%% General linear mixed effects model
% Stim and measured ROI:
% HC= 1; Amg= 2; PCC= 3; ACC= 4
pcc_hc_cod              = out(3,1).cod;            % measured PCC (3) and stim HC (1)
pcc_hc_stimgroup        = cell(size(pcc_hc_cod));  % hc stim gets code 1 
pcc_hc_measgroup        = cell(size(pcc_hc_cod));  % pcc measure gets code 3
[pcc_hc_stimgroup{:}]   = deal('HC');
[pcc_hc_measgroup{:}]   = deal('PCC');
sub_ind1                = out(3,1).subj_ind;        % HC-PCC connection 

acc_hc_cod              = out(4,1).cod;             % measured ACC (4) and stim HC (1)
acc_hc_stimgroup        = cell(size(acc_hc_cod));   % hc stim gets code 1
acc_hc_measgroup        = cell(size(acc_hc_cod));   % acc measure gets code 4
[acc_hc_stimgroup{:}]   = deal('HC');
[acc_hc_measgroup{:}]   = deal('ACC');
sub_ind2                = out(4,1).subj_ind;        % HC-ACC connection 

pcc_amg_cod             = out(3,2).cod;             % measured PCC (3) and stim Amygdala (2)
pcc_amg_stimgroup       = cell(size(pcc_amg_cod));  % amg stim gets code 2
pcc_amg_measgroup       = cell(size(pcc_amg_cod));  % pcc measure gets code 3
[pcc_amg_stimgroup{:}]  = deal('AMG');
[pcc_amg_measgroup{:}]  = deal('PCC');       
sub_ind3                = out(3,2).subj_ind;        % Amg-PCC connection

acc_amg_cod             = out(4,2).cod;             % measured ACC (4) and stim Amygdala (2)
acc_amg_stimgroup       = cell(size(acc_amg_cod));  % amg stim gets code 2
acc_amg_measgroup       = cell(size(acc_amg_cod));  % acc measure gets code 4
[acc_amg_stimgroup{:}]  = deal('AMG');
[acc_amg_measgroup{:}]  = deal('ACC');
sub_ind4                = out(4,2).subj_ind;        % Amg-ACC connection 

% prepare to fit LME model
sub_ind         = [sub_ind1; sub_ind2; sub_ind3; sub_ind4];                                                 % take the four connections
y               = [pcc_hc_cod(:); acc_hc_cod(:); pcc_amg_cod(:); acc_amg_cod(:)];                           % prepare matrix with cod
stim_site       = [pcc_hc_stimgroup(:); acc_hc_stimgroup(:); pcc_amg_stimgroup(:); acc_amg_stimgroup(:)];   % grup by stim ROI
measure_site    = [pcc_hc_measgroup(:); acc_hc_measgroup(:); pcc_amg_measgroup(:); acc_amg_measgroup(:)];   % group by measured ROI

tbl = table(y,categorical(stim_site),categorical(measure_site),categorical(sub_ind),...
    'VariableNames',{'ccep_val','stim_site','measure_site','sub_ind'});
%% test different LME models for best fit
% Random intercept model with a fixed slope
lme = fitlme(tbl,'ccep_val ~ 1 + stim_site*measure_site + (1|sub_ind)');
anova(lme);                             % perform F-test that all fixed-effects coefficients are zero
anova(lme,'DFMethod','satterthwaite')   % get degrees of freedom w Satterthwaite method. Produces smaller denominator degrees of freedom and slightly larger p-values
disp('Random intercept model with a fixed slope');

% Random intercept and slope, with possible correlation between them
lme2 = fitlme(tbl,'ccep_val ~ 1 + stim_site*measure_site + (1 + stim_site*measure_site|sub_ind)');
anova(lme2);                            % perform F-test that all fixed-effects coefficients are zero
anova(lme2,'DFMethod','satterthwaite')  % get degrees of freedom w Satterthwaite method. Produces smaller denominator degrees of freedom and slightly larger p-values
disp('Random intercept and slope, with possible correlation between them');

% Independent random effects terms for intercept and slope
lme3 = fitlme(tbl,'ccep_val ~ 1 + stim_site*measure_site + (1|sub_ind) + (-1 + stim_site*measure_site|sub_ind)');
anova(lme3);                            % perform F-test that all fixed-effects coefficients are zero
anova(lme3,'DFMethod','satterthwaite')  % get degrees of freedom w Satterthwaite method. Produces smaller denominator deg
disp('Independent random effects terms for intercept and slope')

%% Compare LME models
% %% Run a Theoretical Likelihood Ratio Test to choose the model to use
% %% Check nested conditions between models with 1000 replications     
compare(lme, lme2,'CheckNesting',true,'NSim',1000)  % compare smaller with larger model
compare(lme3,lme2,'CheckNesting',true,'NSim',1000)  % compare complex models
compare(lme, lme3,'CheckNesting',true,'NSim',1000)  % compare smaller and large model

%% Test LME further for stim ROI
% now test only for PCC
sub_ind     = [sub_ind1; sub_ind3];
y           = [pcc_hc_cod(:); pcc_amg_cod(:)];
stim_site   = [pcc_hc_stimgroup(:); pcc_amg_stimgroup(:)]; % stim group
tbl         = table(y,categorical(stim_site),categorical(sub_ind),'VariableNames',{'ccep_val','stim_site','sub_ind'});

lme_pcc     = fitlme(tbl,'ccep_val ~ 1 + stim_site + (stim_site|sub_ind)');
anova(lme_pcc);                           % perform F-test that all fixed-effects coefficients are zero
anova(lme_pcc,'DFMethod','satterthwaite') % degrees of freedom w Satterthwaite method
disp('Stim site interaction in PCC ')

% now test only for ACC
sub_ind     = [sub_ind2; sub_ind4];
y           = [acc_hc_cod(:); acc_amg_cod(:)];
stim_site   = [acc_hc_stimgroup(:); acc_amg_stimgroup(:)]; % stim group
tbl         = table(y,categorical(stim_site),categorical(sub_ind),'VariableNames',{'ccep_val','stim_site','sub_ind'});

lme_acc     = fitlme(tbl,'ccep_val ~ 1 + stim_site + (stim_site|sub_ind)'); 
anova(lme_acc);                           % perform F-test that all fixed-effects coefficients are zero
anova(lme_acc,'DFMethod','satterthwaite') % degrees of freedom w Satterthwaite method
disp('Stim site interaction in ACC ')

%% Distribution plot
sub_ind     = [sub_ind1; sub_ind2; sub_ind3; sub_ind4];
y           = [pcc_hc_cod(:); acc_hc_cod(:); pcc_amg_cod(:); acc_amg_cod(:)];
stim_site   = [pcc_hc_stimgroup(:); acc_hc_stimgroup(:); ...
    pcc_amg_stimgroup(:); acc_amg_stimgroup(:)]; % grup by stim ROI
measure_site = [pcc_hc_measgroup(:); acc_hc_measgroup(:); ...
    pcc_amg_measgroup(:); acc_amg_measgroup(:)]; % grup by measured ROI

% prepare groups to test
plot_groups = zeros(size(stim_site));
plot_groups(ismember(stim_site,'HC') & ismember(measure_site,'ACC'))    = 1;   
plot_groups(ismember(stim_site,'AMG') & ismember(measure_site,'ACC'))   = 2;
plot_groups(ismember(stim_site,'HC') & ismember(measure_site,'PCC'))    = 3;
plot_groups(ismember(stim_site,'AMG') & ismember(measure_site,'PCC'))   = 4;

% plot distribution 
figure('Position',[0 0 250 400]), hold on;
distributionPlot(y,'groups',plot_groups,'color',[.7 .7 1],'distWidth',1,'showMM',0);

plotSpread(y, 'distributionIdx', plot_groups, 'xValues',1:4,'spreadWidth',.5, 'distributionColors',[0.5, 0.5, 0.5]);

plot([1:4],[median(y(plot_groups==1)) median(y(plot_groups==2)) median(y(plot_groups==3)) median(y(plot_groups==4))],'k.','MarkerSize',30);

xlim([0 5]),ylim([-0.6 1.4]),ylabel('Explained Variance (R^2)')
set(gca,'YTick',[0:.2:1],'XTick',[1:4],'XTickLabel',{'H-PCC','A-PCC','H-ACC','A-ACC'})


