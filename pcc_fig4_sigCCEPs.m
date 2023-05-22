clearvars, close all, clc
% startup

%% Load stats across subs
% cd to HAP_wave repository
addpath(genpath(pwd))
% set local path to your BIDs directory:

% load the meta data
all_subjects = {'01','02','03','04','05','06','07'}; % 
all_hemi = {'r','r','r','l','r','l','l'};
all_runs = {'01','01','01','01','01','01','01',};

t_win_cod = [0.015 .5];

for ss = 1:length(all_subjects)
    bids_sub = all_subjects{ss};
    bids_ses = 'ieeg01';
    bids_task = 'ccep';
    bids_run = all_runs{ss};

    % set filenames: mefd, events, channels and electrodes
    fileName = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
        ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_ieeg.mefd']);
    events_tsv_name = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
        ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_events.tsv']);
    channels_tsv_name = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
        ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_channels.tsv']);
    electrodes_tsv_name = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
        ['sub-' bids_sub '_ses-' bids_ses '_electrodes.tsv']);

    % load events, channels and electrodes
    events_table = readtable(events_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'}); % events table
    channels_table = readtable(channels_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'}); % channels table
    electrodes_table = readtable(electrodes_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'}); % electrodes table
    % find good sEEG channels
    good_channels = find(ismember(channels_table.type,{'ECOG','SEEG'}) & ismember(channels_table.status,'good'));

    % load stats/average CCEPs from outputName
    outputName = fullfile(localDataPath,'derivatives','stats',['sub-' bids_sub],...
        ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_stats.mat']);
    load(outputName,'average_ccep','average_ccep_names','tt','srate','cross_ccep_t','cross_ccep_p');

    % FDR correction for significance
    cross_ccep_sig = NaN(size(cross_ccep_p));
    [temp_ccep_sig,p_th] = ccepPCC_fdr(cross_ccep_p(good_channels,:),0.05);
    cross_ccep_sig(good_channels,:) = temp_ccep_sig;
    
    all_out(ss).average_ccep = average_ccep;
    all_out(ss).average_ccep_names = average_ccep_names;
    all_out(ss).tt = tt;
    all_out(ss).srate = srate;
    all_out(ss).cross_ccep_t = cross_ccep_t;
    all_out(ss).cross_ccep_p = cross_ccep_p;
    all_out(ss).sig_ccep = cross_ccep_sig;
    
    % label areas of stim and measured electrodes
    all_out(ss).nr_channels = height(channels_table); % get number of channels
    all_out(ss).channel_names = channels_table.name; % list of channel names
    % find good sEEG channels
    all_out(ss).good_channels = find(ismember(channels_table.type,{'ECOG','SEEG'}) & ismember(channels_table.status,'good'));
    all_out(ss).bad_channels = find(~ismember(channels_table.type,{'ECOG','SEEG'}) | ~ismember(channels_table.status,'good'));

    % which hemisphere are we interested
    % hemi = 'r'
    hemi = all_hemi{ss};
    % loop through channel names to identify channel areas
    if hemi == 'l'
        areas_interest_L = [17 18 11106 11107 11108 11109 11110 11123 10];
        area_names = {'hippocampus','amygdala','ant cing','mid ant cing','mid post cing','post dorsal cing','post ventral cing','parahippocampal gyrus','thalamus'};
        channel_areas = zeros(size(all_out(ss).channel_names));
        stim_areas = zeros(length(all_out(ss).average_ccep_names),2);
        % loop through channel names to get channel areas
        for kk = 1:length(all_out(ss).channel_names)
            % which number has this channel in the electrodes.tsv file?
            thisElPos = find(ismember(electrodes_table.name,all_out(ss).channel_names{kk}));
            if ismember(electrodes_table.Destrieux_label(thisElPos),areas_interest_L)
                channel_areas(kk) = find(ismember(areas_interest_L,electrodes_table.Destrieux_label(thisElPos)));
            end
        end
    else
        areas_interest_R = [53 54 12106 12107 12108 12109 12110 12123 49];
        area_names = {'hippocampus','amygdala','ant cing','mid ant cing','mid post cing','post dorsal cing','post ventral cing','parahippocampal gyrus','thalamus'};
        channel_areas = zeros(size(all_out(ss).channel_names));
        stim_areas = zeros(length(all_out(ss).average_ccep_names),2);
        % loop through channel names to get channel areas
        for kk = 1:length(all_out(ss).channel_names)
            % which number has this channel in the electrodes.tsv file?
            thisElPos = find(ismember(electrodes_table.name,all_out(ss).channel_names{kk}));
            if ismember(electrodes_table.Destrieux_label(thisElPos),areas_interest_R)
                channel_areas(kk) = find(ismember(areas_interest_R,electrodes_table.Destrieux_label(thisElPos)));
            end
        end
    end

    % loop through stim pairs
    for kk = 1:length(all_out(ss).average_ccep_names)
        if sum(ismember(all_out(ss).average_ccep_names{kk},'-'))==1 % 1 -
            stim_areas(kk,1) = channel_areas(ismember(all_out(ss).channel_names,extractBefore(all_out(ss).average_ccep_names{kk},'-')));
            stim_areas(kk,2) = channel_areas(ismember(all_out(ss).channel_names,extractAfter(all_out(ss).average_ccep_names{kk},'-')));
        else % assume the second - if there is a - in the channel name
            dash_in_name = find(ismember(all_out(ss).average_ccep_names{kk},'-'));
            el1_name = all_out(ss).average_ccep_names{kk}(1:dash_in_name(2)-1);
            el2_name = all_out(ss).average_ccep_names{kk}(dash_in_name(2)+1:end);
            stim_areas(kk,1) = channel_areas(ismember(channel_names,el1_name));
            stim_areas(kk,2) = channel_areas(ismember(channel_names,el2_name));
        end
    end
    stim_area = max(stim_areas,[],2); % take the maximum area between the two stim pairs, ignoring zeros, preferring second area
    all_out(ss).stim_area = stim_area;
    all_out(ss).channel_areas = channel_areas;
    %%% --> now we have channel_areas and stim_areas
    clear stim_area channel_area average_ccep average_ccep_names tt srate cod_ccep cod_ccep_p cod_ccep_quart
end

%% Plot significant CCEPs across subs in different colors
sub_color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
sub_axis = {[-250 650],[-100 300],[-350 300],[-300 2100],[-200 350],[-500 1000],[-300 300]};% blue, orange, green, purple, mustard, celeste, wine

figure('Position',[0 0 400 800]), hold on

stim_names = {'Hipp','Amyg','Thal'};
stim_area = 1; % 1 for hippocampus, 2 for amygdala, 3 for thalamus
rec_names = {'PCC','ACC'};
rec_area = 1; % 1 for PCC, 2 for ACC

out_plot_responses_norm = [];
out_subj_ind = [];
resp_counter = 0;

for ss = 1:7 
    subplot(7,1,(ss)), hold on
    % Get sites that belong to the measurement ROI (rec_area)
    if rec_area == 1
        these_measured_sites = find(all_out(ss).channel_areas==5 | all_out(ss).channel_areas==6 | all_out(ss).channel_areas==7); % posterior cing
    elseif rec_area==2
        these_measured_sites = find(all_out(ss).channel_areas==3 | all_out(ss).channel_areas==4); % anterior cing
    end
    
    % Get sites that belong to the stimulated ROI (stim_area)
    if stim_area == 1
        these_stim_sites = find(all_out(ss).stim_area == 1 | all_out(ss).stim_area == 8); % hippocampal formation 
    elseif stim_area==2
        these_stim_sites = find(all_out(ss).stim_area == 2); % amygdala
    elseif stim_area==3
        these_stim_sites = find(all_out(ss).stim_area == 9); % thalamus
    end
    
    sign_resp = all_out(ss).sig_ccep==1; %  p<0.05 FDR corrected
    
    % Look for coverage in the stimulated and/or measurement ROI
    if isempty(these_measured_sites)
        % no coverage sites of interest 
        % draw brown dotted line             / Meaning: No electrodes implanted in ROI for this sub
        yline(ss,':', 'Color', [.9 .3 .3],'LineWidth', .5);
        
    elseif isempty(these_stim_sites) 
        % no stimulated ROI coverage        / Meaning: No SPES of interest was delivered ):
        % draw gray dotted line if no stim_sites
        yline(ss,':', 'Color', [.3 .3 .3], 'LineWidth', .5); 

    else
        % Coverage in both measurement and stimulated ROI
        % black line if no sign responses     / Meaning: CCEPs were not significant under SEPS of interest ):
        plot([-.2 .8], [ss, ss], 'Color', 'k', 'LineWidth', .7); 

        % loop over measured sites
        for kk = 1:length(these_measured_sites)
            % loop over the stimulated pairs
            for ll = 1:length(these_stim_sites)
                % test significance
                if sign_resp(these_measured_sites(kk), these_stim_sites(ll)) == 1
%                    subj_resp_sign(ss) = subj_resp_sign(ss) + 1;
                   
                   plot_responses = squeeze(all_out(ss).average_ccep(these_measured_sites(kk), these_stim_sites(ll), :));
                   plot_responses(all_out(ss).tt > -0.010 & all_out(ss).tt < 0.010) = NaN;
                   plot(all_out(ss).tt, ss + plot_responses, 'Color',sub_color{ss}, 'LineWidth', .5);
                   xlim([-0.2 .6]), ylim([sub_axis{ss}]);
                   xlabel('Time (s)'), ylabel('Voltage (uV)')
                   
                    % save outputs
                    resp_counter = resp_counter + 1;
                    out_sub_responses(resp_counter, :) = plot_responses;%_norm;
                    out_plot_responses(resp_counter, :) = plot_responses;%_norm;
                    out_subj_ind(resp_counter, :) = ss;

                end % done testing significance
            end % done looping through stim pairs
        end  % done looping through stim pairs
    end
end
