
clearvars, close all, clc
% startup

%% Get stats across subs
% Dependencies: Van den Boom matmef repository: https://github.com/MaxvandenBoom/matmef

addpath(genpath(pwd));

% set local path to your BIDS directory:
myPath = setLocalDataPath(1);
localDataPath = myPath.input;

% load the meta data
all_subjects = {'01','02','03','04','05','06','07','08'};
all_hemi = {'r','r','r','l','r','l','l','r'};
all_runs = {'01','01','01','01','01','01','01','01'};

for ss = 1:length(all_subjects)
    bids_sub = all_subjects{ss};
    bids_ses = 'ieeg01';
    bids_task = 'ccep';
    bids_run = all_runs{ss};

    % set var names: mefd, events, channels and electrodes
    fileName = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
        ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_ieeg.mefd']);
    events_tsv_name = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
        ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_events.tsv']);
    channels_tsv_name = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
        ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_channels.tsv']);
    electrodes_tsv_name = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
        ['sub-' bids_sub '_ses-' bids_ses '_electrodes.tsv']);

    % load events, channels and electrodes
    [metadata] = readMef3(fileName); % metadata table
    events_table = readtable(events_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'}); % events table
    channels_table = readtable(channels_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'}); % channels table
    electrodes_table = readtable(electrodes_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'}); % electrodes table

    % --------------------------------------------------------------------
    % get necessary parameters from data
    srate = metadata.time_series_metadata.section_2.sampling_frequency;  % sampling frequency
    nr_channels = length(metadata.time_series_channels); % number of channels
    % list of channel names
    channel_names = cell(nr_channels,1);
    for kk = 1:nr_channels
        channel_names{kk} = metadata.time_series_channels(kk).name;
    end
    % find good sEEG/ECoG channels
    good_channels = find(ismember(channels_table.type,{'ECOG','SEEG'}) & ismember(channels_table.status,'good'));
    

    % load event data, clip other stim amplitudes
    events_table_clipped = bids_clipEvents(events_table,'electrical_stimulation_current', {'4.0 mA', '6.0 mA'}); % keep only events with stim current == 4.0 or 6.0 mA

    % see which stim pairs are in areas of interest
    areas_interest = [17 18 11106 11107 11108 11109 11110 11123 10 53 54 12106 12107 12108 12109 12110 12123 49]; % Destrieux_label
    % get areas fo each channel name and each ccep stim pair name
    channel_names = channels_table.name;
    channel_areas = zeros(size(channel_names));
    ccep_stim_names = unique(events_table_clipped.electrical_stimulation_site);
    ccep_stim_areas = zeros(length(ccep_stim_names),2);

    % loop through channel names to get channel areas
    for kk = 1:length(channel_names)
        % which number has this channel in the electrodes.tsv file?
        thisElPos = find(ismember(electrodes_table.name,channel_names{kk}));
        if ismember(electrodes_table.Destrieux_label(thisElPos),areas_interest)
            channel_areas(kk) = areas_interest(ismember(areas_interest,electrodes_table.Destrieux_label(thisElPos)));
        end
    end

    % get areas for each stim pair name 
    for kk = 1:length(ccep_stim_names)
        if sum(ismember(ccep_stim_names{kk},'-'))==1 % 1 -
            ccep_stim_areas(kk,1) = channel_areas(ismember(channel_names,extractBefore(ccep_stim_names{kk},'-')));
            ccep_stim_areas(kk,2) = channel_areas(ismember(channel_names,extractAfter(ccep_stim_names{kk},'-')));
        else % assume the second - if there is a - in the channel name
            dash_in_name = find(ismember(ccep_stim_names{kk},'-'));
            el1_name = ccep_stim_names{kk}(1:dash_in_name(2)-1);
            el2_name = ccep_stim_names{kk}(dash_in_name(2)+1:end);
            ccep_stim_areas(kk,1) = channel_areas(ismember(channel_names,el1_name));
            ccep_stim_areas(kk,2) = channel_areas(ismember(channel_names,el2_name));
        end
    end

    % only keep events where we stimulated the limbic network
    ccep_stim_areas_limbic = ccep_stim_areas(sum(ccep_stim_areas,2)>0,:);
    ccep_stim_names_limbic = ccep_stim_names(sum(ccep_stim_areas,2)>0,:);
    events_table_clipped = bids_clipEvents(events_table_clipped,'electrical_stimulation_site', ccep_stim_names_limbic); % keep only events with stim current == 4.0 or 6.0 mA

    % remove events if less than 4 of the same type
    remove_stim_pair = {};
    [uni_stim_sets,bb,cc] = unique(events_table_clipped.electrical_stimulation_site);
    for kk = 1:length(uni_stim_sets)
        if length(find(ismember(events_table_clipped.electrical_stimulation_site,uni_stim_sets{kk})))<4
            disp(['less then 4 events for sub ' int2str(ss) ' ' uni_stim_sets{kk}])
            remove_stim_pair =[remove_stim_pair uni_stim_sets(kk)];
        end
    end
    events_table_clipped(ismember(events_table_clipped.electrical_stimulation_site,remove_stim_pair),:) = [];

    % stuff
    baseline_t = [-0.5 -0.05]; 
    t_win_crp = [0.015 1];
    % only calculate stats for channel_areas>0
    [average_ccep,average_ccep_names,tt,srate,crp_out] = ...
        ccepPCC_loadAverageSubset(fileName,events_table_clipped,good_channels, channel_areas, baseline_t,t_win_crp,1);
    
    % now get stim areas, incorrect order:
    [~,x_ind] = ismember(average_ccep_names,ccep_stim_names_limbic); % indiced into stim_names from ccep_names
    average_ccep_areas = ccep_stim_areas_limbic(x_ind,:);
    clear ccep_stim_areas_limbic ccep_stim_names_limbic % these were only used to clip the events table, not for computations

    outputName = fullfile(localDataPath,'derivatives','stats',['sub-' bids_sub],...
        ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_crp.mat']);
    save(outputName,'average_ccep','average_ccep_names','average_ccep_areas','tt','srate','crp_out','channel_names','channel_areas');
    clear average_ccep average_ccep_names tt srate crp_out
    disp(['sub ' bids_sub ' stats saved'])
end



