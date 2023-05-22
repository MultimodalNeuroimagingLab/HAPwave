
clearvars, close all, clc

%% Plot pannels with CCEPs coming from single-pair
addpath(genpath(pwd))

% set local path to your BIDS directory:
myPath = setLocalDataPath(1);
localDataPath = myPath.input;

% load the meta-data
all_subjects = {'01','02','03','04','05','06','07','08','09'};
all_hemi = {'r','r','r','l','r','l','l','r','r'};
all_runs = {'01','01','01','01','01','01','01','01','01'};

ss = 1 ;%:lenght(all_subjects)
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

baseline_t = [-0.5 -0.05];

% make a stim_pair vector to average epochs
%   - group F01-F02 with F02-F01
%   - stim_pair_nr & stim_pair_name have length of the events table
%   - stim_pair_nr indicates a condition number, starts at 1
%   - stim_pair_name contains the stim_pair (F02-F01 is switched to F01-F02)

% clip only HC outcoing CCEPs
events_table_clipped = bids_clipEvents(events_table,'electrical_stimulation_site', {'RC1-RC2','RC2-RC3','RC3-RC4','RC4-RC5'});

% initialize stim_pair name and number (conditions vectors)
stim_pair_nr = NaN(height(events_table_clipped),1); % number of pair (e.g. 1)
stim_pair_name = cell(height(events_table_clipped),1); % name of pair (e.g. LTG1-LTG2)

% get all stim + electrodes
stimEl1 = extractBefore(events_table_clipped.electrical_stimulation_site,'-');
% get all stim - electrodes
stimEl2 = extractAfter(events_table_clipped.electrical_stimulation_site,'-');

condition_type_counter = 0;
for kk = 1:height(events_table_clipped)

    % which electrodes are stimulated
    el1 = stimEl1{kk};
    el2 = stimEl2{kk};

    % is this trial a stimulation trial: does el1 have content?
    if ~isempty(el1) && ~isempty(el2)
        % if this trial type does not exist yet & is a stimulation trial
        if sum(strcmp(stim_pair_name,[el1 '-' el2]))==0 && ... % does el1-el2 already exist?
                sum(strcmp(stim_pair_name,[el2 '-' el1]))==0 % group el2-el1 with el1-el2
            condition_type_counter = condition_type_counter+1;

            % find all trials with el1 & el2 | el2 & el1
            theseTrials = strcmp(stimEl1,el1) & strcmp(stimEl2,el2) | ...
                strcmp(stimEl2,el1) & strcmp(stimEl1,el2);
            trial_nrs = find(theseTrials==1); % number of trials of this type
            for ll = 1:sum(theseTrials)
                stim_pair_name{trial_nrs(ll),1} = [el1 '-' el2];
            end
            stim_pair_nr(theseTrials) = condition_type_counter;
        end
    end
end
clear el1 el2 trial_nrs theseTrials epoch_type_counter


%% plot one condition/stim pair to check things

el_stim = {'RC1-RC2','RC2-RC3','RC3-RC4'};%,'RC4-RC5'};
el_record = 'RY2';
ref_chan = 'RX8'; % leave [] if not using 

% do not include bad epochs 
epochs_include = ones(size(stim_pair_nr));
epochs_include(~ismember(events_table_clipped.status,'good')) = 0;

% set epoch parameters
epoch_length = 4; % in seconds, -1:3
epoch_prestim_length = 1; % in seconds
ttt = (1:epoch_length*srate)/srate - epoch_prestim_length;

figure('Position',[0 0 550 800]), hold on;
for kk = 1:length(el_stim)

    % epochs of this condition
    these_epochs = find(ismember(stim_pair_name,el_stim(kk)) & epochs_include==1);

    % for this condition number (kk), load each epoch
    all_start = round((events_table_clipped.onset(these_epochs)-epoch_prestim_length)*srate);
    all_end = round((events_table_clipped.onset(these_epochs)+epoch_length-epoch_prestim_length)*srate);
    epoch_ranges = [all_start all_end];
    [~,signaldata] = readMef3(fileName,[], el_record, 'samples', epoch_ranges);  % read all channels, samples 0-1000
    signaldata = permute(signaldata,[1 3 2]);

    if ~isempty(ref_chan)
        [~,refdata] = readMef3(fileName,[], ref_chan, 'samples', epoch_ranges);  % read all channels, samples 0-1000
        refdata = permute(refdata,[1 3 2]);
        refdata = mean(refdata,1); % average across reference channels
        signaldata = signaldata-repmat(refdata,size(signaldata,1),1,1);
    end

    % baseline subtract
    samples_base = find(ttt>-1 & ttt<-0.1);
    data_epoch = ieeg_baselinesubtract(signaldata,samples_base,'median');
    
    % plot single and average responses
    subplot(3,1,kk), hold on
%     plot_responses = (squeeze(data_epoch));
    plot_responses = squeeze(signaldata);
    plot(ttt(ttt >0.010 & ttt <2),zeros(size(ttt(ttt >0.010 & ttt <2))),'k:') % plot zero line
    plot_responses(ttt > -0.010 & ttt < 0.010) = NaN; % do not plot stim artifact period
    plot(ttt,plot_responses,'Color',[.7 .7 .7 .3]) % single responses in gray
    plot(ttt,mean(plot_responses),'k','LineWidth',1.5) % mean in thicker black line
% %     Confidence interval doesn't work
%     plotCurvConf(ttt, plot_responses,'k');% plots 95% Confidende interval in gray with 50% transparency
    
    % set axis and labels
    xlabel('Time(s)'), xlim([-.2 2]);
    ylabel('Amplitude (uV)'),ylim([-400 1000]);
%     title((['Sub-' bids_sub ' Hippocampal ' el_stim(kk) ' in ' el_record]),'FontSize',3);
end


%% Plot responses from sub 8
% load the meta-data
all_subjects = {'01','02','03','04','05','06','07','08','09'};
all_hemi = {'r','r','r','l','r','l','l','r','r'};
all_runs = {'01','01','01','01','01','01','01','01','01'};

ss = 8;%:lenght(all_subjects)
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

baseline_t = [-0.5 -0.05];

% make a stim_pair vector to average epochs
%   - group F01-F02 with F02-F01
%   - stim_pair_nr & stim_pair_name have length of the events table
%   - stim_pair_nr indicates a condition number, starts at 1
%   - stim_pair_name contains the stim_pair (F02-F01 is switched to F01-F02)

% clip only 55 trials stim pair for sub 8
events_table_clipped = bids_clipEvents(events_table,'electrical_stimulation_site', {'RC1-RC2'});

% initialize stim_pair name and number (conditions vectors)
stim_pair_nr = NaN(height(events_table_clipped),1); % number of pair (e.g. 1)
stim_pair_name = cell(height(events_table_clipped),1); % name of pair (e.g. LTG1-LTG2)

% get all stim + electrodes
stimEl1 = extractBefore(events_table_clipped.electrical_stimulation_site,'-');
% get all stim - electrodes
stimEl2 = extractAfter(events_table_clipped.electrical_stimulation_site,'-');

condition_type_counter = 0;
for kk = 1:height(events_table_clipped)

    % which electrodes are stimulated
    el1 = stimEl1{kk};
    el2 = stimEl2{kk};

    % is this trial a stimulation trial: does el1 have content?
    if ~isempty(el1) && ~isempty(el2)
        % if this trial type does not exist yet & is a stimulation trial
        if sum(strcmp(stim_pair_name,[el1 '-' el2]))==0 && ... % does el1-el2 already exist?
                sum(strcmp(stim_pair_name,[el2 '-' el1]))==0 % group el2-el1 with el1-el2
            condition_type_counter = condition_type_counter+1;

            % find all trials with el1 & el2 | el2 & el1
            theseTrials = strcmp(stimEl1,el1) & strcmp(stimEl2,el2) | ...
                strcmp(stimEl2,el1) & strcmp(stimEl1,el2);
            trial_nrs = find(theseTrials==1); % number of trials of this type
            for ll = 1:sum(theseTrials)
                stim_pair_name{trial_nrs(ll),1} = [el1 '-' el2];
            end
            stim_pair_nr(theseTrials) = condition_type_counter;
        end
    end
end
clear el1 el2 trial_nrs theseTrials epoch_type_counter

%% plot one condition/stim pair to check things

el_stim = {'RC1-RC2'}; % stim site with 55 trials
el_record = {'RZ1','RPO1','RQ1'}; % all PCC sites
ref_chan = 'RX7'; % reference to white matter. Leave [] if not using

use_CAR = 1; % reference to adjCAR

% do not include bad epochs 
epochs_include = ones(size(stim_pair_nr));
% % epochs_include(events_table.status~=1) = 0;
epochs_include(~ismember(events_table.status,'good')) = 0;

% set epoch parameters
epoch_length = 5; % in seconds, -1:3
epoch_prestim_length = 2; % in seconds
ttt = (1:epoch_length*srate)/srate - epoch_prestim_length;

figure('Position',[0 0 550 800]), hold on;
for kk = 1:length(el_record)

    % epochs of this condition
    these_epochs = find(ismember(stim_pair_name,el_stim) & epochs_include==1);

    % for this condition number (kk), load each epoch
    all_start = round((events_table_clipped.onset(these_epochs)-epoch_prestim_length)*srate);
    all_end = round((events_table_clipped.onset(these_epochs)+epoch_length-epoch_prestim_length)*srate);
    epoch_ranges = [all_start all_end];
    [~,signaldata] = readMef3(fileName,[], el_record(kk), 'samples', epoch_ranges);  % read all channels, samples 0-1000
      
    % reference to a 'CAR' electrode
    if use_CAR == 1
%         [signaldata] = ccepPCC_CAR64blocks(signaldata,tt,good_channels); % old CAR
        perc_channels = 0.2;
        car_timeint = [0.015 0.500];
        [signaldata] = ccep_CAR64blocks_percent(signaldata,ttt,good_channels,perc_channels,car_timeint);
    end

%  signaldata = permute(signaldata,[1 3 2]);
%     if ~isempty(ref_chan)
%         [~,refdata] = readMef3(fileName,[], ref_chan, 'samples', epoch_ranges);  % read all channels, samples 0-1000
%         refdata = permute(refdata,[1 3 2]);
%         refdata = mean(refdata,1); % average across reference channels
%         signaldata = signaldata-repmat(refdata,size(signaldata,1),1,1);
%     end
% 
    % baseline subtract
    samples_base = find(ttt>-1 & ttt<-0.1);
    data_epoch = ieeg_baselinesubtract(signaldata,samples_base,'median');
    
    % plot single and average responses
    subplot(3,1,kk), hold on
    plot_responses = (squeeze(data_epoch));
    plot(ttt(ttt>.010 & ttt<2),zeros(size(ttt(ttt>.010 & ttt<2))),'k:') % plot zero line
    plot_responses(ttt>-0.010 & ttt<0.010) = NaN; % do not plot stim artifact period
    plot(ttt,plot_responses,'Color',[.7 .7 .7 .3]) % single responses in gray
    plot(ttt,mean(plot_responses),'k','LineWidth',1.5) % mean in thicker black line
% %     Confidence interval doesn't work
%     plotCurvConf(ttt, plot_responses,'k');% plots 95% Confidende interval in gray with 50% transparency

    % set axis and labels
    xlabel('Time(s)'), xlim([-.2 2]);
    ylabel('Amplitude (uV)'),ylim([-400 920]);
%     title((['Sub-' bids_sub ' Hippocampal ' el_stim(kk) ' in ' el_record]),'FontSize',3);
end
