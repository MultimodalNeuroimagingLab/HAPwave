function [average_ccep,average_ccep_names,tt,srate,crp_out,single_trials] = ...
    ccepPCC_loadAverageSubset(fileName,events_table,good_channels, channel_areas ,baseline_t,t_win_cod,use_CAR,varargin)
%
% input:
%   fileName
%   events_tsv_name
%   good_channels
%   baseline_t: in seconds 
%   t_win_cod: for stats
%   use_CAR: default  = 1
%   channel_areas: only calculate stats for channel areas >0
%
% output:
%   average_ccep
%   average_ccep_names
%   tt
%   srate
%   cross_ccep_t: 
%   cross_ccep_p:
% 
% DH, HH and GO, Multimodal Neuroimaging Lab, Mayo Clinic 2021

if isempty(use_CAR)
    use_CAR = 1; % default is to do CAR
end

if isempty(varargin)
    return_single_trials = 0;
else 
    if length(varargin)>=1
        if ~isempty(varargin{1})
            return_single_trials = varargin{1};
        else
            return_single_trials = 0;
        end
    end
end

% load metadata
[metadata]      = readMef3(fileName);

%% get necessary parameters from data

% sampling frequency
srate           = metadata.time_series_metadata.section_2.sampling_frequency;

% number of channels
nr_channels     = length(metadata.time_series_channels);

% list of channel names
channel_names   = cell(nr_channels,1);
for kk = 1:nr_channels
    channel_names{kk} = metadata.time_series_channels(kk).name;
end

%% make a stim_pair vector to average epochs
%   - group F01-F02 with F02-F01
%   - stim_pair_nr & stim_pair_name have length of the events table
%   - stim_pair_nr indicates a condition number, starts at 1
%   - stim_pair_name contains the stim_pair (F02-F01 is switched to F01-F02)

% initialize stim_pair name and number (conditions vectors)
stim_pair_nr    = NaN(height(events_table),1);                   % number of pair (e.g. 1)
stim_pair_name  = cell(height(events_table),1);                  % name of pair (e.g. LTG1-LTG2) 

% get all stim + electrodes
stimEl1         = extractBefore(events_table.electrical_stimulation_site,'-');
% get all stim - electrodes
stimEl2         = extractAfter(events_table.electrical_stimulation_site,'-');

condition_type_counter = 0;
for kk = 1:height(events_table)

    % which electrodes are stimulated
    el1 = stimEl1{kk};
    el2 = stimEl2{kk};
        
    % is this trial a stimulation trial: does el1 have content?
    if ~isempty(el1) && ~isempty(el2)
        % if this trial type does not exist yet & is a stimulation trial
        if sum(strcmp(stim_pair_name,[el1 '-' el2]))==0 && ...  % does el1-el2 already exist?
                sum(strcmp(stim_pair_name,[el2 '-' el1]))==0    % group el2-el1 with el1-el2
            condition_type_counter = condition_type_counter+1;

            % find all trials with el1 & el2 | el2 & el1
            theseTrials = strcmp(stimEl1,el1) & strcmp(stimEl2,el2) | ...
                strcmp(stimEl2,el1) & strcmp(stimEl1,el2);
            trial_nrs = find(theseTrials==1);                   % number of trials of this type
            for ll = 1:sum(theseTrials)
                stim_pair_name{trial_nrs(ll),1} = [el1 '-' el2];
            end
            stim_pair_nr(theseTrials) = condition_type_counter;
        end
    end
end
clear el1 el2 trial_nrs theseTrials epoch_type_counter

%% load data for each condition and put average in matrix
% average_ccep: matrix with average cceps (measured electrodes X stim pairs X time)
% average_ccep_names: names of stimulated channels

% do not include bad epochs: interictal activity & wrong onset
epochs_include          = ones(size(stim_pair_nr));
epochs_include(~ismember(events_table.status,'good')) = 0;

% set epoch parameters
epoch_length            = 5; % in seconds, -1:3
epoch_prestim_length    = 2; % in seconds

tt = (1:epoch_length*srate)/srate - epoch_prestim_length;

% initialize output
average_ccep            = NaN(nr_channels,max(stim_pair_nr),epoch_length*srate);
average_ccep_names      = cell(max(stim_pair_nr),1);
crp_out                 = [];

for kk = 1:max(stim_pair_nr) % condition number
    disp(['loading data for condition ' int2str(kk) ' out of ' int2str(max(stim_pair_nr))])
    
    % epochs of this condition
    these_epochs        = find(stim_pair_nr==kk & epochs_include==1);
    
    if ~isempty(these_epochs)
    
        % save name of the current epoch
        average_ccep_names{kk} = stim_pair_name{these_epochs(1)};

        % for this condition number (kk), load each epoch    
        all_start           = round((events_table.onset(these_epochs)-epoch_prestim_length)*srate);
        all_end             = round((events_table.onset(these_epochs)+epoch_length-epoch_prestim_length)*srate);
        epoch_ranges        = [all_start all_end];
        [~,signaldata]      = readMef3(fileName,[], [], 'samples', epoch_ranges);  % read all channels, samples 0-1000 
        
        % exclude stimulated channels from the good channels
        stimEl1_nr          = find(ismember(channel_names,extractBefore(average_ccep_names{kk},'-')));
        stimEl2_nr          = find(ismember(channel_names,extractAfter(average_ccep_names{kk},'-')));

        good_channels_car   = setdiff(good_channels,[stimEl1_nr stimEl2_nr]);

        % adjusted Common Average Reference (CAR)
        if use_CAR == 1
            perc_channels   = 0.2;
            car_timeint     = [0.015 0.500];
            [signaldata]    = ccep_CAR64blocks_percent(signaldata,tt,good_channels_car,perc_channels,car_timeint);
        end

        signaldata          = permute(signaldata,[1 3 2]);

        % baseline subtract
        samples_base        = find(tt>baseline_t(1) & tt<baseline_t(2));
        data_epoch          = ieeg_baselinesubtract(signaldata,samples_base,'median');
        
        % run CRP for channels in the limbic network
        for ii = 1:size(channel_areas,1)
            if channel_areas(ii)>1  && ismember(ii,good_channels_car) % && size(data_epoch,2)>3% only limbic 
                crp_out(ii,kk).data = squeeze(data_epoch(ii,:,:));
                crp_out(ii,kk).tt   = tt;
                V                   = squeeze(data_epoch(ii,:,tt>t_win_cod(1) & tt<t_win_cod(2)));
                t_win               = tt(tt>t_win_cod(1) & tt<t_win_cod(2));
                [crp_parms, crp_projs] = CRP_method(V',t_win);
                crp_out(ii,kk).crp_parms = crp_parms;
                crp_out(ii,kk).crp_projs = crp_projs;
            else
                crp_out(ii,kk).data = [];
                crp_out(ii,kk).tt   = [];
                crp_out(ii,kk).crp_parms = [];
                crp_out(ii,kk).crp_projs = [];
            end
        end

        % put average in matrix
        average_ccep(:,kk,:)        = squeeze(nanmean(data_epoch,2));

    else
        % save name of the current epoch
        average_ccep_names{kk}      = stim_pair_name{find(stim_pair_nr==kk,1)};
    end
    if return_single_trials         == 1 && max(stim_pair_nr) == 1
        single_trials = data_epoch;
    elseif return_single_trials     == 1 && max(stim_pair_nr) > 1
        disp('can not return single trials, only works for 1 stim pair')
        single_trials = [];
    end
    clear these_epochs_data data_epoch ll_start ll_end
end

%% make stimulated electrodes a NaN

for kk = 1:size(average_ccep,2) % epochs
    % stimulated electrode names
    el1     = extractBefore(average_ccep_names{kk},'-');
    el2     = extractAfter(average_ccep_names{kk},'-');
    
    % stimulated electrode index in the data
    el1_nr  = ismember(channel_names,el1);
    el2_nr  = ismember(channel_names,el2);
    
    % set to NaN
    average_ccep(el1_nr==1,kk,:) = NaN;
    average_ccep(el2_nr==1,kk,:) = NaN;
    
    clear el1 el2 el1_nr el2_nr% housekeeping
end

