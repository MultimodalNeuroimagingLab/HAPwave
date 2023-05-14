
function [events_table,channels_table,electrodes_table,all_out] = pcc_loadAveragesStats(localDataPath,bids_sub,bids_ses,bids_task,bids_run)

% just making sure we always load stats in the same way across scripts.

% set filenames: mefd, events, channels and electrodes
events_tsv_name = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
    ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_events.tsv']);
channels_tsv_name = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
    ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_channels.tsv']);
electrodes_tsv_name = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
    ['sub-' bids_sub '_ses-' bids_ses '_electrodes.tsv']);

% load events, channels and electrodes
events_table = readtable(events_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'}); % events table
channels_table = readtable(channels_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'}); % channels table
electrodes_table = sortElectrodes(electrodes_tsv_name, channels_tsv_name, false); % electrodes table sorted according to channels

% get distances from electrodes to pial and gray/white
gL_pial = gifti(fullfile(localDataPath,'derivatives','freesurfer',['sub-' bids_sub],'pial.L.surf.gii'));
gR_pial = gifti(fullfile(localDataPath,'derivatives','freesurfer',['sub-' bids_sub],'pial.R.surf.gii'));
gL_white = gifti(fullfile(localDataPath,'derivatives','freesurfer',['sub-' bids_sub],'white.L.surf.gii'));
gR_white = gifti(fullfile(localDataPath,'derivatives','freesurfer',['sub-' bids_sub],'white.R.surf.gii'));
distLim = 6;
[dist_gw_info] = ieeg_eldist2pial2white(electrodes_table,gL_pial,gR_pial,gL_white,gR_white, distLim);


% find good sEEG channels
good_channels = find(ismember(channels_table.type,{'ECOG','SEEG'}) & ismember(channels_table.status,'good'));

% load stats/average CCEPs from outputName
crpFile = fullfile(localDataPath,'derivatives','stats',['sub-' bids_sub],...
    ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_crp.mat']);
load(crpFile,'average_ccep','average_ccep_names','average_ccep_areas','tt','srate','crp_out','channel_names','channel_areas');

all_out.average_ccep = average_ccep;
all_out.average_ccep_names = average_ccep_names;
all_out.average_ccep_areas = average_ccep_areas;

all_out.tt = tt;
all_out.srate = srate;
all_out.crp_out = crp_out;

% label areas of stim and measured electrodes
all_out.nr_channels = height(channels_table); % get number of channels
all_out.channel_names = channels_table.name; % list of channel names
all_out.channel_areas = channel_areas;

% find good sEEG channels
all_out.good_channels = find(ismember(channels_table.type,{'ECOG','SEEG'}) & ismember(channels_table.status,'good'));
all_out.bad_channels = find(~ismember(channels_table.type,{'ECOG','SEEG'}) | ~ismember(channels_table.status,'good'));

all_out.elec_relDist = dist_gw_info.rel_dist; % distance from electrode to gray/white, negative if in white matter
all_out.cortex_thick = dist_gw_info.dist_pialwhite; % distance from pial to gray/white


clear stim_area channel_area average_ccep average_ccep_names tt srate 