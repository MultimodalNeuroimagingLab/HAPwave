clearvars, clc, close all
addpath(genpath(pwd))

% set local path to your BIDS directory:
myPath = setLocalDataPath(1);
localDataPath = myPath.input;

% load the meta data
all_subjects = {'01','02','03','04','05','06','07','08'}; % 
all_hemi = {'r','r','r','l','r','l','l','r'};
all_runs = {'01','01','01','01','01','01','01','01'};


out_chaninfo = zeros(length(all_subjects),4);
% load the meta data
for ss = 1:length(all_subjects)
    bids_sub = all_subjects{ss};
    bids_ses = 'ieeg01';
    bids_task = 'ccep';
    bids_run = all_runs{ss};

    [events_table,channels_table,electrodes_table,sub_out] = pcc_loadAveragesStats(localDataPath,bids_sub,bids_ses,bids_task,bids_run);
    all_out(ss) = sub_out;
    channels_out(ss).chan = channels_table;
    
    
    nr_good_chan = length(find(ismember(channels_table.type,'SEEG') & ismember(channels_table.status,'good')));
    nr_bad_chan = length(find(ismember(channels_table.type,'SEEG') & ismember(channels_table.status,'bad')));
    prop_bad_chan = nr_bad_chan/(nr_good_chan + nr_bad_chan);

    sprintf(['subject ' int2str(ss) '\n'...
        'good channels: '  int2str(nr_good_chan) '\n'...
        'bad channels: '  int2str(nr_bad_chan) '\n'...
        'bad channel proportion: ' num2str(prop_bad_chan)])

    out_chaninfo(ss,1) = ss;
    out_chaninfo(ss,2) = nr_good_chan;
    out_chaninfo(ss,3) = nr_bad_chan;
    out_chaninfo(ss,4) = prop_bad_chan;

end



