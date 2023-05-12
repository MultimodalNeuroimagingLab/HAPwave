clearvars, clc, close all
% startup 

%% Load stats across subs
% Dependencies: AGOV pcc_project repository: http://birrepos.mayo.edu/multimodal-neuroimaging-lab/pcc_project
%               MNL ieeg basics repository
% cd to HAPwave repository
addpath(genpath(pwd))
% addpath(genpath('/Users/M219978/Documents/git/mnl_ieegBasics'))
% set local path to your BIDs directory:

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
   
    % find good sEEG channels
    good_channels = find(ismember(channels_table.type,{'ECOG','SEEG'}) & ismember(channels_table.status,'good'));

    % load stats/average CCEPs from outputName
    crpFile = fullfile(localDataPath,'derivatives','stats',['sub-' bids_sub],...
        ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_crp.mat']);
    load(crpFile,'average_ccep','average_ccep_names','average_ccep_areas','tt','srate','crp_out','channel_names','channel_areas');
    
    all_out(ss).average_ccep = average_ccep;
    all_out(ss).average_ccep_names = average_ccep_names;
    all_out(ss).average_ccep_areas = average_ccep_areas;

    all_out(ss).tt = tt;
    all_out(ss).srate = srate;
    all_out(ss).crp_out = crp_out;
    
    % label areas of stim and measured electrodes
    all_out(ss).nr_channels = height(channels_table); % get number of channels
    all_out(ss).channel_names = channels_table.name; % list of channel names
    all_out(ss).channel_areas = channel_areas;

    % find good sEEG channels
    all_out(ss).good_channels = find(ismember(channels_table.type,{'ECOG','SEEG'}) & ismember(channels_table.status,'good'));
    all_out(ss).bad_channels = find(~ismember(channels_table.type,{'ECOG','SEEG'}) | ~ismember(channels_table.status,'good'));

    clear stim_area channel_area average_ccep average_ccep_names tt srate 
end



%% Load stim and measure areas across subjects
sub_color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
% blue, orange, green, purple, mustard, celeste, wine

area_names = {'Hipp','Amyg','PCC','ACC'};   
area_codes_r = {[12123 53],[54],[12108 12109 12110],[12106 12107]}; % right
area_codes_l = {[11123 17],[18],[11108 11109 11110],[11106 11107]}; % left

sub_hemi = {'r','r','r','l','r','l','l','r'};

%
crp_var = 1; %1 for alpha prime, 4 for COD


for ss = 8%:8
   

    if isequal(sub_hemi{ss},'l')
        area_codes = area_codes_l;
    elseif isequal(sub_hemi{ss},'r')
        area_codes = area_codes_r;
    end

    figure
    
    pcc_els = find(ismember(all_out(ss).channel_areas,area_codes{3}));
    pcc_els = setdiff(pcc_els,all_out(ss).bad_channels); % only good channels
    
    acc_els = find(ismember(all_out(ss).channel_areas,area_codes{4}));
    acc_els = setdiff(acc_els,all_out(ss).bad_channels); % only good channels
    
    resp_ampl_pcc = NaN(length(all_out(ss).average_ccep_names), length(pcc_els), 4);
    resp_ampl_acc = NaN(length(all_out(ss).average_ccep_names), length(acc_els), 4);
    
    for kk = 1: length(all_out(ss).average_ccep_names) 
    
        if ~isempty(all_out(ss).crp_out(pcc_els(1),kk).crp_parms) % there are sufficient trials
    
            if ismember(all_out(ss).average_ccep_areas(kk,1),area_codes{1}) || ...
                ismember(all_out(ss).average_ccep_areas(kk,2),area_codes{1})% hippocampal stim
        
                % get response info from ACC and PCC
                for ii = 1:length(pcc_els)
                    resp_ampl_pcc(kk,ii,1) = mean(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.al_p); % mean across trials
                    resp_ampl_pcc(kk,ii,2) = mean(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.Vsnr); % mean across trials | S1 mean 0.0765
                    resp_ampl_pcc(kk,ii,3) = mean(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.expl_var); % mean across trials
                    resp_ampl_pcc(kk,ii,4) = mean(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.cod); % mean across trials
                    
                    subplot(2,3,1), hold on
    %                 plot(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.parms_times, all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.C,'b')
    %                 plot(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.parms_times, mean(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.V_tR,2),'b')
                    plot(all_out(ss).tt, squeeze(all_out(ss).average_ccep(pcc_els(ii),kk,:)),'b')
    
                    xlim([-.5 1]), title(['S' int2str(ss) ' HP -> PCC'])
                    ylim([-200 600])
                end

                for ii = 1:length(acc_els)
                    resp_ampl_acc(kk,ii,1) = mean(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.al_p); % mean across trials
                    resp_ampl_acc(kk,ii,2) = mean(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.Vsnr); % mean across trials
                    resp_ampl_acc(kk,ii,3) = mean(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.expl_var); % mean across trials
                    resp_ampl_acc(kk,ii,4) = mean(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.cod); % mean across trials
        
%                     if mean(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.cod) > 1
%                         ss
%                         all_out(ss).channel_names{acc_els(ii)}
%                         all_out(ss).average_ccep_names{kk}
%                     end


                    subplot(2,3,2), hold on
    %                 plot(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.parms_times, all_out(ss).crp_out(acc_els(ii),kk).crp_parms.C,'r')
    %                 plot(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.parms_times, mean(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.V_tR,2),'r')
                    plot(all_out(ss).tt, squeeze(all_out(ss).average_ccep(acc_els(ii),kk,:)),'r')
                    
                    xlim([-.5 1]), title('HP -> ACC')
                    ylim([-200 600])
                end

        
            elseif ismember(all_out(ss).average_ccep_areas(kk,1),area_codes{2}) || ...
                ismember(all_out(ss).average_ccep_areas(kk,2),area_codes{2})% amygdala stim
        
                % get response info from ACC and PCC
                for ii = 1:length(pcc_els)
                    resp_ampl_pcc(kk,ii,1) = mean(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.al_p); % mean across trials
                    resp_ampl_pcc(kk,ii,2) = mean(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.Vsnr); % mean across trials
                    resp_ampl_pcc(kk,ii,3) = mean(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.expl_var); % mean across trials
                    resp_ampl_pcc(kk,ii,4) = mean(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.cod); % mean across trials
        
                    subplot(2,3,4), hold on
    %                 plot(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.parms_times, all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.C,'b')
    %                 plot(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.parms_times, mean(all_out(ss).crp_out(pcc_els(ii),kk).crp_parms.V_tR,2),'b')
                    plot(all_out(ss).tt, squeeze(all_out(ss).average_ccep(pcc_els(ii),kk,:)),'b')
                    
                    xlim([-.5 1]), title('AMYG -> PCC')
                    ylim([-200 600])
        
                end

                for ii = 1:length(acc_els)
                    resp_ampl_acc(kk,ii,1) = mean(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.al_p); % mean across trials
                    resp_ampl_acc(kk,ii,2) = mean(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.Vsnr); % mean across trials
                    resp_ampl_acc(kk,ii,3) = mean(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.expl_var); % mean across trials
                    resp_ampl_acc(kk,ii,4) = mean(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.cod); % mean across trials

                    subplot(2,3,5), hold on
    %                 plot(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.parms_times, all_out(ss).crp_out(acc_els(ii),kk).crp_parms.C,'r')
    %                 plot(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.parms_times, mean(all_out(ss).crp_out(acc_els(ii),kk).crp_parms.V_tR,2),'r')
                    plot(all_out(ss).tt, squeeze(all_out(ss).average_ccep(acc_els(ii),kk,:)),'r')
                    
                    xlim([-.5 1]), title('AMYG -> ACC')
                    ylim([-200 600])
        
                end
        
            end
        end
    end
    
    % find all hipp sites and compare response amplitude in PCC and ACC
    hc_stims = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{1}) | ...
            ismember(all_out(ss).average_ccep_areas(:,2),area_codes{1})); % hippocampal stim
    amg_stims = find(ismember(all_out(ss).average_ccep_areas(:,1),area_codes{2}) | ...
            ismember(all_out(ss).average_ccep_areas(:,2),area_codes{2})); % amygdala stim
           
    disp('CRP parameters');
    subplot(2,3,3), hold on % hc stim
    plot(1, nanmean( resp_ampl_pcc(hc_stims,:,crp_var), 2), '*'),  % average pcc response to hipp 
    plot(2, nanmean( resp_ampl_acc(hc_stims,:,crp_var), 2), '*') % average acc response to hipp
    xlim([0 3])
    % display crp_var average
    title(['HC stim'])
    pcc_mean = nanmean(nanmean( resp_ampl_pcc(hc_stims,:,crp_var), 2))
    text(0.5,min(nanmean( resp_ampl_pcc(hc_stims,:,crp_var), 2)),{pcc_mean},'Rotation',90);%,'VerticalAlignment','baseline');
    acc_mean = nanmean(nanmean( resp_ampl_acc(hc_stims,:,crp_var), 2))
    text(1.5,min(nanmean( resp_ampl_acc(hc_stims,:,crp_var), 2)),{acc_mean},'Rotation',90);%,'VerticalAlignment','baseline');


    subplot(2,3,6), hold on % amg stim
    plot(1, nanmean( resp_ampl_pcc(amg_stims,:,crp_var), 2), '*') % average pcc response to amyg
    plot(2, nanmean( resp_ampl_acc(amg_stims,:,crp_var), 2), '*') % average acc response to amyg
    xlim([0 3])
    % add & display crp_var average
    title(['Amg stim'])
    pcc_mean = nanmean(nanmean( resp_ampl_pcc(amg_stims,:,crp_var), 2))
    text(0.5,min(nanmean( resp_ampl_pcc(amg_stims,:,crp_var), 2)),{pcc_mean},'Rotation',90);%,'VerticalAlignment','baseline');
    acc_mean = nanmean(nanmean( resp_ampl_acc(amg_stims,:,crp_var), 2))
    text(1.5,min(nanmean( resp_ampl_acc(amg_stims,:,crp_var), 2)),{acc_mean},'Rotation',90);%,'VerticalAlignment','baseline');

end

%% now show distribution plots
addpath(genpath('/Users/M219978/Documents/git/measures-of-effect-size-toolbox'))

% list all CRP amplitudes in one column
% pick a variable from the output
crp_var = 1; %1 for alpha prime, 4 for COD

% 1,1 list all PCC measured HC stim
pcc_hc_amp = resp_ampl_pcc(hc_stims,:,crp_var);
% generate two groups for this: one for stim, one for measured
pcc_hc_stimgroup = ones(size(pcc_hc_amp)); % hc stim gets code 1
pcc_hc_measgroup = 3*ones(size(pcc_hc_amp)); % pcc measure gets code 3

% 1,2 list all ACC measured HC stim
acc_hc_amp = resp_ampl_acc(hc_stims,:,crp_var);
% generate two groups for this: one for stim, one for measured
acc_hc_stimgroup = ones(size(acc_hc_amp)); % hc stim gets code 1
acc_hc_measgroup = 4*ones(size(acc_hc_amp)); % acc measure gets code 4

% 2,1 list all PCC measured AMG stim
pcc_amg_amp = resp_ampl_pcc(amg_stims,:,crp_var);
% generate two groups for this: one for stim, one for measured
pcc_amg_stimgroup = 2*ones(size(pcc_amg_amp)); % amg stim gets code 2
pcc_amg_measgroup = 3*ones(size(pcc_amg_amp)); % pcc measure gets code 3

% 2,2 list all ACC measured AMG stim
acc_amg_amp = resp_ampl_acc(amg_stims,:,crp_var);
% generate two groups for this: one for stim, one for measured
acc_amg_stimgroup = 2*ones(size(acc_amg_amp)); % amg stim gets code 2
acc_amg_measgroup = 4*ones(size(acc_amg_amp)); % acc measure gets code 4

X = [pcc_hc_amp(:); acc_hc_amp(:); pcc_amg_amp(:); acc_amg_amp(:)];

group1 = [pcc_hc_stimgroup(:); acc_hc_stimgroup(:); pcc_amg_stimgroup(:); acc_amg_stimgroup(:)]; % stim group
group2 = [pcc_hc_measgroup(:); acc_hc_measgroup(:); pcc_amg_measgroup(:); acc_amg_measgroup(:)]; % measurement group

% to plot
% mesdplot(X,groupIx,nSample,factor,isDep,fName,contrast)
[stats,sumTable] = mes2way(X,[group1 group2],{'omega2','eta2'},'doDataPlot',true);%,'cWeight',[1 -1;-1 1]);

plot_groups = zeros(size(group1));
plot_groups(group1==1 & group2==3) = 1;
plot_groups(group1==2 & group2==3) = 2;
plot_groups(group1==1 & group2==4) = 3;
plot_groups(group1==2 & group2==4) = 4;

X = [X; NaN; NaN; NaN; NaN];
plot_groups = [plot_groups; 1; 2; 3; 4];

% figure,hold on
figure('Position',[0 0 250 400]), hold on;
distributionPlot(X,'groups',plot_groups,'addSpread',1,'showMM',5);%
str = sumTable{5,6};
text(2.5,max(X),{str},'HorizontalAlignment','center','VerticalAlignment','cap');%,'Rotation',90,'VerticalAlignment','baseline');

if sumTable{5,6} < 0.05 % if p <-0.05
    plot(2.5,max(X),'r*', 'MarkerSize',20)
end
xlim([0 5])
