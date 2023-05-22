
clearvars, close all, clc
startup

%% Plot figure showing a delay in CCEP depending on stim site
addpath(genpath(pwd));

% set local path to your BIDS directory:
myPath = setLocalDataPath(1);
localDataPath = myPath.input;

% Load variables and stats
all_subjects = {'01','02','03','04','05','06','07'}; % S4, S5 & perhaps S3 show temporal delay
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
    % find good sEEG/ECoG channels
    all_out(ss).good_channels = find(ismember(channels_table.type,{'ECOG','SEEG'}) & ismember(channels_table.status,'good'));
    all_out(ss).bad_channels = find(~ismember(channels_table.type,{'ECOG','SEEG'}) | ~ismember(channels_table.status,'good'));

    % which hemisphere are we interested
    % hemi = 'r'
    hemi = all_hemi{ss};
    if hemi == 'l'
        areas_interest_L = [17 18 11106 11107 11108 11109 11110 11123 10];
        area_names = {'hippocampus','amygdala','ant cing','mid ant cing','mid post cing','post dorsal cing','post ventral cing','parahippocampal gyrus', 'thalamus'};
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


%% Sub-04, plot temporal delay
% rgb_color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
% blue, orange, green, purple, mustard, celeste, wine
figure('Position',[0 0 500 350]), hold on; %('Position',[0 0 600 200]), hold on

stim_area = 1; % 1 for hippocampus, 2 for amygdala, 3 for thalamus

out_plot_responses_norm = [];
out_subj_ind = [];
resp_counter = 0;
for ss = 4 % clear temporal shift, two hipp probes
    these_measured_sites = find(all_out(ss).channel_areas==5 | all_out(ss).channel_areas==6 | all_out(ss).channel_areas==7); % posterior cing
    if stim_area==1
        these_stim_sites = find(all_out(ss).stim_area==1 | all_out(ss).stim_area==8); % hippocampal formation 
    elseif stim_area==2
        these_stim_sites = find(all_out(ss).stim_area==2); % amygdala
    elseif stim_area==3
        these_stim_sites = find(all_out(ss).stim_area==9); % thalamus
    end
  
    plot(all_out(ss).tt, zeros(size(all_out(ss).tt)),'k')
    
    sign_resp = all_out(ss).sig_ccep==1; % plot p<0.05 FDR corrected
    
    % walk through pairs, test significance and plot if so
    for kk = 1:length(these_measured_sites)
        for ll = 1:length(these_stim_sites)
            resp_peak_late = [];
            resp_peak_early = [];

            if sign_resp(these_measured_sites(kk),these_stim_sites(ll))==1
                plot_responses = squeeze(all_out(ss).average_ccep(these_measured_sites(kk),these_stim_sites(ll),:));
                plot_responses(all_out(ss).tt>-0.010 & all_out(ss).tt<0.010) = NaN;

                % scaling to L2 norm (https://en.wikipedia.org/wiki/Feature_scaling)
                % unit length taken in same window as stats
                response_vector_length = sum(plot_responses(all_out(ss).tt>t_win_cod(1) &  all_out(ss).tt<t_win_cod(2)).^2).^.5;
                plot_responses_norm = plot_responses./(response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial

                % we looked at the data and visually split out 2 stim pairs
                % that produce different latencies
                if ll<=2 % stim pair 1 and 2 from LB probe
                    plot(all_out(ss).tt, plot_responses_norm,'color',[0.6350 0.0780 0.1840])
                    disp(['Red, later' all_out(ss).average_ccep_names(these_stim_sites(ll))])
                    resp_peak_late = [];

                elseif ll>=3 % stim pair 3:5 from LC probe
                    plot(all_out(ss).tt, plot_responses_norm,'color',[0.4660 0.6740 0.1880])
                    disp(['Green, early' all_out(ss).average_ccep_names(these_stim_sites(ll))])
                    resp_peak_early = [];
                end

                % save outputs
                resp_counter = resp_counter + 1;
                out_plot_responses_norm(resp_counter,:) = plot_responses_norm;
                out_subj_ind(resp_counter,:) = ss;

                % find peaks
                mx_peaks = max(out_plot_responses_norm);
                peak_time = find(mx_peaks);
                
            end
        end
    end
end
title('Green=LC1-LC5; Red=LB3-LB5')
xlim([-0.2 .6])%, ylim([-2000 4000]);
xlabel('Time (s)'),ylabel('Normilized Voltage') 

set(gcf,'PaperPositionMode','auto')
% print(fullfile('./local',figName),'-dpng');%,'-r700';
% print(fullfile('./local',figName),'-painters','-depsc')%,'-r500',)

%% Sub-05, plot temporal delay
% rgb_color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840],[0.4940 0.1840 0.5560]};
% blue, orange, green, purple, mustard, celeste, wine, purple
figure('Position',[0 0 500 350]), hold on; %('Position',[0 0 600 200]), hold on

stim_area = 1; % 1 for hippocampus, 2 for amygdala, 3 for thalamus

out_plot_responses_norm = [];
out_subj_ind = [];
resp_counter = 0;
for ss = 5 % S5 with two parahipp probes
    these_measured_sites = find(all_out(ss).channel_areas==5 | all_out(ss).channel_areas==6 | all_out(ss).channel_areas==7); % posterior cing
    if stim_area==1
        these_stim_sites = find(all_out(ss).stim_area==1 | all_out(ss).stim_area==8); % hippocampal formation
    elseif stim_area==2
        these_stim_sites = find(all_out(ss).stim_area==2); % amygdala
    elseif stim_area==3
        these_stim_sites = find(all_out(ss).stim_area==9); % thalamus
    end

    plot(all_out(ss).tt, zeros(size(all_out(ss).tt)),'k')

    sign_resp = all_out(ss).sig_ccep==1; % plot p<0.05 FDR corrected

    % walk through pairs, test significance and plot if so
    for kk = 1:length(these_measured_sites)
        for ll = 1:length(these_stim_sites)
            if sign_resp(these_measured_sites(kk),these_stim_sites(ll))==1
                plot_responses = squeeze(all_out(ss).average_ccep(these_measured_sites(kk),these_stim_sites(ll),:));
                plot_responses(all_out(ss).tt>-0.010 & all_out(ss).tt<0.010) = NaN;

                % scaling to unit length(https://en.wikipedia.org/wiki/Feature_scaling)
                % unti length taken in same window as stats
                response_vector_length = sum(plot_responses(all_out(ss).tt>t_win_cod(1) &  all_out(ss).tt<t_win_cod(2)).^2).^.5;
                plot_responses_norm = plot_responses./(response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial
                
                % we looked at the data and visually split out 3 stim pairs
                % that produce different latencies
                if ll<=3% First 3 pairs in parahipp,  from RB probe
                    plot(all_out(ss).tt, plot_responses_norm,'color',[0.6350 0.0780 0.1840])
                    disp(['Red, later' all_out(ss).average_ccep_names(these_stim_sites(ll))])
                elseif ll>=4 % Las pairs in parahipp, from RC probe
                    plot(all_out(ss).tt, plot_responses_norm,'color',[0.4660 0.6740 0.1880])
                    disp(['Green, early' all_out(ss).average_ccep_names(these_stim_sites(ll))])
%                 elseif ll==5 % stim pair 5 is the earliest
%                     plot(all_out(ss).tt, ss *.2 + plot_responses_norm,'color',[0.4940 0.1840 0.5560])
%                     disp(['earlier purple' all_out(ss).average_ccep_names(these_stim_sites(ll))])
                end
                
                % save outputs
                resp_counter = resp_counter+1;
                out_plot_responses_norm(resp_counter,:) = plot_responses_norm;
                out_subj_ind(resp_counter,:) = ss;
            end
        end
    end
end

title('Green=RC probe; Red=RB probe')
xlim([-0.2 .6])%, ylim([-2000 4000]);
xlabel('Time (s)'), ylabel('Normalized Voltage')

set(gcf,'PaperPositionMode','auto')
% print(fullfile('./local',figName),'-dpng');%,'-r700';
% print(fullfile('./local',figName),'-painters','-depsc')%,'-r500',)

%% Sub-03, plot temporal delay
% rgb_color = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.4660 0.6740 0.1880],[0.4940 0.1840 0.5560],[0.9290 0.6940 0.1250],[0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840]};
% blue, orange, green, purple, mustard, celeste, wine
figure('Position',[0 0 500 350]), hold on; %('Position',[0 0 600 200]), hold on

stim_area = 1; % 1 for hippocampus, 2 for amygdala, 3 for thalamus

out_plot_responses_norm = [];
out_subj_ind = [];
resp_counter = 0;
for ss = 3 % S3 with one probe
    these_measured_sites = find(all_out(ss).channel_areas==5 | all_out(ss).channel_areas==6 | all_out(ss).channel_areas==7); % posterior cing
    if stim_area==1
        these_stim_sites = find(all_out(ss).stim_area==1 | all_out(ss).stim_area==8); % hippocampal formation 
    elseif stim_area==2
        these_stim_sites = find(all_out(ss).stim_area==2); % amygdala
    elseif stim_area==3
        these_stim_sites = find(all_out(ss).stim_area==9); % thalamus
    end
  
    plot(all_out(ss).tt, zeros(size(all_out(ss).tt)),'k')
    
    sign_resp = all_out(ss).sig_ccep==1; % plot p<0.05 FDR corrected
    
    % walk through pairs, test significance and plot if so
    for kk = 1:length(these_measured_sites)
        for ll = 1:length(these_stim_sites)
            if sign_resp(these_measured_sites(kk),these_stim_sites(ll))==1
                plot_responses = squeeze(all_out(ss).average_ccep(these_measured_sites(kk),these_stim_sites(ll),:));
                plot_responses(all_out(ss).tt>-0.010 & all_out(ss).tt<0.010) = NaN;

                % scaling to L2 norm (https://en.wikipedia.org/wiki/Feature_scaling)
                % unit length taken in same window as stats
                response_vector_length = sum(plot_responses(all_out(ss).tt>t_win_cod(1) &  all_out(ss).tt<t_win_cod(2)).^2).^.5;
                plot_responses_norm = plot_responses./(response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial

                % we looked at the data and visually split out 2 stim pairs
                % that produce different latencies
                if ll<=2 % anterior sites, RB1-3
                    plot(all_out(ss).tt, plot_responses_norm,'color',[0.6350 0.0780 0.1840])
                    disp(['Red, later' all_out(ss).average_ccep_names(these_stim_sites(ll))])
                elseif ll>=3 % posterior sites, from RB3-5
                    plot(all_out(ss).tt, plot_responses_norm,'color',[0.4660 0.6740 0.1880])
                    disp(['Green, early' all_out(ss).average_ccep_names(these_stim_sites(ll))])
                end

                % save outputs
                resp_counter = resp_counter + 1;
                out_plot_responses_norm(resp_counter,:) = plot_responses_norm;
                out_subj_ind(resp_counter,:) = ss;
            end
        end
    end
end
title('Green=RB3-RB5; Red=RB1-RB3')
xlim([-0.2 .6])%, ylim([-2000 4000]);
xlabel('time (s)'),ylabel('Normilized Voltage') 

set(gcf,'PaperPositionMode','auto')
% print(fullfile('./local',figName),'-dpng');%,'-r700';
% print(fullfile('./local',figName),'-painters','-depsc')%,'-r500',)
