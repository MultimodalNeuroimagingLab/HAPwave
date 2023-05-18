clearvars, clc, close all
% startup

%% Load stats across subs
% cd to pcc_project repository
addpath(genpath(pwd))
% set local path to your BIDs directory:

% load the meta data
all_subjects = {'01','02','03','04','05','06','07','08'}; % 
all_hemi = {'r','r','r','l','r','l','l','r'};
all_runs = {'01','01','01','01','01','01','01','01'};

t_win_cod = [0.015 .5];

for ss = 1:7%length(all_subjects)
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
    electrodes_table = readtable(electrodes_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'}); % electrodes table
    % find good sEEG channels
    good_channels = find(ismember(channels_table.type,{'ECOG','SEEG'}) & ismember(channels_table.status,'good'));

    % load stats/average CCEPs from outputName
    statsFile = fullfile(localDataPath,'derivatives','stats',['sub-' bids_sub],...
        ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_stats.mat']);
    load(statsFile,'average_ccep','average_ccep_names','tt','srate','cross_ccep_t','cross_ccep_p');

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

%% Get significant CCEPs across subs 

stim_names = {'Hipp','Amyg','Thal'};
stim_area = 2; % 1 for hippocampus, 2 for amygdala, 3 for thalamus
rec_names = {'PCC','ACC'};
rec_area = 1; % 1 for PCC, 2 for ACC

out_plot_responses_norm = [];
out_subj_ind = [];
resp_counter = 0;

subj_resp_total = zeros(7,1); % stim-->measured pair
subj_resp_sign = zeros(7,1); % stim-->measured pair significant

for ss = 1:7 
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
    if ~isempty(these_measured_sites)
        % loop over measured sites
        for kk = 1:length(these_measured_sites)
            % loop over the stimulated pairs
            for ll = 1:length(these_stim_sites)
                subj_resp_total(ss) = subj_resp_total(ss) + 1;
                % test significance
                if sign_resp(these_measured_sites(kk), these_stim_sites(ll)) == 1
                   subj_resp_sign(ss) = subj_resp_sign(ss) + 1;
                   
                   plot_responses = squeeze(all_out(ss).average_ccep(these_measured_sites(kk), these_stim_sites(ll), :));
                   plot_responses(all_out(ss).tt > -0.010 & all_out(ss).tt < 0.010) = NaN;

%                     % Normalize. Scaling to unit length (Euclidean lenght): https://en.wikipedia.org/wiki/Feature_scaling
%                     % unit length taken in same window as stats
                    response_vector_length = sum(plot_responses(all_out(ss).tt > t_win_cod(1) &  all_out(ss).tt < t_win_cod(2)) .^ 2) .^ .5;
                    plot_responses_norm = plot_responses ./ (response_vector_length*ones(size(plot_responses))); % normalize (L2 norm) each trial

                    % save outputs
                    resp_counter = resp_counter + 1;
                    out_plot_responses_norm(resp_counter, :) = plot_responses_norm;
                    out_subj_ind(resp_counter, :) = ss;

                end % done testing significance
            end % done looping through stim pairs
        end  % done looping through stim pairs
    end
end

%% create input to MDS across subjects
V_allsubs = [];
labels_allsubs = [];
for ss = 1:7

    tt = all_out(ss).tt;
%     t_win = find(tt>0.015 & tt<0.600);
    t_win = find(tt>0.000 & tt<1.500);
    good_channels = all_out(ss).good_channels;
    
    all_data = all_out(ss).average_ccep(good_channels,:,t_win);
    X = reshape(all_data,size(all_data,1)*size(all_data,2),size(all_data,3));

    % do the same to channel and stim area labels
    conn_label = cat(3,repmat(all_out(ss).channel_areas(good_channels),1,length(all_out(ss).stim_area),1),...
        repmat(all_out(ss).stim_area',length(all_out(ss).channel_areas(good_channels)),1));
    % column 1 now has channel areas, column 2 has stim areas
    conn_label_R = reshape(conn_label,size(conn_label,1)*size(conn_label,2),2);

    % remove empty rows from both
    conn_label_R(isnan(X(:,1)),:) = [];
    X(isnan(X(:,1)),:) = [];

    % only use connections in limbic network
    use_these = find(conn_label_R(:,1)>0 & conn_label_R(:,2)>0);
    XX = X(use_these,:);
    use_these_labels = conn_label_R(use_these,:);
    
    V_allsubs = [V_allsubs; XX];
    labels_allsubs = [labels_allsubs; use_these_labels];
end


%% plot group level MDS
xx_allsubs = V_allsubs;

% normalization does not change much here, since done later
vector_length = sum(xx_allsubs.^ 2,2) .^ .5;
xx_allsubs = xx_allsubs ./ repmat(vector_length,1,size(xx_allsubs,2)); % normalize (L2 norm) each trial

% xx_allsubs = abs(xx_allsubs);

% construct distance matrix
D_length = sqrt(sum(xx_allsubs.^2,2));
D_dot = (xx_allsubs*xx_allsubs');
D = D_dot./(D_length*D_length'); 
D(D>1) = 1;
D_mat = rad2deg(acos(D));
[Y,eigvals] = cmdscale(round(D_mat),2);
% this is the same as the cosine distance:
% D_mat = pdist(xx_allsubs,'spearman'); % cosine or spearman
% [Y,eigvals] = cmdscale((D_mat),2);

% D_mat = pdist(xx_allsubs,'cosine'); % may like this, or cityblock?
% [Y,eigvals] = cmdscale((D_mat),2);

% figure,plot(eigvals)

% D_mat = pdist(xx_allsubs,'correlation');
% [Y] = mdscale(round(D_mat*100),2); % nonclassical

figure('Position',[0 0 800 300]),hold on 

subplot(1,2,1),hold on
hipp2pcc = (ismember(labels_allsubs(:,1),[5 6 7]) & ismember(labels_allsubs(:,2),[1 8]));
h1 = scatter(Y(hipp2pcc,1),Y(hipp2pcc,2),80,[.2 .2 1],'filled','MarkerFaceAlpha',.5); 

hipp2acc = (ismember(labels_allsubs(:,1),[3 4]) & ismember(labels_allsubs(:,2),[1 8]));
h2 = scatter(Y(hipp2acc,1),Y(hipp2acc,2),50,[1 1 0],'filled','MarkerFaceAlpha',.5); 

amyg2acc = (ismember(labels_allsubs(:,1),[3 4]) & ismember(labels_allsubs(:,2),[2]));
h3 = scatter(Y(amyg2acc,1),Y(amyg2acc,2),80,[1 .3 0],'filled','MarkerFaceAlpha',.5); 

ant2pcc = (ismember(labels_allsubs(:,1),[5 6 7]) & ismember(labels_allsubs(:,2),[9]));
h4 = scatter(Y(ant2pcc,1),Y(ant2pcc,2),80,[0 1 1],'filled','MarkerFaceAlpha',.8); 

amyg2pcc = (ismember(labels_allsubs(:,1),[5 6 7]) & ismember(labels_allsubs(:,2),[2]));
h5 = scatter(Y(amyg2pcc,1),Y(amyg2pcc,2),50,[0 1 .1],'filled','MarkerFaceAlpha',.5); 

scatter(Y(~hipp2pcc & ~amyg2pcc & ~hipp2acc & ~amyg2acc & ~ant2pcc,1),Y(~hipp2pcc & ~hipp2acc & ~amyg2pcc & ~amyg2acc & ~ant2pcc,2),...
    10,[.5 .5 .5],'filled','MarkerFaceAlpha',.5); 


subplot(1,2,2),hold on
% An ellipse can be defined as the locus of all points that satisfy the equations
% x = a cos t % y = b sin t
% where:
% x,y are the coordinates of any point on the ellipse, % a, b are the radius on the x and y axes respectively,
t = linspace(0,2*pi) ;

hipp2pcc = (ismember(labels_allsubs(:,1),[5 6 7]) & ismember(labels_allsubs(:,2),[1 8]));
cc = mean(Y(hipp2pcc,:),1);
xy_std = std(Y(hipp2pcc,:),[],1);
x = xy_std(1)*cos(t) ; y = xy_std(2)*sin(t) ;
plot(cc(1)+x,cc(2)+y,'r')

hipp2acc = (ismember(labels_allsubs(:,1),[3 4]) & ismember(labels_allsubs(:,2),[1 8]));
cc = mean(Y(hipp2acc,:),1);
xy_std = std(Y(hipp2acc,:),[],1);
x = xy_std(1)*cos(t) ; y = xy_std(2)*sin(t) ;
plot(cc(1)+x,cc(2)+y,'b')

amyg2acc = (ismember(labels_allsubs(:,1),[3 4]) & ismember(labels_allsubs(:,2),[2]));
cc = mean(Y(amyg2acc,:),1);
xy_std = std(Y(amyg2acc,:),[],1);
x = xy_std(1)*cos(t) ; y = xy_std(2)*sin(t) ;
plot(cc(1)+x,cc(2)+y,'Color',[0 .3 0])

ant2pcc = (ismember(labels_allsubs(:,1),[5 6 7]) & ismember(labels_allsubs(:,2),[9]));
cc = mean(Y(ant2pcc,:),1);
xy_std = std(Y(ant2pcc,:),[],1);
x = xy_std(1)*cos(t) ; y = xy_std(2)*sin(t) ;
plot(cc(1)+x,cc(2)+y,'m')

amyg2pcc = (ismember(labels_allsubs(:,1),[5 6 7]) & ismember(labels_allsubs(:,2),[2]));
cc = mean(Y(amyg2pcc,:),1);
xy_std = std(Y(amyg2pcc,:),[],1);
x = xy_std(1)*cos(t) ; y = xy_std(2)*sin(t) ;
plot(cc(1)+x,cc(2)+y,'c')


axis equal