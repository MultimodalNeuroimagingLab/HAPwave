
clearvars, clc%, close all

%% Plot brain with estimation of electrodes, fornix & cingulum bundles 
% Dependencies: MNL ieeg basics, vistasoft, and AFQ & Along-tract statistics github repositories
% cd to HAPwave repository
addpath(genpath(pwd))

% set local path to your BIDS directory:
myPath = setLocalDataPath(1);
localDataPath = myPath.input;
 
%% load subject 02 - all tracks
sub_label = '02';
sub_hemi = 'r';
bids_path = localDataPath;

dwi_file = fullfile(bids_path,'derivatives','qsiprep',['sub-' sub_label],'ses-compact3T01','dwi',['sub-' sub_label '_ses-compact3T01_acq-diadem_space-T1w_desc-preproc_dwi.nii.gz']);
ni_dwi = niftiRead(dwi_file);

all_trks = {'Fornix_R',...
    'Cingulum_Parolfactory_R',...
    'Cingulum_Frontal_Parietal_R',...
    'Cingulum_Parahippocampal_Parietal_R',...
    'Cingulum_Parahippocampal_R'};...

fg_fromtrk = [];

for ss = 1:length(all_trks)

    tkr_name = all_trks{ss};
    trk_file = fullfile(bids_path,'derivatives','qsiprep',['sub-' sub_label],'ses-compact3T01','dwi','dsi_studio',[tkr_name '.trk']);
    [header,tracks] = trk_read(trk_file);

    % correct header
    header.vox_to_ras = ni_dwi.qto_xyz;

    transf_mat = header.vox_to_ras;
    for ii = 1:3
        transf_mat(:,ii) = transf_mat(:, ii)./header.voxel_size(ii);
    end

    fg_fromtrk(ss).name = tkr_name;
    fg_fromtrk(ss).colorRgb = [20 90 200];
    fg_fromtrk(ss).thickness = 0.5;
    fg_fromtrk(ss).visible = 1;
    fg_fromtrk(ss).seeds = [];
    fg_fromtrk(ss).seedRadius = 0;
    fg_fromtrk(ss).fibers = cell(length(tracks),1);
    for kk = 1:length(tracks)
        this_strm = transf_mat*[tracks(kk).matrix ones(length(tracks(kk).matrix),1)]';
        fg_fromtrk(ss).fibers{kk} = this_strm(1:3,:);
        clear this_strm
    end
    clear header tracks
end

electrodes_tsv = fullfile(bids_path,'derivatives','qsiprep',['sub-' sub_label],['sub-' sub_label '_ses-ieeg01_space-T1w_desc-qsiprep_electrodes.tsv']);
loc_info = readtable(electrodes_tsv,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

t1_name = fullfile(bids_path,'derivatives','qsiprep',['sub-' sub_label],'anat',['sub-' sub_label '_desc-preproc_T1w.nii.gz']);
t1 = readFileNifti(t1_name);
t1.sto_xyz = t1.qto_xyz;
t1.sto_ijk = t1.qto_ijk;

% only contacts in ROIs
hippEls = {'RB1','RB2','RB3','RB4','RB5'};
hipp_inds = find(ismember(loc_info.name,hippEls));
pcEls = {'RZ1','RY1'};
pc_inds = find(ismember(loc_info.name,pcEls));
antEls = {'ROP1'};
ant_inds = find(ismember(loc_info.name,antEls));
elecmatrix = [loc_info.x loc_info.y loc_info.z];

%% load subject 07 - all tracks
sub_label = '07';
sub_hemi = 'l';
bids_path = localDataPath;

dwi_file = fullfile(bids_path,'derivatives','qsiprep',['sub-' sub_label],'ses-compact3T01','dwi',['sub-' sub_label '_ses-compact3T01_acq-diadem_space-T1w_desc-preproc_dwi.nii.gz']);
ni_dwi = niftiRead(dwi_file);

all_trks = {'Fornix_L',...
    'Cingulum_Parolfactory_L',...
    'Cingulum_Frontal_Parietal_L',...
    'Cingulum_Parahippocampal_Parietal_L',...
    'Cingulum_Parahippocampal_L'};%

fg_fromtrk = [];

for ss = 1:length(all_trks)

    tkr_name = all_trks{ss};
    trk_file = fullfile(bids_path,'derivatives','qsiprep',['sub-' sub_label],'ses-compact3T01','dwi','dsi_studio',[tkr_name '.trk']);
    [header,tracks] = trk_read(trk_file);

    % correct header
    header.vox_to_ras = ni_dwi.qto_xyz;

    transf_mat = header.vox_to_ras;
    for ii = 1:3
        transf_mat(:,ii) = transf_mat(:, ii)./header.voxel_size(ii);
    end

    fg_fromtrk(ss).name = tkr_name;
    fg_fromtrk(ss).colorRgb = [20 90 200];
    fg_fromtrk(ss).thickness = 0.5;
    fg_fromtrk(ss).visible = 1;
    fg_fromtrk(ss).seeds = [];
    fg_fromtrk(ss).seedRadius = 0;
    fg_fromtrk(ss).fibers = cell(length(tracks),1);
    for kk = 1:length(tracks)
        this_strm = transf_mat*[tracks(kk).matrix ones(length(tracks(kk).matrix),1)]';
        fg_fromtrk(ss).fibers{kk} = this_strm(1:3,:);
        clear this_strm
    end
    clear header tracks
end

electrodes_tsv = fullfile(bids_path,'derivatives','qsiprep',['sub-' sub_label],['sub-' sub_label '_ses-ieeg01_space-T1w_desc-qsiprep_electrodes.tsv']);
loc_info = readtable(electrodes_tsv,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

t1_name = fullfile(bids_path,'derivatives','qsiprep',['sub-' sub_label],'anat',['sub-' sub_label '_desc-preproc_T1w.nii.gz']);
t1 = readFileNifti(t1_name);
t1.sto_xyz = t1.qto_xyz;
t1.sto_ijk = t1.qto_ijk;

% only contacts in ROIs
hippEls = {'LB1','LB2','LB3','LB4','LC1','LC2'};
hipp_inds = find(ismember(loc_info.name,hippEls));
pcEls = {'LZ1','LZ4'};
pc_inds = find(ismember(loc_info.name,pcEls));
antEls = {'LK1'};
ant_inds = find(ismember(loc_info.name,antEls));
elecmatrix = [loc_info.x loc_info.y loc_info.z];

%% Render + electrodes
% Assign colors
trk_color = {[0.4660 0.6740 0.1880],[0 0.4470 0.7410],[0.9290 0.6940 0.1250],[0.8500 0.3250 0.0980],[0.7000 0.5000 0.7000],[0.6350 0.0780 0.1840],[0.4940 0.1840 0.5560]}; % color codes
trk_color_label = {'green','blue','yellow','orange','pink','red','purple'}; % color labels

if isequal(sub_hemi,'l')
    h_rotation = 270;
    h_plane = -.5;
    h_light = 'left';
elseif isequal(sub_hemi,'r')
     h_rotation = 90;
     h_plane = .5;
     h_light = 'right';
end

% fornix 
AFQ_RenderFibers(fg_fromtrk(1),'numfibers',500,'color',trk_color{1},'alpha',.2);% light yellow
% cingulate
AFQ_RenderFibers(fg_fromtrk(2),'numfibers',700,'color',trk_color{2},'alpha',.15,'newfig',false);
AFQ_RenderFibers(fg_fromtrk(3),'numfibers',300,'color',trk_color{3},'alpha',.1,'newfig',false); 
AFQ_RenderFibers(fg_fromtrk(4),'numfibers',300,'color',trk_color{4},'alpha',.1,'newfig',false);
AFQ_RenderFibers(fg_fromtrk(5),'numfibers',200,'color',trk_color{4},'alpha',.1,'newfig',false);

% Select ROIs:
ieeg_elAdd(elecmatrix(ismember(loc_info.name,hippEls),:),'w',50)    % Hipp contacts
ieeg_elAdd(elecmatrix(ismember(loc_info.name,pcEls),:),'w',50)      % PCC contacts
ieeg_elAdd(elecmatrix(ismember(loc_info.name,antEls),:),'w',50)     % ANT contacts
% ieeg_elAdd(elecmatrix,[.99 .99 .99],20) % plot all electrodes

% Set the sagittal plane to plot
AFQ_AddImageTo3dPlot(t1,[(h_plane), 0, 0],[],0,[],[0 3000]) %  R = +, L = -

view(h_rotation,30)
axis image
% camlight(h_light);

%% save figure
imName = fullfile(bids_path,'derivatives','matlabOut','dsi',['sub-' sub_label '_ses-compact3T01_dMRIwEls.png']);
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',imName)

