
clearvars, clc, close all

%% Plot brain with estimation of electrodes, fornix & cingulum bundles 
% DEPENDENCIES: Vistasoft, AFQ, and Along-tract statistics github repositories
% cd to HAPwave repository
addpath(genpath(pwd))

% set local path to your BIDS directory:
myPath          = setLocalDataPath(1);
localDataPath   = myPath.input;
 
%% load subject - all tracks
% Subjects 2 and 7 have dMRI data
sub_label   = '02'; 
% bids_path   = localDataPath;

ss = sub_label;

switch ss
    case {'02'}
        % set hemisphere
        sub_hemi    = 'r'; 
        % set contacts in ROIs
        hippEls     = {'RB1','RB2','RB3','RB4','RB5'};
        pcEls       = {'RZ1','RY1'};
        antEls      = {'ROP1'};
        % set name of bundles
        all_trks    = {'Fornix_R',...
                       'Cingulum_Parolfactory_R',...
                       'Cingulum_Frontal_Parietal_R',...
                       'Cingulum_Parahippocampal_Parietal_R',...
                       'Cingulum_Parahippocampal_R'};...
    case {'07'}
        % set hemisphere
        sub_hemi    = 'l';
        % only contacts in ROIs
        hippEls     = {'LB1','LB2','LB3','LB4','LC1','LC2'};
        pcEls       = {'LZ1','LZ4'};
        antEls      = {'LK1'};
        % set name of bundles
        all_trks    = {'Fornix_L',...
                       'Cingulum_Parolfactory_L',...
                       'Cingulum_Frontal_Parietal_L',...
                       'Cingulum_Parahippocampal_Parietal_L',...
                       'Cingulum_Parahippocampal_L'};%
end

% load dMRI-brain
dwi_file    = fullfile(localDataPath,'derivatives','dMRI',['sub-' sub_label],'dwi',['sub-' sub_label '_ses-compact3T01_acq-diadem_space-T1w_desc-preproc_dwi.nii.gz']);
ni_dwi      = niftiRead(dwi_file);
fg_fromtrk  = [];

for jj          = 1:length(all_trks)

    tkr_name        = all_trks{jj};
    trk_file        = fullfile(localDataPath,'derivatives','dMRI',['sub-' sub_label],'dwi','trks',[tkr_name '.trk']);
    [header,tracks] = trk_read(trk_file);

    % correct header
    header.vox_to_ras = ni_dwi.qto_xyz;

    transf_mat      = header.vox_to_ras;
    for ii = 1:3
        transf_mat(:,ii) = transf_mat(:, ii)./header.voxel_size(ii);
    end

    fg_fromtrk(jj).name         = tkr_name;
    fg_fromtrk(jj).colorRgb     = [20 90 200];
    fg_fromtrk(jj).thickness    = 0.5;
    fg_fromtrk(jj).visible      = 1;
    fg_fromtrk(jj).seeds        = [];
    fg_fromtrk(jj).seedRadius   = 0;
    fg_fromtrk(jj).fibers       = cell(length(tracks),1);
    for kk = 1:length(tracks)
        this_strm = transf_mat*[tracks(kk).matrix ones(length(tracks(kk).matrix),1)]';
        fg_fromtrk(jj).fibers{kk} = this_strm(1:3,:);
        clear this_strm
    end
    clear header tracks
end
% load electrodes in dMRI-brain space
electrodes_tsv  = fullfile(localDataPath,'derivatives','dMRI',['sub-' sub_label],['sub-' sub_label '_ses-ieeg01_space-T1w_desc-qsiprep_electrodes.tsv']);
loc_info    = readtable(electrodes_tsv,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

% get XYZs of ROI
hipp_inds   = find(ismember(loc_info.name,hippEls));
pc_inds     = find(ismember(loc_info.name,pcEls));
ant_inds    = find(ismember(loc_info.name,antEls));
elecmatrix  = [loc_info.x loc_info.y loc_info.z];

t1_name     = fullfile(localDataPath,'derivatives','dMRI',['sub-' sub_label],'anat',['sub-' sub_label '_desc-preproc_T1w.nii.gz']);
t1          = readFileNifti(t1_name);
t1.sto_xyz  = t1.qto_xyz;
t1.sto_ijk  = t1.qto_ijk;

    
%% Render + electrodes
% Assign colors to use in bondles
trk_color   = {[0.466 0.674 0.188],... % green
               [0.000 0.447 0.741],... % blue
               [0.929 0.694 0.125],... % yellow
               [0.850 0.325 0.098],... % orange
               [0.700 0.500 0.700],... % pink
               [0.635 0.078 0.184],... % red
               [0.494 0.184 0.556]};   % purple

% set some plotting params depending on hemi
if isequal(sub_hemi,'l')
    h_rotation  = 270;
    h_plane     = -.5;
    h_light     = 'left';
elseif isequal(sub_hemi,'r')
     h_rotation = 90;
     h_plane    = .5;
     h_light    = 'right';
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
AFQ_AddImageTo3dPlot(t1,[(h_plane), 0, 0],[],0,[],[]) %  R = +, L = -

view(h_rotation,30)
axis image
camlight(h_light);
