clearvars, clc, close all
% startup

%% Get native XYZ into MNI space
% GOV, DH, HH 2022 Multimodal Neuroimaging Lab
% Dependencies: vistasoft and spm12 github repositories

% % cd to HAPwave repository
addpath(genpath(pwd))

% set local path to your BIDS directory:
myPath = setLocalDataPath(1);
localDataPath = myPath.input;
outDataPath = myPath.output;

%% convert all subjects and save MNI positions 
% load the meta data
all_subjects    = {'01','02','03','04','05','06','07','08'}; 
all_hemi        = {'r','r','r','l','r','l','l','r'};

for ss = 1:length(all_subjects)
    bids_sub        = all_subjects{ss};
    bids_ses        = 'ieeg01';
    bids_task       = 'ccep';
    bids_run        = '01';

    % set .nii and channels, & electrodes filename
    niiPath             = fullfile(localDataPath,'derivatives','freesurfer',['sub-' bids_sub],['sub-' bids_sub '_ses-mri01_T1w_deFaced.nii']);
    electrodes_tsv_name = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
                                   ['sub-' bids_sub '_ses-' bids_ses '_electrodes.tsv']);
    channels_tsv_name   = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
                                   ['sub-' bids_sub '_ses-' bids_ses '_task-' bids_task '_run-' bids_run '_channels.tsv']);
    % load channels & electrodes
    channels_table      = readtable(channels_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'}); % channels table
    electrodes          = readtable(electrodes_tsv_name, 'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    xyz                 = [electrodes.x, electrodes.y, electrodes.z];

    % Set where electrode nifti images are saved to & read from
    spmT1wseg           = fullfile(localDataPath,'derivatives','MNI',['sub-' bids_sub]);

    % run the function to get the xyz positions in MNI 
    % y_spmOutputFile.nii MUST be in the same folder as nativeT1.nii
    xyzMni              = ieeg_getXyzMni(xyz, niiPath, spmT1wseg);
    % find XYZ coordinates in new table
    t_new               = electrodes;
    t_new.x             = xyzMni(:,1);
    t_new.y             = xyzMni(:,2);
    t_new.z             = xyzMni(:,3);

    % write as a new space-MNI152_electrodes.tsv
    [a,b,c]             = fileparts(electrodes_tsv_name);
    mni_elec_tsv_name   = fullfile(localDataPath,'derivatives','MNI',['sub-' bids_sub],['sub-' bids_sub '_ses-' bids_ses '_space-MNI152NLin6Sym_electrodes.tsv']);
    writetable(t_new, mni_elec_tsv_name, 'FileType','text','Delimiter','\t');

    % save in a nifti to check
    mni_templateImage   = fullfile(localDataPath,'derivatives','MNI','mni152.nii');
    ieeg_position2reslicedImage(xyzMni,mni_templateImage)

end

%% Concatenate XYZs across subjects
all_subjects    = {'01','02','03','04','05','06','07','08'}; % 
mni_elec        = table();
MNI_elec_xsub   = fullfile(localDataPath,'derivatives','MNI','sub-all_ses-ieeg01_space-MNI152NLin6Sym_electrodes.tsv');

for ss = 1:length(all_subjects)
    bids_sub        = all_subjects{ss};
    bids_ses        = 'ieeg01';
    bids_task       = 'ccep';
    bids_run        = '01';

    % get electrodes filename
    mni_elec_tsv    = fullfile(localDataPath,'derivatives','MNI', ['sub-' bids_sub],['sub-' bids_sub '_ses-' bids_ses '_space-MNI152NLin6Sym_electrodes.tsv']);
    electrodes      = readtable(mni_elec_tsv, 'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
   
    % add column with subject's code
    ncol            = height(electrodes);
    electrodes.sub_code = repmat({bids_sub}, ncol, 1);
    
    if ~iscell(electrodes.seizure_zone) 
        electrodes.seizure_zone = num2cell(electrodes.seizure_zone);
    end

    % concatenate the subjects
    mni_elec        = vertcat(mni_elec,electrodes);
end
mni_elec        = bids_tsv_nan2na(mni_elec);
writetable(mni_elec, MNI_elec_xsub, 'FileType','text','Delimiter','\t');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IF ONLY PLOTTING start here
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare limbic's XYZs
% set variables
all_subjects    = {'01','02','03','04','05','06','07','08'}; % 
MNI_elec_xsub   = fullfile(localDataPath,'derivatives','MNI','sub-all_ses-ieeg01_space-MNI152NLin6Sym_electrodes.tsv');
mni_elec        = readtable(MNI_elec_xsub,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
mni_elec        = bids_tsv_nan2na(mni_elec); 

% get MNI area number and xyz, easiest if these are number vectors
mni_area        = NaN(size(mni_elec,1),1);
mni_sub         = NaN(size(mni_elec,1),1);
elecmatrix      = NaN(size(mni_elec,1),3);

for kk = 1:size(mni_elec,1)
    if isnumeric(mni_elec.Destrieux_label{kk}) % should not have n/a
        mni_sub(kk)     = mni_elec.sub_code{kk};
        mni_area(kk)    = mni_elec.Destrieux_label{kk};
        elecmatrix(kk,:)= [mni_elec.x{kk} mni_elec.y{kk} mni_elec.z{kk}];

    end
end


%% Render MNI brain with electrodes
% Sort limbic ROIs by hemisphere
limbic_code_L   = [17 11123 18 11106 11107 11108 11109 11110 10];
limbic_code_R   = [53 12123 54 12106 12107 12108 12109 12110 49];

% %% Select hemisphere to render%%
hemi = 'l';


% get render of selected hemisphere
if isequal(hemi,'l')
    g = gifti(fullfile(localDataPath, 'derivatives','MNI','MNI152_left.gii'));
    ldegree = 97; % change viewing angle
    area_code = limbic_code_L;
elseif isequal(hemi,'r')
    g = gifti(fullfile(localDataPath, 'derivatives','MNI','MNI152_right.gii'));
    ldegree = 270; % change viewing angle
    area_code = limbic_code_R;
end
figure
tH = ieeg_RenderGifti(g);
ieeg_viewLight(ldegree,0) % change viewing angle   
set(tH,'FaceAlpha',.40) % make transparent


% select the electrode positions and codes of the sites we want to plot
ROI_inds = find(ismember(mni_area,area_code)); % find gives indices
el_xyz = elecmatrix(ROI_inds,:);            % electrode coordinates
el_codes = mni_area(ROI_inds);              % Destrieux code
el_groups = zeros(size(el_codes));          % initiate here. 1=Hip, 2=Amyg, 3=ACC, 4=PCC, 5=Thal
sub_codes = mni_elec.sub_code(ROI_inds);    % subject labels

limbic_grp_name = {'Hip' 'Amg' 'ACC' 'PCC' 'ANT'};
limbic_grp = {[1 2] [3] [4 5] [6 7 8] [9]}; % color for each area
limbic_color = {[0.8500 0.3250 0.0980],...  orange
                [0.4940 0.1840 0.5560],...  purple
                [0.9290 0.6940 0.1250],...  yellow
                [0 0.4470 0.7410],...       blue
                [0.4660 0.6740 0.1880]}; %  green                                                                                    % 6celeste, 7wine 

msize   = 7;
eledge  = [.2 .2 .2];
el_symb = {'^' 'o' 'p' 'd' 'v' 's' 'h' '>'}; % symbol for each subject

% get electrodes gorups, such that we can easily assign a symbol later
for kk = 1:length(ROI_inds)     % Loop through all limbic electrodes in hemisphere
    % get index into limbic grouping 
    index_into_area_codes = find(ismember(area_code,el_codes(kk)));
    for ll = 1:length(limbic_grp)
        if ismember(index_into_area_codes,limbic_grp{ll})
            el_groups(kk) = ll;
        end
    end
end

for kk = 1:length(ROI_inds) % Loop through all limbic electrodes in hemisphere
    ieeg_eladd_symbol(el_xyz(kk,:),limbic_color{el_groups(kk)},eledge,msize,'o')
%     % Use a different symbol per subject
%     ieeg_eladd_symbol(el_xyz(kk,:),limbic_color{el_groups(kk)},eledge,msize,el_symb{sub_codes{kk}})
end


%%
% filename
figureName = fullfile(outDataPath,['fig2_hemi_' hemi '_coverageMNI']);
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',figureName)
% close all
