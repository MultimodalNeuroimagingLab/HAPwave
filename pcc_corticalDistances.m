
% This is an example script where we calculate the distance for electrodes
% to the closest point on the pial and gray/white matter surface.

%% Dependencies: MNL ieeg basics & vistasoft github repositories
addpath(genpath(pwd))

% set local path to your BIDS directory:
myPath          = setLocalDataPath(1);
localDataPath   = myPath.input;

bids_subjects   = {'01','02','03','04','05','06','07','08'};

for ss = 1:length(bids_subjects)
    bids_sub    = bids_subjects{ss};
    bids_ses    = 'ieeg01';

    if isequal(bids_sub,'01')
        elsPCC_set = {'RY1','RY2','RQ1','RZ1'}; % super, deep, super, deep 
        elsACC_set = {'RX1','RX2','RN1','RN2'};
    elseif isequal(bids_sub,'02')
        elsPCC_set = {'RZ1'};
        elsACC_set = {'RX1','RX2','RX3','RX4','RF1','RF2','RF3','RF4','RF5'};
    elseif isequal(bids_sub,'03')
        elsPCC_set = {'RZ1','RZ2'};% super, deep
        elsACC_set = {'RAF2','RFM1','RFM2','RFM3','RFM4','RX1','RN2','RY1','RY2','RY3','RY4','RE2','RE3','RE4',};
    elseif isequal(bids_sub,'04')
        elsPCC_set = {'LY1','LZ1','LZ2'}; 
        elsACC_set = {'LX1','LX2'};
    elseif isequal(bids_sub,'05')
        elsPCC_set = {'RY1','RZ1'}; % super, deep
        elsACC_set = {}; %no ACC coverage
    elseif isequal(bids_sub,'06')
        elsPCC_set = {'LZ1','LQ1'}; 
        elsACC_set = {}; %no ACC coverage
    elseif isequal(bids_sub,'07')
        elsPCC_set = {'LZ1','LZ4'}; 
        elsACC_set = {'LF3','LF4', 'LO2','LO3','LX1','LX2','LE1','LE2','LN1, LN3','LN4','LY1'};
    elseif isequal(bids_sub,'08')
        elsPCC_set = {'RZ1','RPO1','RQ1'};
        elsACC_set = {'RX1','RX2','RF6','RF7'};
    end    
    
    % get electrode positions matrix
    electrodes_tsv_name = fullfile(localDataPath,['sub-' bids_sub],['ses-' bids_ses],'ieeg',...
        ['sub-' bids_sub '_ses-' bids_ses '_electrodes.tsv']);
    loc_info            = readtable(electrodes_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    
    % load pial and gray/white surfaces
    gL_pial     = gifti(fullfile(localDataPath,'derivatives','freesurfer',['sub-' bids_sub],'pial.L.surf.gii'));
    gR_pial     = gifti(fullfile(localDataPath,'derivatives','freesurfer',['sub-' bids_sub],'pial.R.surf.gii'));
    gL_white    = gifti(fullfile(localDataPath,'derivatives','freesurfer',['sub-' bids_sub],'white.L.surf.gii'));
    gR_white    = gifti(fullfile(localDataPath,'derivatives','freesurfer',['sub-' bids_sub],'white.R.surf.gii'));
    
    distLim     = 6; % set distance limit
    [out]       = ieeg_eldist2pial2white(loc_info,gL_pial,gR_pial,gL_white,gR_white, distLim);

    for kk = 1:length(els_set)
        this_ind = find(ismember(out.name,els_set{kk}));
        disp([bids_sub ' ' els_set{kk} ' dist el-white ' num2str(out.rel_dist(this_ind)) ' dist pial-white '  num2str(out.dist_pialwhite(this_ind))])
    end
    
end
