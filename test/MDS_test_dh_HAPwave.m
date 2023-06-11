
% Run anova code first to load data

%% create input to MDS across subjects

% area_names = {'Hipp','Amyg','PCC','ACC'};   
% out(measure_ind,stim_ind)

tt_int = [0.015 0.250]; %reasonable whether normalized or non-normalized
% with spearman
% tt_int = [0.015 0.250]; %reasonable whether normalized or non-normalized
% with spearman

set1 = out(3,1).plot_responses_norm;
sign_reversal = out(3,1).elec_relDist>2.3;
set1(sign_reversal==1,:) = -1*set1(sign_reversal==1,:);
% set1 = diff(set1,[],2);

set2 = out(4,2).plot_responses_norm;
sign_reversal = out(4,2).elec_relDist>2.3;
set2(sign_reversal==1,:) = -1*set2(sign_reversal==1,:);
% set2 = diff(set2,[],2);

% V non-normalized
% V_allsubs = [out(3,1).plot_responses(out(3,1).p<0.05,tt>tt_int(1) & tt<tt_int(2)); ...
%     out(4,2).plot_responses(out(4,2).p<0.05,tt>tt_int(1) & tt<tt_int(2))];
% 
% V normalized 
V_allsubs = [set1(out(3,1).p<0.05,tt>tt_int(1) & tt<tt_int(2)); ...
    set2(out(4,2).p<0.05,tt>tt_int(1) & tt<tt_int(2))];

% C
% V_allsubs = [out(3,1).C(out(3,1).p<0.05,tt>tt_int(1) & tt<tt_int(2)); ...
%     out(4,2).C(out(4,2).p<0.05,tt>tt_int(1) & tt<tt_int(2))];

[u,s,v] = svd(abs(V_allsubs));

labels_allsubs = [zeros(size(out(3,1).subj_ind(out(3,1).p<0.05) )); ...
    1+zeros(size(out(4,2).subj_ind(out(4,2).p<0.05)))];
all_sub_nrs = [out(3,1).subj_ind(out(3,1).p<0.05); ...
    out(4,2).subj_ind(out(4,2).p<0.05)];

%% plot group level MDS
xx_allsubs = V_allsubs;

% normalization does not change much here, since done later
% vector_length = sum(xx_allsubs.^ 2,2) .^ .5;
% xx_allsubs = xx_allsubs ./ repmat(vector_length,1,size(xx_allsubs,2)); % normalize (L2 norm) each trial

% xx_allsubs = abs(xx_allsubs);

% % construct distance matrix
% D_length = sqrt(sum(xx_allsubs.^2,2));
% D_dot = (xx_allsubs*xx_allsubs');
% D = D_dot./(D_length*D_length'); 
% D(D>1) = 1;
% D_mat = rad2deg(acos(D));
% [Y,eigvals] = cmdscale(round(D_mat),2);

% this is the same as the cosine distance:

% D_mat = pdist(xx_allsubs,'spearman'); % cosine or spearman
% [Y,eigvals] = cmdscale((D_mat),2);

D_mat = pdist(xx_allsubs,'cosine'); % may. like this, or cityblock?
[Y,eigvals] = cmdscale((D_mat),2);

% figure,plot(eigvals)

% D_mat = pdist(xx_allsubs,'correlation');
% [Y] = mdscale(round(D_mat*100),6); % nonclassical

%%
sub_symbol = {'>','o','<','d','*','x','+','^','s'};

figure('Position',[0 0 800 300]),hold on 

subplot(1,2,2),hold on
for kk = 1:8
    hipp2pcc = ismember(labels_allsubs(:,1),[0]) & all_sub_nrs==kk;
    h1 = scatter(Y(hipp2pcc,1),Y(hipp2pcc,2),80,[.2 .2 1],'Marker',sub_symbol{kk}); 
%     h1 = scatter(Y(hipp2pcc,1),Y(hipp2pcc,2),80,[.2 .2 1],'filled','MarkerFaceAlpha',.5,'Marker',sub_symbol{kk}); 
end
for kk = 1:8
    amg2acc = ismember(labels_allsubs(:,1),[1]) & all_sub_nrs==kk;
    h2 = scatter(Y(amg2acc,1),Y(amg2acc,2),80,[1 .2 .2],'Marker',sub_symbol{kk}); 
%     h2 = scatter(Y(amg2acc,1),Y(amg2acc,2),80,[1 .2 .2],'filled','MarkerFaceAlpha',.5,'Marker',sub_symbol{kk}); 
end

%% Try zero padding after CRP duration


figure,hold on
hipp2pcc = ismember(labels_allsubs(:,1),[0]);
plot3(Y(hipp2pcc,1),Y(hipp2pcc,2),Y(hipp2pcc,3),'.','MarkerSize',20,'Color',[.2 .2 1]); 

amg2acc = ismember(labels_allsubs(:,1),[1]);
plot3(Y(amg2acc,1),Y(amg2acc,2),Y(amg2acc,3),'.','MarkerSize',20,'Color',[1 .2 .2]);


%%
%%
%% Try wavelet
%%
%%

tt_int = [0.100 0.800]; %reasonable whether normalized or non-normalized
% tt_int = [0.100 0.500]; %reasonable whether normalized or non-normalized
tt_v = tt(tt>tt_int(1) & tt<tt_int(2));

sign_responses = out(3,1).p<0.05;
V1 = out(3,1).plot_responses_norm(sign_responses>0,tt>tt_int(1) & tt<tt_int(2))';
sign_responses = out(4,2).p<0.05;
V2 = out(4,2).plot_responses_norm(sign_responses>0,tt>tt_int(1) & tt<tt_int(2))';

V = cat(2,V1,V2);

labels_allsubs = [zeros(size(V1,2),1); ...
    1+zeros(size(V2,2),1)];
all_sub_nrs = [out(3,1).subj_ind(out(3,1).p<0.05); ...
    out(4,2).subj_ind(out(4,2).p<0.05)];


srate = 2048;

[S, f] = ieeg_getWaveletSpectrogram(V, srate,  [1, 200]);

S_mean = mean(mean(S,3),2);
S = S./repmat(S_mean,1,size(S,2),size(S,3)) - 1;
% S = log10(S) - log10(repmat(S_mean,1,size(S,2),size(S,3)));

%%

load ./loc_colormap.mat

figure
subplot(3,2,1)
imagesc(tt_v,f,mean(S(:,:,labels_allsubs==0),3),[-20 20])
axis xy
title('HC to PCC')
ylabel('Frequency (Hz)')

subplot(3,2,2)
imagesc(tt_v,f,mean(S(:,:,labels_allsubs==1),3),[-20 20])
axis xy
title('AMG to ACC')

this_f = find(f>100 & f<150);
% this_f = find(f>50 & f<200);

%%
subplot(3,2,3)
title('power from 100-150 Hz')
imagesc(tt_v,[],squeeze(mean(S(this_f,:,labels_allsubs==0),1))',[-10 10])
xlabel('Time(s)')
ylabel('CCEP nr')

subplot(3,2,4)
imagesc(tt_v,[],squeeze(mean(S(this_f,:,labels_allsubs==1),1))',[-10 10])

xlabel('Time(s)')

% figure,imagesc(tt_v,f,squeeze(S_shift(10,:,:)))

% shift dimensions for spectral SVD
S_shift = shiftdim(S,2);
S_cat = reshape(S_shift,size(S_shift,1),size(S_shift,2) * size(S_shift,3));
[u,s,v] = svd(S_cat,'econ');

colormap(cm)

% put spectral SVD pattern back in matrix
spec_pat1 = zeros(size(S_shift,2),size(S_shift,3));
spec_pat1(:) = v(:,1);
spec_pat2 = zeros(size(S_shift,2),size(S_shift,3));
spec_pat2(:) = v(:,2);
spec_pat3 = zeros(size(S_shift,2),size(S_shift,3));
spec_pat3(:) = v(:,3);

subplot(3,4,9)
imagesc(tt_v,f,spec_pat1,[-0.05 0.05])
axis xy
ylabel('Frequency (Hz)')
xlabel('Time(s)')
title('PC1')

subplot(3,4,10)
imagesc(tt_v,f,spec_pat2,[-0.05 0.05])
axis xy
xlabel('Time(s)')
title('PC2')

subplot(3,4,11)
imagesc(tt_v,f,spec_pat3,[-0.05 0.05])
axis xy
xlabel('Time(s)')
title('PC3')

subplot(3,4,12),hold on
hipp2pcc = ismember(labels_allsubs(:,1),[0]);
amg2acc = ismember(labels_allsubs(:,1),[1]);
plot(u(hipp2pcc,1),u(hipp2pcc,2),'r.')
plot(u(amg2acc,1),u(amg2acc,2),'b.')
% xlabel('PC1'),ylabel('PC2')

th = 0;
pc_nr = 1;
length(find(u(hipp2pcc,pc_nr)>th))
length(find(u(hipp2pcc,pc_nr)<th))
length(find(u(amg2acc,pc_nr)>th))
length(find(u(amg2acc,pc_nr)<th))


%% Test MSD on average spectral thingy

xx_allsubs = squeeze(mean(S(this_f,:,:),1))';

% normalization does not change much here, since done later
% vector_length = sum(xx_allsubs.^ 2,2) .^ .5;
% xx_allsubs = xx_allsubs ./ repmat(vector_length,1,size(xx_allsubs,2)); % normalize (L2 norm) each trial

% xx_allsubs = abs(xx_allsubs);

% % construct distance matrix
% D_length = sqrt(sum(xx_allsubs.^2,2));
% D_dot = (xx_allsubs*xx_allsubs');
% D = D_dot./(D_length*D_length'); 
% D(D>1) = 1;
% D_mat = rad2deg(acos(D));
% [Y,eigvals] = cmdscale(round(D_mat),2);

% this is the same as the cosine distance:

% D_mat = pdist(xx_allsubs,'spearman'); % cosine or spearman
% [Y,eigvals] = cmdscale((D_mat),2);

% D_mat = pdist(xx_allsubs,'cosine'); % may like this, or cityblock?
% [Y,eigvals] = cmdscale((D_mat),2);

% figure,plot(eigvals)

D_mat = pdist(xx_allsubs,'correlation');
[Y] = mdscale(round(D_mat*100),6); % nonclassical

sub_symbol = {'>','o','<','d','*','x','+','^','s'};

figure('Position',[0 0 800 300]),hold on 

subplot(1,2,2),hold on
for kk = 1:8
    hipp2pcc = ismember(labels_allsubs(:,1),[0]) & all_sub_nrs==kk;
    h1 = scatter(Y(hipp2pcc,1),Y(hipp2pcc,2),80,[.2 .2 1],'Marker',sub_symbol{kk}); 
%     h1 = scatter(Y(hipp2pcc,1),Y(hipp2pcc,2),80,[.2 .2 1],'filled','MarkerFaceAlpha',.5,'Marker',sub_symbol{kk}); 
end
for kk = 1:8
    amg2acc = ismember(labels_allsubs(:,1),[1]) & all_sub_nrs==kk;
    h2 = scatter(Y(amg2acc,1),Y(amg2acc,2),80,[1 .2 .2],'Marker',sub_symbol{kk}); 
%     h2 = scatter(Y(amg2acc,1),Y(amg2acc,2),80,[1 .2 .2],'filled','MarkerFaceAlpha',.5,'Marker',sub_symbol{kk}); 
end

%% bad idea - too make colums
% % labels_allsubs = [zeros(size(out(3,1).subj_ind(out(3,1).p<0.05) )); ...
% %     1+zeros(size(out(4,2).subj_ind(out(4,2).p<0.05)))];
% Y = labels_allsubs;
% X = S_cat;
% 
% W = LDA(X,Y);
% %
% % % Calulcate linear scores for training data
% % L = [ones(25,1) X] * W';
% %
% % % Calculate class probabilities
% % P = exp(L) ./ repmat(sum(exp(L),2),[1 2]);


%% spectra one subject

figure
for s_nr = 1:8

    subplot(8,2,s_nr*2-1)
    imagesc(tt_v,f,mean(S(:,:,labels_allsubs==0 & all_sub_nrs==s_nr),3),[-8 8])
    axis xy
    
    ylabel('Frequency (Hz)')
    
    subplot(8,2,s_nr*2)
    imagesc(tt_v,f,mean(S(:,:,labels_allsubs==1 & all_sub_nrs==s_nr),3),[-8 8])
    axis xy
    
end
colormap(cm)
