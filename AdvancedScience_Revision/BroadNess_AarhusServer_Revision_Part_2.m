%% LBPD functions
% ALWAYS RUN THIS SECTION
%1) Add LBPD functions
% starting up some of the functions for the following analysis
pathl = '/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD'; %path to stored functions (THIS MUST BECOME A FOLDER ON YOUR OWN COMPUTER)
addpath(pathl);
LBPD_startup_D(pathl); % add to path LBPD and other OSL functions
addpath('/projects/MINDLAB2023_MEG-AuditMemDement/scripts/chiaramalvaso/Broadness')
%% CLUSTER CONFIGURATION
% ALWAYS RUN THIS SECTION IF YOU WANT TO USE CLUSTER
%clusterconfig('scheduler', 'cluster');
clusterconfig('scheduler', 'none'); %use this line to run local
clusterconfig('long_running',1);
clusterconfig('wait',1);
clusterconfig('slot',1); %choose the appropriate amount of memory slot

%% POLARITY CORRECTION ON MAIN EFFECT
% ALWAYS RUN THIS SECTION BEFORE ANY OTHER SECTION!!

% This section perform the polarity correction on the main effect

%loading data in 3559-voxel space (8mm)
%getting sign of the voxels based on the aggregated source reconstruction (over participants) - AVERAGED TRIALS
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/sources_main_effects.mat');

%actual computation
%adjusting polarity
timex = 45:52;
vect = zeros(3559,1);
for jj = 1:3559 %over brain voxels
    if squeeze(mean(t_val_s(jj,timex,1),2)) > 0 %if the data in voxel jj is positive during N100 time
        vect(jj,1) = -1; %storing a vector with 1 and -1 to be used for later statistics
    else
        vect(jj,1) = 1;
    end
end
dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
for cc = 1:size(dum,3) %over conditions
    for jj = 1:size(t_val_s,1) %over brain voxels
        dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
        disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
    end
end

%% PCA ON MAIN EFFECT (averaged subjects and conditions)
% ALWAYS RUN THIS SECTION BEFORE ANY OTHER SECTION!!

% 1) Time series using PCA on main effect. Saving a matrix with dimensions (time x component x condition x subjects)                                                                                        
% 2) Thresholded activation patterns for each component


%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');

%creating the output directory - choose the approriate name
dirname = 'AdvancedScience_Revision';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness';
if ~exist([outputdirpath '/' dirname], 'dir')
    mkdir([outputdirpath '/' dirname])
end
outputdir = [outputdirpath '/' dirname];
%%%%%%%%%%%%%%%%%%%%%%% PCA on main effect to get wcoeff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining the necessary input parameters
S = [];
S.H = dum(:,:,1:5); %MAIN EFFECT WITH POLARITY CORRECTION
S.permnum = 1;
S.fig_l = 0; % 1 to plot the results, 0 for no plotting
S.sign_eig = '0';
%S.namenii = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/Test'; %path were you store the .nii images
S.time = time;
S.rand_l = 1;
S.onefig = 0;

[ OUT ] = PCA_LBPD( S );

% extract the values of wcoeffs so you can use them to compute the time series for single subjects
wcoeff = OUT.W;

%%%%%%%%%%%%%%%%%%%%%%%%% Time series for each subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCs = 1:2; % array containing the principal components you want to consider

list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/SUBJ*.mat');
conds = 5; %number of experimental conditions
J = zeros(size(wcoeff,2),PCs(end),conds,length(list));

for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name])
    t_val_s = OUT.sources_ERFs;
    dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
    %adjusting polarity 
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(t_val_s,1) %over brain voxels
            dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            %            disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for iii = 1:conds %over experimental conditions
        J(:,:,iii,ii) = dum(:,1:size(wcoeff,2),iii)' * wcoeff(:,PCs); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
    end
    disp(ii)
end

%% RQA
% INPUT: J = 4-d matrix with the time series, dimensions are (time, components, conditions, subjects)
% 1) Phase space plot
% 2) Recurrence matrix (with and without thresholding)
% 3) Metrics on thresholded recurrence matrics and statistics
% 4) correlation between metrics and behavioural results and statistics

% This code stores the phase space coordinates in a cell array with
% dimensions (condition, participants). Each cell contains the phase space
% coordinate for each subject and experimental condition.


%creating the output directory - choose the approriate name
clear dirname
clear outputdirpath
clear outputdir
dirname = 'PhaseSpace_analysis';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/AdvancedScience_Revision';
if ~exist([outputdirpath '/' dirname], 'dir')
    mkdir([outputdirpath '/' dirname])
end
outputdir = [outputdirpath '/' dirname];

%%%%%%%%%%%%%%%%%%%%%%%%% PHASE SPACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PCs = 1:2; % array with the principal components you want to consider
timemin = 1; % starting time point
timemax = 651; % ending time point

phase_space = cell(size(J,3), size(J,4));
for cond = 1:size(J,3) % over conditions
    for sub = 1:size(J,4) % over subjects
    phase_space_temp = zeros(timemax-timemin+1,length(PCs)); % initializating the matrix that will contain the coordinates in the phase space
        for t = timemin:timemax % over time points
            for cc = 1:length(PCs) % over components
                phase_space_temp(t,cc) = J(t,PCs(cc),cond,sub);
            end
        end
    phase_space{cond,sub} = phase_space_temp;
    clear phase_space_temp
        
    end
end

save([outputdir '/phase_space_plot'], 'phase_space')


%%%%%%%%%%%%%%%%%%%%%% RECURRENCE PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RP = cell(size(phase_space,1), size(phase_space,2));
RP_thresh = cell(size(phase_space,1), size(phase_space,2));
eps = 0.1; %parameter for thresholding (5%)

% Computing the norm between every pair of points in the phase space
for cc = 1:size(phase_space,1) %over conditions
   for sub = 1:size(phase_space,2) %over subjects
       disp(['Condition ' num2str(cc) '- subject ' num2str(sub)])
       RP_temp = zeros(size(phase_space{cc,sub},1));
       bamba = zeros(size(RP_temp,1));
        for ii = 1:size(phase_space{cc,sub},1) %over time-points
            for jj = 1:size(phase_space{cc,sub},1) %over time-points
                RP_temp(ii,jj) = norm(phase_space{cc,sub}(ii,:) - phase_space{cc,sub}(jj,:)); %distance between the n components for time-points ii and jj
            end
        end
   bamba(RP_temp<max(RP_temp(:))*eps) = 1; %recurrent values
   RP{cc,sub} = RP_temp;
   RP_thresh{cc,sub} = bamba;
   clear RP_temp
   clear bamba
   end
end

%%%% plotting the recurrence plot with no thresholding %%%%
for cond = 1:size(RP,1)
    RP_temp = zeros(size(RP{1,1},1), size(RP{1,1},2), size(RP,2));
    for sub=1:size(RP,2) % over subjects
        RP_temp(:,:,sub) = RP{cond,sub};
    end
    % Computing the z score
    RP_temp_av = mean(RP_temp(:,:,:),3);
    RP_temp_std = std(RP_temp(:,:,:),0,3);
    bamba = zeros(size(RP_temp_av,1));
    bamba(RP_temp_av<max(RP_temp_av(:))*eps) = 1; %recurrent values
    figure; 
    imagesc(time(timemin:timemax),time(timemin:timemax), RP_temp_av./RP_temp_std); xlabel('Time (s)'); ylabel('Time (s)'); set(gca,'YDir','normal');
    colorbar
    set(gcf,'Color','w')
%     caxis([0.8 2.8])
%     title(['Condition ' num2str(cond) ' no thresholding'])
     print(gcf, [outputdir '/RP_condition_' num2str(cond) '_no_thresh'], '-dpdf', '-r300')
    clear RP_temp
    clear RP_temp_av
    clear bamba
end



%%%%%%%%%%%%%%%%%%%% Metrics on thresholded RP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('/projects/MINDLAB2017_MEG-LearningBach/scripts/Mattia/crptool'); % path to stored functions

lmin = 2; % minimum consecutive elements on diagonal lines

vmin = 2; % minimum consecutive elements on vertical lines

metrics = zeros(size(RP_thresh,1),size(RP_thresh,2), 8);

for cc = 1:size(RP_thresh,1) % over conditions
    
    for sub = 1:size(RP_thresh,2) % over subjects
        disp(['Condition ' num2str(cc) '- subject ' num2str(sub)])
        %%%%%%%%%%% RECURRENCE RATE (RR) %%%%%%%%%
        metrics(cc,sub,1) = sum(RP_thresh{cc,sub}(:)) / numel(RP_thresh{cc,sub}); %number of recurrence point divided by total number of points
        %%%%%%%%%%% mean diagonal line (L) %%%%%%%%%%%%%%%%
        clear Ldiags
        [~ , Ldiags] = dl(RP_thresh{cc,sub}); %Ldiags contains length of the diagonals found in the plot (showing that something is recurrent (with a time-lag)
        Ldiags(Ldiags<lmin) = [];
        metrics(cc,sub,2) = mean(Ldiags); %mean diagonals length
        %%%%%%%%%%% Determinism (DET) %%%%%%%%%%%%%%%
        if isempty(Ldiags)
            Ldiags = 0;
            warning('No diagonal elements in the RP');
        end
        if sum(RP_thresh{cc,sub}(:)) > 0
            metrics(cc,sub,3) = sum(Ldiags) / (sum(RP_thresh{cc,sub}(:))); %all diagonal summed divided by total number of recurrences = proportion of elements arranged in diagonals compared to all elements
        else
            metrics(cc,sub,3) = NaN;
            warning('No recurrence points in the RP');
        end
        %%%%%%%%%%%% Entropy (ENTR) %%%%%%%%%%%%%%%
        histL = hist(Ldiags(:), [1:min(size(RP_thresh{cc,sub}))]); %distribution of the number of occurrences of the diagonal of different lengths
        %when histograms are similar, high entropy, otherwise low entropy (when a particular size is very clearly emerging)
        metrics(cc,sub,4) = entropy(histL(:));
        %%%%%%%%%%%% Trapping time (TT) %%%%%%%%%%%%%
        clear TTverts
        [~, TTverts] = tt(RP_thresh{cc,sub}); %TTverts contains legnth of the vertical lines in the plot (showing that the plot is trapped in a state)
        TTverts(TTverts < vmin) = [];
        metrics(cc,sub,5) = mean(TTverts); %mean vertical lines length
        %%%%%%%%%%%% Laminarity %%%%%%%%%%%%%%%
        if sum(TTverts)>0
            metrics(cc,sub,6) = sum(TTverts) / (sum(RP_thresh{cc,sub}(:)));
        else
          metrics(cc,sub,6) = NaN;
        end
         
        %%%%%%%%%%% Maximal duration of laminar state, i.e, vertical line (Vmax) %%%%%%%%%%% 
        metrics(cc,sub,7) = max(TTverts);
        %%%%%%%%%%% Divergence (DIV; inverse of maximum diagonal length ) %%%%%%%%%%%
        clear Lmax
        Lmax = max(Ldiags(1:end-1));
        metrics(cc,sub,8) = 1 / Lmax;
    end
end

save([outputdir '/metrics_from_RP'], 'metrics')

%%%%%%%%%%%%%%%%%%%%%%%%%% Statistics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for now I'm performing a t test to find differences between each N condition and the previously memorized one M
pvals = zeros(size(metrics,1)-1, size(metrics,3));
tvals = zeros(size(metrics,1)-1, size(metrics,3));
for mm = 1:size(metrics,3) %over metrics
    for cc = 2:size(metrics,1) % over conditions
        [~,p,~,stats] = ttest(squeeze(metrics(1,:,mm)), squeeze(metrics(cc,:,mm))); % contrasting condition 1 vs all the others
        tvals(cc-1,mm) = stats.tstat;
        pvals(cc-1,mm) = p;
        clear p
    end
end

fdr_siglevel = zeros(size(pvals,1),1);
sig_pvals =zeros(size(pvals,1),size(pvals,2));
for cc = 1:size(pvals,1)
    [siglevel,~,~] = fdr(pvals(cc,:));
    fdr_siglevel(cc) = siglevel;
    for mm = 1:size(pvals,2)
        if pvals(cc,mm) <= fdr_siglevel(cc)
            sig_pvals(cc,mm) = pvals(cc,mm);
        else
            sig_pvals(cc,mm) = 1;
        end
    end
    clear siglevel
    
end

save([outputdir '/significant_pvals_on_metrics'], 'sig_pvals')

%%%%%%%%%%%%%%%%%% Correlation with behavioural data %%%%%%%%%%%%%%%%%%%%%%
% ACCURACY
accuracy = [cell2mat(Block_3_t{2:end,2}),  cell2mat(Block_3_t{2:end,4}), cell2mat(Block_3_t{2:end,6}), cell2mat(Block_3_t{2:end,8}), cell2mat(Block_3_t{2:end,10})];
pvals_accuracy = zeros(8,1);
corr_accuracy = zeros(8,1);
accuracy_av = mean(accuracy,2);
metrics_av = squeeze(mean(metrics,1));
for mm = 1:size(metrics_av,2)
    disp(mm)
    [r,p] = corr(accuracy_av,metrics_av(:,mm), 'Type', 'Spearman');
    pvals_accuracy(mm) = p;
    corr_accuracy(mm) = r;
%     figure
%     scatter(accuracy_av,metrics_av(:,mm))
end
fdr_accuracy = fdr(pvals_accuracy);

%REACTION TIME
rt = cell2mat(Block_3(2:end,14:18));
rt_av = mean(rt,2);
pvals_rt = zeros(8,1);
corr_rt = zeros(8,1);
metrics_av = squeeze(mean(metrics,1));
for mm = 1:size(metrics_av,2)
    disp(mm)
    [r,p] = corr(rt_av,metrics_av(:,mm), 'Type', 'Spearman');
    pvals_rt(mm) = p;
    corr_rt(mm) = r;
%     figure
%     scatter(accuracy_av,metrics_av(:,mm))
end
fdr_rt = fdr(pvals_rt);


pvals_tot = [pvals_accuracy; pvals_rt];
fdr_tot = fdr(pvals_tot);

%% SPATIAL GRADIENT
% 1) Activation patterns
% 2) K means clustering with different number of clusters
% 3) Find the optimal cluster number
% 4) Nifti images for the optimal cluster number
%creating the output directory - choose the approriate name
dirname = 'SpatialGradient_analysis';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/AdvancedScience_Revision';
if ~exist([outputdirpath '/' dirname], 'dir')
    mkdir([outputdirpath '/' dirname])
end
outputdir = [outputdirpath '/' dirname];

%%%%%%%%%%%%%%%%%%%%%%%% Activation Patterns %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wcoeff = OUT.W;
PCs = 1:2;
dum_averaged = mean(dum(:,:,1:5),3);
C = cov(dum_averaged');
SO = zeros(size(wcoeff,1),2);
for ii = 1:PCs(end) %over significant PCs
    act_pattern(:,ii) = OUT.W(:,ii)'*C;
    threshold = mean(abs(wcoeff(:,ii))) + std(abs(wcoeff(:,ii)));
    
    for jj = 1:size(wcoeff,1)
        if abs(wcoeff(jj,ii)) > threshold
            SO(jj,ii) = wcoeff(jj,ii);
        else
            SO(jj,ii) = 0;
        end
    end
end
save([outputdir '/ActivationPattern_thresholded'], 'SO')

%%%%%%%%%%%%%%%%%%%%%%%%% K-means clustering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUM = zeros(20,1);
SO_zscore = zscore(SO);
sil = zeros(20,1);
idx_ks = zeros(3559,20);
for k = 1:20
    disp(k)
    [idx, C, sumD, D]= kmeans(SO_zscore,k,'Replicates', 100);
    idx_ks(:,k) = idx;
%     evalclusters(SO_zscore,'kmeans', 'silhouette', 'KList', 2:20);
    SUM(k) = sum(sumD);
end

save([outputdir '/ElbowPlot_SUM'], 'SUM')
save([outputdir '/Cluster_idx'], 'idx_ks')

%%%%%%%%%%%%%%%%% Finding the optimal number of clusters %%%%%%%%%%%%%%%%%%
% Running evalclusters 1000 times to find the optimanl number of clusters (8 in this case)
eva_optimal = zeros(1000,1);
for ii = 1:1000
    disp(ii)
    eva = evalclusters(SO_zscore,'kmeans', 'silhouette', 'KList', 2:20);
    eva_optimal(ii) = eva.OptimalK;
end

%%%%%%%%%%%%%%%%%% Nifti images and cluster information %%%%%%%%%%%%%%%%%%%
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/AdvancedScience_Revision/SpatialGradient_analysis/MNI152_8mm_coord_dyi.mat')
optimalK = 8; % optimal number of clusters for k means
clusters_idx = squeeze(idx_ks(:,optimalK));
maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
clustered_SO = zeros(size(SO,1), size(SO,2), optimalK);
for ii = 1:optimalK %over clusters
    clear cluster_mask
    cluster_mask = clusters_idx == ii;
    clustered_SO(:,:,ii) = SO.*cluster_mask;
    
end
headers = {'Voxel','X', 'Y', 'Z', 'PC1', 'PC2'};
filename = 'clusters.xlsx';
if exist(filename, 'file') == 2
    delete(filename);
end
for ii = 1:optimalK % over clusters
    clear cluster_activations
    clear cluster_count
    clear cluster_binary
    clear SO_temp
    clear MNI_temp
    clear excel_data
    clear tbl
%     SO_temp = squeeze(clustered_SO(:,:,ii));
    cluster_binary = zeros(3559,1);
    cluster_count = 1;
    for jj = 1:3559 % over brain voxels
            
        if clusters_idx(jj) == ii 
            cluster_count = cluster_count + 1;
            cluster_activations(cluster_count,:) = SO(jj,:);
            MNI_temp(cluster_count,:) = MNI8(jj,:);
            cluster_binary(jj) = 1;
        else
            cluster_binary(jj) = 0;
            
        end
    end
    %%%%% excel file %%%%%%
    sheetname = ['cluster' num2str(ii)];
    excel_data = [(1:size(cluster_activations,1))', MNI_temp, cluster_activations];
    tbl = array2table(excel_data, 'VariableNames', headers);
    writetable(tbl, [outputdir '/' filename], 'Sheet', sheetname, 'Range', 'A1');
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); % no time here so size of 4D = 1
    for jj = 1:size(cluster_binary,1) %over brain sources
        dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = cluster_binary(jj,:); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - component ' num2str(ii)])
    save_nii(nii,[outputdir '/binary_ActPatterns_OPTIMALK_cluster_' num2str(ii) '.nii']); %printing image

end

%% INDEPENDENT COMPONENT ANALYSIS
% the following sections refers to ICA analysis
%% ICA WITH DIFFERENT NUMBER OF COMPONENTS
% This section perform ICA multiple times, forcing each time the extraction
% of M_max components
M_max = 26;
M_min = 2;

load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');

%creating the output directory - choose the approriate name
dirname = 'ICA_analysis';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/AdvancedScience_Revision';
if ~exist([outputdirpath '/' dirname], 'dir')
    mkdir([outputdirpath '/' dirname])
end
outputdir = [outputdirpath '/' dirname];


data = mean(dum(:,:,:),3); % averaging over conditions

% Demean data (Cohen de-means them in video #89)
data_demeaned = data - mean(data,2); % de meaning over time

list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/SUBJ*.mat');
conds = 5; %number of experimental conditions

sorted_comp = cell(M_max -M_min +1,size(J,3),size(J,2)); 
sorted_corr = cell(M_max -M_min +1, size(J,3),size(J,2));

for m = M_min:M_max
%%%%%%%%%%% Performing ICA on averaged subjects and conditions %%%%%%%%%%%%
    disp(['Performing ICA with ' num2str(m) ' components'])
    % Run ICA
    iVecs = jader(data_demeaned,m); % JADE algorithm is deterministic (i.e., same results when re-running)
    % function provided by Cohen (LAN course); read carefully function help
    % rows contain eigenvectors, unlike PCA and GED
    icScores = zeros(size(iVecs,1),size(dum,2),conds,length(list));
    
%%%%%%%%%%%%%%%%%%% Time series for each subject %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for ii = 1:length(list) %over subjects
        disp(['Computing time series for subject ' num2str(ii)])
        load([list(ii).folder '/' list(ii).name])
        t_val_s = OUT.sources_ERFs;
        dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
        %adjusting polarity %TO BE DON WITH ICA?
        for cc = 1:size(dum,3) %over conditions
            for jj = 1:size(t_val_s,1) %over brain voxels
                dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
                %            disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
            end
        end
        for iii = 1:conds %over experimental conditions
            icScores(:,:,iii,ii) = iVecs*dum(:,:,iii); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
        end
    end
    icScores_av = squeeze(mean(icScores(:,1:end-1,:,:),4)); % making ica time series of the same size as PCA time series
    
%%%%%%%%%%%%%%%%%%%%% storing correlation information %%%%%%%%%%%%%%%%%%%%%
    J_av = squeeze(mean(J,4));
 
    for cond = 1:size(J_av,3) % over experimental conditions
        for c_pca = 1:size(J_av,2) % over PCA components
            pca_ts_z = squeeze((J_av(:,c_pca,cond) - mean(J_av(:,c_pca,cond),1))/std(J_av(:,c_pca,cond),1));
            corr_temp = zeros(size(icScores_av,1),1);
            for c_ica = 1:size(icScores_av,1) %over ICA componets
                ica_ts_z = squeeze((icScores_av(c_ica,:,cond) - mean(icScores_av(c_ica,:,cond)))/std(icScores_av(c_ica,:,cond)));            
                corr_temp(c_ica) = abs(corr(pca_ts_z, ica_ts_z'));
            end

           [sorted_c, sorted_idx] = sort(corr_temp, 'descend');
           sorted_comp{m-1,cond,c_pca} = sorted_idx;
           sorted_corr{m-1,cond,c_pca} = sorted_c;
           clear sorted_c
           clear sorted_idx       
        end
    end

end

save([outputdir '/correlations_coeffs_' num2str(M_max) 'maxcomponents'] , 'sorted_corr');
save([outputdir '/correlations_index_' num2str(M_max) 'maxcomponents'], 'sorted_comp');

%% ICA WITH A CHOSEN NUMBER OF COMPONENTS
% 1) Perform ICA forcing the algorithm to extract a certain number (m) of
% components. m correspond to the number of principal components extracted
% with PCA that explain more than 95% (or 99%) of the total variance.
% 2) Compute the time series
% 3)Check the correlation of the ICA time series with the one obtained with
% PCA


%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');

%creating the output directory - choose the approriate name
dirname = 'ICA_analysis';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/AdvancedScience_Revision';
if ~exist([outputdirpath '/' dirname], 'dir')
    mkdir([outputdirpath '/' dirname])
end
outputdir = [outputdirpath '/' dirname];

%%%%%%%%%%%%%%%%%%%%%%%%%%% actual ICA computation %%%%%%%%%%%%%%%%%%%%%%%%

data = mean(dum(:,:,:),3); % averaging over conditions

% Demean data (Cohen de-means them in video #89)
data_demeaned = data - mean(data,2); % de meaning over time

% Finding the appropriate number of components to extract using PCA
total_varexp = 0.95;
if isempty(total_varexp)
    m = 2; % forcing the alghoritm to extract 2 components
else
    [~,~,latent] = pca(data_demeaned');
    explained = cumsum(latent) / sum (latent);
    m = find(explained >= total_varexp, 1);
end
% Run ICA
iVecs = jader(data_demeaned,m); % JADE algorithm is deterministic (i.e., same results when re-running)
% function provided by Cohen (LAN course); read carefully function help
% rows contain eigenvectors, unlike PCA and GED

% compute IC scores (activation timeseries)
% reconstructing the time series independently for each condition
%%%%%%%%%%%%%%%%%%%%%%%%% Time series for each subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/SUBJ*.mat');
conds = 5; %number of experimental conditions
icScores = zeros(size(iVecs,1),size(dum,2),conds,length(list));

for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name])
    t_val_s = OUT.sources_ERFs;
    dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
    %adjusting polarity %TO BE DON WITH ICA?
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(t_val_s,1) %over brain voxels
            dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            %            disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for iii = 1:conds %over experimental conditions
        icScores(:,:,iii,ii) = iVecs*dum(:,:,iii); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
    end
    disp(ii)
end

save([outputdir '/all_ICA_components_' num2str(m)], 'icScores')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjust the size of the data and select time points corresponding to the time interval that's relevant for the statistic
starting_time = 0.350; %starting time IN SECONDS
ending_time = 2.50;    %final time IN SECONDS
reduced_index = find(time >= starting_time & time <= ending_time);


%dataa = permute(J, [2,1,3,4]);
dataa = icScores; % time series coming from ICA have already the proper dimensions
dataa_lesstime = dataa(:,reduced_index,:,:);
lesstime = time(reduced_index);

comparisons = 4;
PCs = size(icScores,1); % number of components I want to see

%check the values to see if they make sense
if comparisons > (size(dataa_lesstime,3)-1)
    error(['ERROR: you want to see ' num2str(comparisons) ' components but you only have ' num2str(size(dataa_lesstime,3)) ' conditions'])
end

if PCs > size(dataa_lesstime,1)
    error(['ERROR: you want to see ' num2str(PCs) ' components but you only computed ' num2str(size(dataa_lesstime,1)) ' components'])
end


% J = (time-points, principal components,conditions,subjects)
P = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T-test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for ii = 1:size(dataa_lesstime,1) %over principal components
    for jj = 1:size(dataa_lesstime,2) %ove time-points
        for cc = 1:(size(dataa_lesstime,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa_lesstime(ii,jj,1,:)),squeeze(dataa_lesstime(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear FDR
fdr_siglevel = zeros(PCs,  comparisons);
for comp = 1:PCs %iterate over principal components
    for cond = 1:comparisons % iterate over conditions
        [fdr_siglevel(comp,cond),~,~] = fdr(P(comp,:,cond));
    end
end

fdr_results = cell(PCs,comparisons);

for pp = 1:PCs
    for cc = 1:comparisons
         fdr_matrix = cell(4,length(lesstime)); %initializing the matrix that will be saved on excel
         fdr_matrix(1,:) = num2cell(lesstime);
         clear cluster;
         clear mask;
         mask =double(P(pp,:,cc) < fdr_siglevel(pp,cc));
         cluster = bwconncomp(mask);
         fdr_matrix(2,:) = num2cell(mask);
         for tt = 1:length(lesstime)
             if mask(tt) == 1
                 fdr_matrix{3,tt} = P(pp,tt,cc);
                 fdr_matrix{4,tt} = T(pp,tt,cc);
             else
                 fdr_matrix{3,tt} = 'N.S.';
                 fdr_matrix{4,tt} = 'N.S.';
             end
         end

             filename = [outputdir '/FDR_allICA_network' num2str(pp) '_comparison_1VS' num2str(cc+1) '.xlsx'];
             writetable(cell2table(fdr_matrix),filename);
         S = {};
         temp_sigwindows = cell(cluster.NumObjects,1);
         for n = 1:cluster.NumObjects
              idx = cluster.PixelIdxList{n};
             S(n).size = length(idx);
             S(n).time_interval = [lesstime(idx(1)) lesstime(idx(end))]; 
             S(n).min_pval = min(P(pp,idx,cc));
             [absmax, indexmax] = max(abs(T(pp,idx,cc)));
             S(n).Tmax = T(pp,idx(1)-1+indexmax,cc); %storing the Tvalue with the sign
             temp_sigwindows{n,1} = [lesstime(idx(1)) lesstime(idx(end))];
         end
          fdr_results{pp,cc} = temp_sigwindows;
    end
end


save([outputdir '/fdr_results_allICA_' num2str(m) 'components'], 'fdr_results');

%% SAVING DATA ONLY FOR BEST ICA COMPONENTS
% This section saves data for the main plot, so the one with the components
% that correlate more with the principal components from PCA. It's the same
% as the previous section, but in this case saves the outputs only for the
% ICA componenents that correlate with the PCA components

%creating the output directory - choose the approriate name
dirname = 'ICA_analysis';
outputdirpath = '/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/Broadness/AdvancedScience_Revision';
if ~exist([outputdirpath '/' dirname], 'dir')
    mkdir([outputdirpath '/' dirname])
end
outputdir = [outputdirpath '/' dirname];

%%%%%%%%%%%%%%%%%%%%%%% actual ICA computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loading time
load('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/time_normal.mat');

data = mean(dum(:,:,:),3); % averaging over conditions

% Demean data (Cohen de-means them in video #89)
data_demeaned = data - mean(data,2); % de meaning over time

% Finding the appropriate number of components to extract using PCA
total_varexp = [];
if isempty(total_varexp)
    m = 2; % forcing the alghoritm to extract 2 components
else
    [~,~,latent] = pca(data_demeaned');
    explained = cumsum(latent) / sum (latent);
    m = find(explained >= total_varexp, 1);
end
% Run ICA
iVecs = jader(data_demeaned,m); % JADE algorithm is deterministic (i.e., same results when re-running)
% function provided by Cohen (LAN course); read carefully function help
% rows contain eigenvectors, unlike PCA and GED

% compute IC scores (activation timeseries)
% reconstructing the time series independently for each condition
%%%%%%%%%%%%%%%%%%%%%%%%% Time series for each subject %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


list = dir('/scratch7/MINDLAB2023_MEG-AuditMemDement/chiaramalvaso/single_trial_0/SUBJ*.mat');
conds = 5; %number of experimental conditions
icScores = zeros(size(iVecs,1),size(dum,2),conds,length(list));

for ii = 1:length(list) %over subjects
    load([list(ii).folder '/' list(ii).name])
    t_val_s = OUT.sources_ERFs;
    dum = zeros(size(t_val_s,1),size(t_val_s,2),size(t_val_s,3));
    %adjusting polarity %TO BE DON WITH ICA?
    for cc = 1:size(dum,3) %over conditions
        for jj = 1:size(t_val_s,1) %over brain voxels
            dum(jj,:,cc) = t_val_s(jj,:,cc) .* vect(jj,1); %reversing (or not)..
            %            disp(['condition ' num2str(cc) ' - source ' num2str(jj)])
        end
    end
    for iii = 1:conds %over experimental conditions
        icScores(:,:,iii,ii) = iVecs*dum(:,:,iii); %matrix multiplication for getting a timeseries obtained by multiplying, for each time-point, each voxel activation by its corresponding load
    end
    disp(ii)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Thresholded brain images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotting weights in the brain (nifti images)
dum_averaged = mean(dum(:,:,1:5),3);
C = cov(dum_averaged');

maskk = load_nii('/projects/MINDLAB2017_MEG-LearningBach/scripts/Leonardo_FunctionsPhD/External/MNI152_8mm_brain_diy.nii.gz'); %getting the mask for creating the figure
for ii = 1:size(iVecs,1) %over significant PCs
    wcoeff(:,ii) = iVecs(ii,:)*C;
    threshold = mean(abs(wcoeff(:,ii))) + std(abs(wcoeff(:,ii)));
    SO = zeros(size(wcoeff,1),1);
    for jj = 1:size(wcoeff,1)
        if abs(wcoeff(jj,ii)) > threshold
            SO(jj) = wcoeff(jj,ii);
        else
            SO(jj) = 0;
        end
    end
    %building nifti image
    SS = size(maskk.img);
    dumimg = zeros(SS(1),SS(2),SS(3),1); %no time here so size of 4D = 1
    for jj = 1:size(SO,1) %over brain sources
        dumm = find(maskk.img == jj); %finding index of sources ii in mask image (MNI152_8mm_brain_diy.nii.gz)
        [i1,i2,i3] = ind2sub([SS(1),SS(2),SS(3)],dumm); %getting subscript in 3D from index
        dumimg(i1,i2,i3,:) = SO(jj,:); %storing values for all time-points in the image matrix
    end
    nii = make_nii(dumimg,[8 8 8]); %making nifti image from 3D data matrix (and specifying the 8 mm of resolution)
    nii.img = dumimg; %storing matrix within image structure
    nii.hdr.hist = maskk.hdr.hist; %copying some information from maskk
    disp(['saving nifti image - component ' num2str(ii)])
    save_nii(nii,[outputdir '/thresholded_actpattern_allICA_' num2str(m) '_network'  num2str(ii) '.nii']); %printing image
end



%%%%%%%%%%%% ordering the time series based on correlation %%%%%%%%%%%%%%%%%%%
icScores_av = squeeze(mean(icScores(:,1:end-1,:,:),4)); % making ica time series of the same size as PCA time series
J_av = squeeze(mean(J,4));
sorted_comp = cell(size(J_av,3),size(J_av,2)); 
sorted_corr = cell(size(J_av,3),size(J_av,2)); 
for cond = 1:size(J_av,3) % over experimental conditions
    for c_pca = 1:size(J_av,2) % over PCA components
        pca_ts_z = squeeze((J_av(:,c_pca,cond) - mean(J_av(:,c_pca,cond),1))/std(J_av(:,c_pca,cond),1));
        corr_temp = zeros(size(icScores_av,1),1);
        for c_ica = 1:size(icScores_av,1) %over ICA componets
            ica_ts_z = squeeze((icScores_av(c_ica,:,cond) - mean(icScores_av(c_ica,:,cond)))/std(icScores_av(c_ica,:,cond)));            
            corr_temp(c_ica) = abs(corr(pca_ts_z, ica_ts_z'));
        end
        
       [sorted_c, sorted_idx] = sort(corr_temp, 'descend');
       sorted_comp{cond,c_pca} = sorted_idx;
       sorted_corr{cond,c_pca} = sorted_c;
       clear sorted_c
       clear sorted_idx       
    end
end

best_ICA = zeros(2, size(icScores,2),size(icScores,3),size(icScores,4));
for cc = 1:size(best_ICA,1)
    for conds = 1:5
        best_ICA(cc,:,conds,:) = icScores(sorted_comp{conds,cc}(1),:,conds,:);
    end
end



save([outputdir '/best_ICA_components_' num2str(m)], 'best_ICA')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Statistic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%adjust the size of the data and select time points corresponding to the time interval that's relevant for the statistic
starting_time = 0.350; %starting time IN SECONDS
ending_time = 2.50;    %final time IN SECONDS
reduced_index = find(time >= starting_time & time <= ending_time);


%dataa = permute(J, [2,1,3,4]);
dataa = best_ICA; % time series coming from ICA have already the proper dimensions
dataa_lesstime = dataa(:,reduced_index,:,:);
lesstime = time(reduced_index);

comparisons = 4;
PCs = size(best_ICA,1); % number of components I want to see

%check the values to see if they make sense
if comparisons > (size(dataa_lesstime,3)-1)
    error(['ERROR: you want to see ' num2str(comparisons) ' components but you only have ' num2str(size(dataa_lesstime,3)) ' conditions'])
end

if PCs > size(dataa_lesstime,1)
    error(['ERROR: you want to see ' num2str(PCs) ' components but you only computed ' num2str(size(dataa_lesstime,1)) ' components'])
end


% J = (time-points, principal components,conditions,subjects)
P = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)
T = zeros(size(dataa_lesstime,1),size(dataa_lesstime,2),(size(dataa_lesstime,3))-1); %PCs x time-points x contrasts (every NewTX versus Old)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T-test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:size(dataa_lesstime,1) %over principal components
    for jj = 1:size(dataa_lesstime,2) %over time-points
        for cc = 1:(size(dataa_lesstime,3)-1) %over the 4 NewTX
            %codes t-test
            [~,p,~,stats] = ttest(squeeze(dataa_lesstime(ii,jj,1,:)),squeeze(dataa_lesstime(ii,jj,cc+1,:))); %contrasting cond1 (Old) with cond X (depending on vectcond it can be either NewT1 or NewT2 or NewT3 or NewT4
            P(ii,jj,cc) = p;
            T(ii,jj,cc) = stats.tstat;
        end
        disp([num2str(ii) ' - ' num2str(jj)])
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FDR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear FDR
fdr_siglevel = zeros(PCs,  comparisons);
for comp = 1:PCs %iterate over principal components
    for cond = 1:comparisons % iterate over conditions
        [fdr_siglevel(comp,cond),~,~] = fdr(P(comp,:,cond));
    end
end

fdr_results = cell(PCs,comparisons);

for pp = 1:PCs
    for cc = 1:comparisons
         fdr_matrix = cell(4,length(lesstime)); %initializing the matrix that will be saved on excel
         fdr_matrix(1,:) = num2cell(lesstime);
         clear cluster;
         clear mask;
         mask =double(P(pp,:,cc) < fdr_siglevel(pp,cc));
         cluster = bwconncomp(mask);
         fdr_matrix(2,:) = num2cell(mask);
         for tt = 1:length(lesstime)
             if mask(tt) == 1
                 fdr_matrix{3,tt} = P(pp,tt,cc);
                 fdr_matrix{4,tt} = T(pp,tt,cc);
             else
                 fdr_matrix{3,tt} = 'N.S.';
                 fdr_matrix{4,tt} = 'N.S.';
             end
         end

             filename = [outputdir '/FDR_bestICA_components' num2str(m) '_network' num2str(pp) '_comparison_1VS' num2str(cc+1) '.xlsx'];
             writetable(cell2table(fdr_matrix),filename);
         S = {};
         temp_sigwindows = cell(cluster.NumObjects,1);
         for n = 1:cluster.NumObjects
              idx = cluster.PixelIdxList{n};
             S(n).size = length(idx);
             S(n).time_interval = [lesstime(idx(1)) lesstime(idx(end))]; 
             S(n).min_pval = min(P(pp,idx,cc));
             [absmax, indexmax] = max(abs(T(pp,idx,cc)));
             S(n).Tmax = T(pp,idx(1)-1+indexmax,cc); %storing the Tvalue with the sign
             temp_sigwindows{n,1} = [lesstime(idx(1)) lesstime(idx(end))];
         end
          fdr_results{pp,cc} = temp_sigwindows;
    end
end


save([outputdir '/fdr_results_bestICA_' num2str(m) 'components'], 'fdr_results');



%%