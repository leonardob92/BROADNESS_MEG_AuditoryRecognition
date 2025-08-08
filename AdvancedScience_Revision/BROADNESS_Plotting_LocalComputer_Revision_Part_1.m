
%% LOADING DATA AND ADDING PATH

addpath('/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code')

load('/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/Data.mat')

%%

%% PLOTTING FOR ADVANCED SCIENCE REVISION - BEHAVIORAL MEASURES CORRELATED WITH THE TIME SERIES OF THE BRAIN NETWORKS (no RTs which are instead simply reported in the supplementary tables and in the main manuscript

%% 1) MEMORISED - brain network 1 - accuracy

sigIdx1 = SIGN{1,1,1,2}; %brain network 1 - Old - accuracy
% Apply to each cell in the input
win1 = cluster_sig_indices(sigIdx1, time);
SignificantWindows = {
    win1;
};

ConditionXGroupColors = {
    [0 0 0];
};
% Use the function
BROADNESS_Plot_ActivationTimeseries(data(1,:,:,1), time, ...
    'ConditionLabels', {'Mem', 'NT1', 'NT2', 'NT3', 'NT4'}, ...
    'ConditionXGroupColors', ConditionXGroupColors, ...
    'SignificantWindows', SignificantWindows, ...
    'SignificanceStyle', 'line', ...
    'SignificanceColors', [1 0 0], ...
    'STEStyle', 2, ...
    'Transparency', 0.3, ...
    'SignificanceLineOffset', 0.01, ...
    'SignificanceLinePadding', 0.05, ...
    'SignificanceLinePosition', 'below', ...
    'XLimits', [-0.1 3]);

exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/BN1_Mem_accuracy.pdf'],'Resolution',300)

%% 2) MEMORISED - brain network 2 - accuracy

sigIdx1 = SIGN{2,1,1,2}; %brain network 1 - Old - accuracy
% Apply to each cell in the input
win1 = cluster_sig_indices(sigIdx1, time);
SignificantWindows = {
    win1;
};

ConditionXGroupColors = {
    [0 0 0];
};
% Use the function
BROADNESS_Plot_ActivationTimeseries(data(2,:,:,1), time, ...
    'ConditionLabels', {'Mem', 'NT1', 'NT2', 'NT3', 'NT4'}, ...
    'ConditionXGroupColors', ConditionXGroupColors, ...
    'SignificantWindows', SignificantWindows, ...
    'SignificanceStyle', 'line', ...
    'SignificanceColors', [0 0 1], ...
    'STEStyle', 2, ...
    'Transparency', 0.3, ...
    'SignificanceLineOffset', 0.03, ...
    'SignificanceLinePadding', 0.05, ...
    'SignificanceLinePosition', 'below', ...
    'XLimits', [-0.1 3]);

exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/BN2_Mem_accuracy.pdf'],'Resolution',300)

%% 3) NT2 - brain network 2 - musical expertise

sigIdx1 = SIGN{2,3,3,2};
% Apply to each cell in the input
win1 = cluster_sig_indices(sigIdx1, time);
SignificantWindows = {
    win1;
};

ConditionXGroupColors = {
    [0 0 0];
};
% Use the function
BROADNESS_Plot_ActivationTimeseries(data(2,:,:,3), time, ...
    'ConditionLabels', {'NT1', 'NT2', 'NT3', 'NT4'}, ...
    'ConditionXGroupColors', ConditionXGroupColors, ...
    'SignificantWindows', SignificantWindows, ...
    'SignificanceStyle', 'line', ...
    'SignificanceColors', [0 0 1; 1 0 0], ...
    'STEStyle', 2, ...
    'Transparency', 0.3, ...
    'SignificanceLineOffset', 0.03, ...
    'SignificanceLinePadding', 0.05, ...
    'SignificanceLinePosition', 'below', ...
    'XLimits', [-0.1 3]);

exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/BN2_NT2_music.pdf'],'Resolution',300)


%% 4) NT3 - brain network 2 - musical expertise

sigIdx1 = SIGN{2,4,3,2};
% Apply to each cell in the input
win2 = cluster_sig_indices(sigIdx1, time);

SignificantWindows = {
    win2;
};

ConditionXGroupColors = {
    [0 0 0];
};
% Use the function
BROADNESS_Plot_ActivationTimeseries(data(2,:,:,4), time, ...
    'ConditionLabels', {'NT3', 'NT4'}, ...
    'ConditionXGroupColors', ConditionXGroupColors, ...
    'SignificantWindows', SignificantWindows, ...
    'SignificanceStyle', 'line', ...
    'SignificanceColors', [0 0 1], ...
    'STEStyle', 2, ...
    'Transparency', 0.3, ...
    'SignificanceLineOffset', 0.03, ...
    'SignificanceLinePadding', 0.05, ...
    'SignificanceLinePosition', 'below', ...
    'XLimits', [-0.1 3]);

exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/BN2_NT3_music.pdf'],'Resolution',300)

%%

%% HISTOGRAM FOR PERMUTATIONS

% Define bins centered on 1, 2, 3, 4
edges = 0.5:1:4.5;
% Plot histogram
figure;
histogram(T_cluster, edges,'FaceColor', [0 0 0.3]);
% Set white background
set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
% Add axis labels
xlabel('Cluster size');
ylabel('Number of occurrencies');
% Set x-axis ticks to integer values
xticks(1:4);
% Optional: Add title
title('Discrete Histogram');

exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/Hist.pdf'],'Resolution',300)

%%

%% RANDOMIZATIONS - ADDITIONAL METRICS

%% loding data

load('/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/Data2.mat');

%% actual computation of randomizations metrics

%%% USER INPUTS %%%
perm_n = 100; %number of permutations
perm_l = 2; %1 = randomization over time; 2 = randomizationover space (only randomizing voxels without disrupting the temporal order of their time series
%%%

% computing permutations
data = zeros(size(averaged_data,1),size(averaged_data,2),perm_n + 1);
averaged_data = mean(averaged_data,3); %average over experimental conditions
data(:,:,1) = averaged_data; %original data
if perm_l == 1
    for permi = 1:perm_n
        data_reshaped = zeros(size(averaged_data,1),size(averaged_data,2)); %preallocating matrix for randomized data
        for sourci = 1:size(averaged_data,1) %over brain sources
            idx_dummy = randperm(size(averaged_data,2)); %creating a permuted array of the indices of the temporal dimension
            data_reshaped(sourci,:) = averaged_data(sourci,idx_dummy); %actual reshaped data
        end
        data(:,:,permi + 1) = data_reshaped;
        disp(permi)
    end
else
    for permi = 1:perm_n
        idx_dummy = randperm(size(averaged_data,1)); %creating a permuted array of the indices of the spatial dimension
        data_reshaped = averaged_data(idx_dummy,:); %actual reshaped data
        data(:,:,permi + 1) = data_reshaped;
        disp(permi)
    end
end
data = permute(data,[2 1 3]); %permuting since I need data to be time x sources x datasets (original data + randomizations)

%extracting information
M = size(data, 3);
N = size(data, 2);

% Step 1: Compute FC matrices
FC = zeros(N, N, M);
for i = 1:M
    FC(:, :, i) = corr(squeeze(data(:, :, i)));  % Pearson correlation
end

% Step 2: zeroing diagonal 
% Get size
[N2, ~, ~] = size(FC);
% Create mask
diagIdx = repmat(eye(N2) == 1, 1, 1, M);
% Set diagonals to 0
FC(diagIdx) = 0;

% Step 3: Fisher z-transform
zFC = atanh(FC);  % Fisher's z

% Plot settings
% clim = [-1 1];  % Set color limits for comparability (adjust as needed)
% --- Plot original z-transformed FC ---
figure;
set(gcf, 'Color', 'w');
imagesc(zFC(:,:,1));  % original
colorbar;
title('Original z-FC Matrix', 'FontSize', 14);
xlabel('Region', 'FontSize', 12);
ylabel('Region', 'FontSize', 12);
box on;
axis square;
set(gca, 'FontSize', 10);
exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/Figure_Perm/Perm_' num2str(perm_l) '_Extra_4.pdf'],'Resolution',300)

% --- Plot first permutation z-transformed FC ---
figure;
set(gcf, 'Color', 'w');
imagesc(zFC(:,:,2));  % first permutation
colorbar;
title('Permutation 1 z-FC Matrix', 'FontSize', 14);
xlabel('Region', 'FontSize', 12);
ylabel('Region', 'FontSize', 12);
box on;
axis square;
set(gca, 'FontSize', 10);
exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/Figure_Perm/Perm_' num2str(perm_l) '_Extra_5.pdf'],'Resolution',300)

% Step 4: Compute similarity to original AND all pairwise similarities
% Define upper triangle mask (excluding diagonal)
N = size(zFC, 1);
M = size(zFC, 3);
mask = triu(true(N), 1);  % mask for upper triangle only
% Preallocate matrix: each row = vectorized FC matrix (upper triangle only)
zFC_vecs = zeros(M, nnz(mask));  % M FCs, each as a 1D vector of unique edges
% Vectorize all zFC matrices
for i = 1:M
    temp = zFC(:,:,i);
    temp = temp(mask);  % extract upper triangle
    zFC_vecs(i, :) = temp;
end
% Compute similarity to original FC (row 1)
similarity_to_original = zeros(M-1, 1);
for i = 2:M
    validIdx = isfinite(zFC_vecs(1,:)) & isfinite(zFC_vecs(i,:));
    similarity_to_original(i-1) = corr(zFC_vecs(1,validIdx)', zFC_vecs(i,validIdx)');
end
% Compute full pairwise similarity matrix
% Rows = FCs, Columns = FCs
similarity_matrix = corr(zFC_vecs', 'rows', 'pairwise');  % handles NaNs
similarity_matrix(logical(eye(M))) = 0;

% --- Plot histogram ---
figure;
set(gcf, 'Color', 'w');
histogram(similarity_to_original, 20, 'FaceColor', [0.3 0.6 0.8], 'EdgeColor', 'k');
xlabel('Correlation with Original FC', 'FontSize', 12);
ylabel('Count', 'FontSize', 12);
title('Similarity of Randomized FCs to Original', 'FontSize', 14);
box on;
grid on;
exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/Figure_Perm/Perm_' num2str(perm_l) '_Extra_1.pdf'],'Resolution',300)

% --- Plot similarity matrix ---
figure;
set(gcf, 'Color', 'w');
imagesc(similarity_matrix);
colorbar;
title('Pairwise FC Similarities (Fisher z)', 'FontSize', 14);
xlabel('FC Index', 'FontSize', 12);
ylabel('FC Index', 'FontSize', 12);
box on;
exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/Figure_Perm/Perm_' num2str(perm_l) '_Extra_2.pdf'],'Resolution',300)

%%

% Step 5: PCA on FCs
zFC_vecs_square = zeros((N*(N-1))/2, M);
for i = 1:M
    zFC_vecs_square(:, i) = squareform(zFC(:, :, i), 'tovector');
end
zFC_vecs_centered = zFC_vecs_square - mean(zFC_vecs_square, 1);
[coeff, score, latent] = pca(zFC_vecs_centered');

% --- PCA scatter plot ---
figure;
set(gcf, 'Color', 'w');
scatter(score(:,1), score(:,2), 60, 1:size(score,1), 'filled');
xlabel('PC1', 'FontSize', 12);
ylabel('PC2', 'FontSize', 12);
title('Cov-STATIS PCA of FC Matrices', 'FontSize', 14);
colormap(parula(M));
colorbar;
box on;
grid on;
exportgraphics(gcf,['/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/Figure_Perm/Perm_' num2str(perm_l) '_Extra_3.pdf'],'Resolution',300)

%%

%%

%%

%% EXTRACTING INFORMATION ABOUT THE RELATIONSHIP BETWEEN BEHAVIOURAL MEASURES AND BRAIN NETWORKS TIME SERIES AND PRINTING IT ON EXCEL FILES

%% 1) MEMORISED - brain network 1 - accuracy

sigIdx1 = SIGN{1,1,1,2}; %brain network 1 - Old - accuracy
% Apply to each cell in the input
win1 = cluster_sig_indices(sigIdx1, time);
SignificantWindows = {
    win1;
};

ConditionXGroupColors = {
    [0 0 0];
};
% Use the function
BROADNESS_Plot_ActivationTimeseries(data(1,:,:,1), time, ...
    'ConditionLabels', {'Mem', 'NT1', 'NT2', 'NT3', 'NT4'}, ...
    'ConditionXGroupColors', ConditionXGroupColors, ...
    'SignificantWindows', SignificantWindows, ...
    'SignificanceStyle', 'line', ...
    'SignificanceColors', [1 0 0], ...
    'STEStyle', 2, ...
    'Transparency', 0.3, ...
    'SignificanceLineOffset', 0.01, ...
    'SignificanceLinePadding', 0.05, ...
    'SignificanceLinePosition', 'below', ...
    'XLimits', [-0.1 3]);

sigIdx1 = SIGN{1,1,1,3}; %brain network 1 - Old - accuracy
writematrix(sigIdx1, '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/1.xlsx', 'Sheet', 1, 'Range', 'A1');


%% 2) MEMORISED - brain network 2 - accuracy AND RTS

sigIdx1 = SIGN{2,1,1,2}; %brain network 1 - Old - accuracy
% Apply to each cell in the input
win1 = cluster_sig_indices(sigIdx1, time);
sigIdx2 = SIGN{2,1,2,2}; %brain network 1 - Old - RTs
% Apply to each cell in the input
win2 = cluster_sig_indices(sigIdx2, time);

SignificantWindows = {
    win1; win2;
};

ConditionXGroupColors = {
    [0 0 0];
};
% Use the function
BROADNESS_Plot_ActivationTimeseries(data(2,:,:,1), time, ...
    'ConditionLabels', {'Mem', 'NT1', 'NT2', 'NT3', 'NT4'}, ...
    'ConditionXGroupColors', ConditionXGroupColors, ...
    'SignificantWindows', SignificantWindows, ...
    'SignificanceStyle', 'line', ...
    'SignificanceColors', [0 0 1; 1 0 0], ...
    'STEStyle', 2, ...
    'Transparency', 0.3, ...
    'SignificanceLineOffset', 0.03, ...
    'SignificanceLinePadding', 0.05, ...
    'SignificanceLinePosition', 'below', ...
    'XLimits', [-0.1 3]);

sigIdx1 = SIGN{2,1,1,3}; %brain network 1 - Old - accuracy
writematrix(sigIdx1, '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/2.xlsx', 'Sheet', 1, 'Range', 'A1');
sigIdx2 = SIGN{2,1,2,3}; %brain network 1 - Old - accuracy
writematrix(sigIdx2, '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/3.xlsx', 'Sheet', 1, 'Range', 'A1');


%% 3) NT1 - brain network 1 - RTS

sigIdx1 = SIGN{1,2,2,2};
% Apply to each cell in the input
win1 = cluster_sig_indices(sigIdx1, time);
SignificantWindows = {
    win1;
};

ConditionXGroupColors = {
    [0 0 0];
};
% Use the function
BROADNESS_Plot_ActivationTimeseries(data(1,:,:,2), time, ...
    'ConditionLabels', {'NT1', 'NT2', 'NT3', 'NT4'}, ...
    'ConditionXGroupColors', ConditionXGroupColors, ...
    'SignificantWindows', SignificantWindows, ...
    'SignificanceStyle', 'line', ...
    'SignificanceColors', [0 0 1], ...
    'STEStyle', 2, ...
    'Transparency', 0.3, ...
    'SignificanceLineOffset', 0.03, ...
    'SignificanceLinePadding', 0.05, ...
    'SignificanceLinePosition', 'below', ...
    'XLimits', [-0.1 3]);

sigIdx1 = SIGN{1,2,2,3}; %brain network 1 - Old - accuracy
writematrix(sigIdx1, '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/4.xlsx', 'Sheet', 1, 'Range', 'A1');

%% 4) NT2 - brain network 2 - musical expertise

sigIdx1 = SIGN{2,3,3,2};
% Apply to each cell in the input
win1 = cluster_sig_indices(sigIdx1, time);
SignificantWindows = {
    win1;
};

ConditionXGroupColors = {
    [0 0 0];
};
% Use the function
BROADNESS_Plot_ActivationTimeseries(data(2,:,:,3), time, ...
    'ConditionLabels', {'NT1', 'NT2', 'NT3', 'NT4'}, ...
    'ConditionXGroupColors', ConditionXGroupColors, ...
    'SignificantWindows', SignificantWindows, ...
    'SignificanceStyle', 'line', ...
    'SignificanceColors', [0 0 1; 1 0 0], ...
    'STEStyle', 2, ...
    'Transparency', 0.3, ...
    'SignificanceLineOffset', 0.03, ...
    'SignificanceLinePadding', 0.05, ...
    'SignificanceLinePosition', 'below', ...
    'XLimits', [-0.1 3]);

sigIdx1 = SIGN{2,3,3,3}; %brain network 1 - Old - accuracy
writematrix(sigIdx1, '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/6.xlsx', 'Sheet', 1, 'Range', 'A1');


%% 5) NT3 - brain network 2 - RTs and musical expertise

sigIdx1 = SIGN{2,4,2,2};
% Apply to each cell in the input
win1 = cluster_sig_indices(sigIdx1, time);
sigIdx1 = SIGN{2,4,3,2};
% Apply to each cell in the input
win2 = cluster_sig_indices(sigIdx1, time);

SignificantWindows = {
    win1;win2;
};

ConditionXGroupColors = {
    [0 0 0];
};
% Use the function
BROADNESS_Plot_ActivationTimeseries(data(2,:,:,4), time, ...
    'ConditionLabels', {'NT3', 'NT4'}, ...
    'ConditionXGroupColors', ConditionXGroupColors, ...
    'SignificantWindows', SignificantWindows, ...
    'SignificanceStyle', 'line', ...
    'SignificanceColors', [1 0 0; 0 0 1], ...
    'STEStyle', 2, ...
    'Transparency', 0.3, ...
    'SignificanceLineOffset', 0.03, ...
    'SignificanceLinePadding', 0.05, ...
    'SignificanceLinePosition', 'below', ...
    'XLimits', [-0.1 3]);

sigIdx1 = SIGN{2,4,3,3}; %brain network 1 - Old - accuracy
writematrix(sigIdx1, '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/7.xlsx', 'Sheet', 1, 'Range', 'A1');

sigIdx1 = SIGN{2,4,2,3}; %brain network 1 - Old - accuracy
writematrix(sigIdx1, '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/8.xlsx', 'Sheet', 1, 'Range', 'A1');


%% 6) NT4 - brain network 2 - RTs 

sigIdx1 = SIGN{2,5,2,2};
% Apply to each cell in the input
win1 = cluster_sig_indices(sigIdx1, time);
SignificantWindows = {
    win1;
};

ConditionXGroupColors = {
    [0 0 0];
};
% Use the function
BROADNESS_Plot_ActivationTimeseries(data(2,:,:,5), time, ...
    'ConditionLabels', {'NT4'}, ...
    'ConditionXGroupColors', ConditionXGroupColors, ...
    'SignificantWindows', SignificantWindows, ...
    'SignificanceStyle', 'line', ...
    'SignificanceColors', [0 0 1], ...
    'STEStyle', 2, ...
    'Transparency', 0.3, ...
    'SignificanceLineOffset', 0.03, ...
    'SignificanceLinePadding', 0.05, ...
    'SignificanceLinePosition', 'below', ...
    'XLimits', [-0.1 3]);

sigIdx1 = SIGN{2,5,2,3}; %brain network 1 - Old - accuracy
writematrix(sigIdx1, '/Users/au550322/Documents/AarhusUniversitet/MattiaRosso/Paper_BROADNESS_PCA/AdvancedScience/Revision_I/Code/10.xlsx', 'Sheet', 1, 'Range', 'A1');

%%

%% TESTING THE FUNCTION - NOT DIRECTLY RELEVANT FOR THE MANUSCRIPT, BUT MEANINGFUL FOR CODING DEVELOPMENT ASSOCIATED TO IT

BROADNESS_Plot_ActivationTimeseries(data, time)

%% groups with user-defined colors

% Define a color for each condition × group
% Each cell corresponds to a condition
% Each row inside the matrix corresponds to a group
% Format: { [group1_color; group2_color], ... }

% Define color scheme for each condition × group
ConditionXGroupColors = {
    [0.8 0.2 0.2; 0.5 0.1 0.1];   % Condition 1: light red / dark red
    [0.2 0.4 0.8; 0.1 0.2 0.5];   % Condition 2: light blue / dark blue
    [0.4 0.7 0.2; 0.2 0.5 0.1];   % Condition 3: light green / dark green
    [0.7 0.3 0.7; 0.4 0.2 0.5];   % Condition 4: light purple / dark purple
    [0.2 0.7 0.7; 0.1 0.5 0.5];   % Condition 5: light cyan / dark cyan
};

% Call the function
BROADNESS_Plot_ActivationTimeseries(data, time, ...
    'Groups', {[1:30], [34:81]}, ...
    'GroupLabels', {'Group A', 'Group B'}, ...
    'ConditionLabels', {'Old', 'New', 'Cond 3', 'Cond 4', 'Cond 5'}, ...
    'ConditionXGroupColors', ConditionXGroupColors, ...
    'SignificantWindows', {
        {[0.4 0.6], [0.9 1.1]};                  % Level 1
        {[0.2 0.5]};                             % Level 2
        {[0.3 0.4], [0.9 1.3], [1.2 1.5]}        % Level 3
    }, ...
    'SignificanceColors', [
        1 0 0;       % red      → for level 1
        0 0.5 1;     % blue     → for level 2
        0 0.8 0.2    % green    → for level 3
    ], ...
    'SignificanceStyle', 'line', ...
    'SignificanceLineOffset', 0.03, ...
    'SignificanceLinePadding', 0.05, ...
    'SignificanceLinePosition', 'above', ...
    'STEStyle', 2, ...
    'Transparency', 0.2, ...
    'XLimits', [-0.1 3]);

%%

