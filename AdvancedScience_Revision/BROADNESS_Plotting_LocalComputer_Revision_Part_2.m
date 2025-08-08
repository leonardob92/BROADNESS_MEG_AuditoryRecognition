%% 
%% RQA- PHASE SPACE PLOT
clear all
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\PhaseSpace_analysis\phase_space_plot.mat')
export_l =0;
timemin = 1; % starting time point
timemax = 651; % ending time point
% prova in locale, per ora qui in pdf
scatsize = 30;
nscaz = size(phase_space{1,1},1);
cmap = jet(nscaz);
X = 30; %time-steps
for cond = 1:5 % over conditions
    phase_temp = zeros(size(phase_space{1,1},1),size(phase_space{1,1},2), size(phase_space,2));
    for sub = 1:size(phase_space,2)
        phase_temp(:,:,sub) = phase_space{cond,sub}(:,:);
    end

    figure
    scatter(squeeze(mean(phase_temp(:,1,:),3)),squeeze(mean(phase_temp(:,2,:),3)),scatsize,cmap,'filled')
    xlim([-2000 1000])
    ylim([-2000 1000])
    set(gcf,'color','w')
    % legend('show')
    grid minor
    box on
    c = colorbar;
    colormap('jet')
    bumba = time_sel(timemin:X:timemax);
    c.Ticks = linspace(0,1,length(bumba));
    c.TickLabels = bumba;
    xlabel('Comp 1')
    ylabel('Comp 2')
    title(['Condition ' num2str(cond)])
    if export_l ==1
    exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\PhaseSpace_analysis\phase_space_condition' num2str(cond) '.pdf'],'Resolution',300)
    end
end
%% RQA- BOX AND VIOLIN PLOT
clear all
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\PhaseSpace_analysis\metrics_from_RP.mat')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\PhaseSpace_analysis\significant_pvals_on_metrics.mat')
labels = {'recurrence_rate','mean_diag_line', 'determinism', 'entropy', 'trapping_time', 'laminarity', 'max_laminar_state', 'divergence' };
export_l =0;
% Define color scheme for each condition × group
ConditionXGroupColors = {
    [1 0 0];   % Condition 1: light red / dark red
    [0.3686 0.6314 0.7412];   % Condition 2: light blue / dark blue
    [0.1882 0.4902 0.8118];   % Condition 3: light green / dark green
    [ 0.0784 0.1569 0.5255];   % Condition 4: light purple / dark purple
    [0 0 0];   % Condition 5: light cyan / dark cyan
};
for cc=1:size(sig_pvals,1) % over comparisons
    for mm = 1:size(sig_pvals,2) % over metrics
        if sig_pvals(cc,mm) ~= 1
            metric = [metrics(1,:,mm) metrics(cc+1,:,mm)];
            cond = [ones(1, 83), (cc+1) * ones(1, 83)];
            figure
            boxchart(cond', metric', 'BoxFaceColor', [0.1882 0.4902 0.8118], 'MarkerColor', [0.1882 0.4902 0.8118], 'LineWidth', 1.5);  % 'flat' = colori personalizzati
            grid minor
            unique_conds = unique(cond);  % es. [1 5]
            xticks(unique_conds)
            set(gcf,'color','w')
            box on
            if export_l ==1
                exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\PhaseSpace_analysis\boxplot_' labels{mm} '_cond1VS' num2str(cc+1) '.pdf'],'Resolution',300)
            end
                        % VIOLIN PLOT
            figure
            h_violin = violinplot(cond', metric');
            % Imposta colori e line width simili
            for v = 1:length(h_violin)
                h_violin(v).ViolinColor = [0.1882 0.4902 0.8112];
                h_violin(v).MedianLineWidth = 1.5;
                h_violin(v).BoxColor = [0 0 0]; % contorno box nero (opzionale)
                h_violin(v).EdgeColor = [0 0 0]; % contorno violino nero (opzionale)
            end
            
            grid minor
            xticks(unique_conds)
            xticklabels(arrayfun(@num2str, unique_conds, 'UniformOutput', false))
            set(gcf,'color','w')
            box on
        end
    end
end

%% SPATIAL GRADIENT - ELBOW PLOT
clear all
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\SpatialGradient_analysis\ElbowPlot_SUM.mat')
export_l = 1; % 1 to export the fig
colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];
figure;
plot(SUM, 'Color',colorline(3,:),'Linewidth', 3);
hold on;
grid minor
set(gcf,'color','w')
box on

if export_l ==1
    exportgraphics(gcf,'C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\SpatialGradient_analysis\ElbowPlot_20k.pdf','Resolution',300)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SPATIAL GRADIENT - SCATTER PLOT
clear all
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\SpatialGradient_analysis\ActivationPattern_thresholded.mat')
export_l = 0; % 1 to save the figure
nvoxels = size(SO,1);
scatsize = 30; % dimension of scatter points
cmap = jet(nvoxels);
figure
scatter(SO(:,1)', SO(:,2)', scatsize, cmap, 'filled')
set(gcf, 'color', 'w')
grid minor
box on
c = colorbar;
colormap('jet')
numTicks = 10;
c.Ticks = linspace(0,1,numTicks);
c.TickLabels = round(linspace(1, nvoxels, numTicks));% corresponding voxel indices

if export_l ==1
    exportgraphics(gcf,'C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\SpatialGradient_analysis\ScatterPlot_component1VS2.pdf','Resolution',300)
end

%% SPATIAL GRADIENT - BEST K HISTOGRAM
clear all
export_l = 1;
colorline = [1 0 0; 0.3686 0.6314 0.7412; 0.1882 0.4902 0.8118; 0.0784 0.1569 0.5255; 0 0 0];
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\SpatialGradient_analysis\eval_optimal_K_1000_repetitions.mat')
figure
histogram(eva_optimal, 'FaceColor',  colorline(3,:), 'FaceAlpha',0.8)
set(gcf, 'color', 'w')
grid minor
box on

if export_l ==1
    exportgraphics(gcf,'C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\SpatialGradient_analysis\EvaOptimal_histogram_1000repetitions.pdf','Resolution',300)
end
%% SPATIAL GRADIENT - CLUSTER PLOT WITH K = 8
clear all
outputdir = 'C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\SpatialGradient_analysis';
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\SpatialGradient_analysis\ActivationPattern_thresholded.mat')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\SpatialGradient_analysis\Cluster_idx.mat')
export_l = 1; % 1 to save the figure
clusterN = 20;
optimalK = 8; %optimal number of clusters
for clustern = 1:clusterN
    IDX = idx_ks(:,clustern);
    scatsize = 30;
    bumba =clustern;
    cmap = jet(bumba);
    cmap = cmap.*0.9;
    figure
    for ii = 1:clustern
    
        scatter(SO(IDX == ii,1)',SO(IDX == ii,2)', scatsize, 'MarkerFaceColor', cmap(ii,:), 'MarkerEdgeColor', cmap(ii,:))
        hold on
        % xlim([-2000 1000])
        % ylim([-2000 1000])
    
    end
    set(gcf,'color','w')
    % legend('show')
    grid minor
    box on
    % Set discrete colormap
    colormap(cmap)
    c = colorbar;
    % colormap('jet')
    % IMPORTANT: map colorbar limits to number of clusters
    clim([0.5 clustern + 0.5])         % <-- fundamental: sets range of color values
    c.Ticks = 1:clustern;
    c.TickLabels = 1:clustern;
    % c.Ticks = linspace(0,1,bumba);
    % c.TickLabels = 1:bumba;
    if export_l ==1
        if clustern == optimalK
            exportgraphics(gcf,[outputdir '\ScatterPlot_component1VS2_CLUSTERS_OptimalK_' num2str(clustern) '.pdf'],'Resolution',300)
        else
            exportgraphics(gcf,[outputdir '\ScatterPlot_component1VS2_CLUSTERS_K_' num2str(clustern) '.pdf'],'Resolution',300)
        end
    end
end

%% ICA-CORRELATION COEFFICIENTS PLOT
% loading data
clear all
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\correlations_coeffs_26maxcomponents.mat')
export_l = 0; %1 to save the plot
% Numero di curve
n = size(sorted_corr,1); % number of ICA components


% Genera 25 tonalità di rosso da chiaro a scuro o viceversa
reds = [ linspace(1, 0.2, n)', zeros(n,1), zeros(n,1)]; % da rosso scuro a chiaro
for pca_comps = 1:size(sorted_corr,3) % over PCA components
    figure;
    for ica_comps = 1:n % over ICA components 
        temp_corr = zeros(ica_comps + 1, size(sorted_corr,2));
        for cond = 1: size(sorted_corr,2)
            temp_corr(:,cond) = sorted_corr{ica_comps, cond, pca_comps};
        end
        temp_corr_av = mean(temp_corr, 2);

        plot(temp_corr_av, 'Color',reds(ica_comps,:),'Linewidth', 1.5)
        hold on
    end
    xlim([0 28])
    ylim([0 1])
    grid minor
    set(gcf,'color','w')
    box on
    if export_l ==1
    exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\Correlation_VS_components_network' num2str(pca_comps) '.pdf'],'Resolution',300)
    end
end
%% TIME SERIES WITH SIGNIFICANCE INTERVAL - TWO NETWORKS, OUT OF 14 TOTAL COMPONENTS
% loading data
clear all
export_l =1;
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\best_ICA_components_14.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\fdr_results_bestICA_14components.mat');
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time
time = time_sel(1:1026);
data_temp = permute(best_ICA,[1 2 4 3]);
for cc = 1:size(data_temp,1)
    data = data_temp(cc,:,:,:);
    
    % Define color scheme for each condition × group
    ConditionXGroupColors = {
        [1 0 0];   % Condition 1: light red / dark red
        [0.3686 0.6314 0.7412];   % Condition 2: light blue / dark blue
        [0.1882 0.4902 0.8118];   % Condition 3: light green / dark green
        [ 0.0784 0.1569 0.5255];   % Condition 4: light purple / dark purple
        [0 0 0];   % Condition 5: light cyan / dark cyan
    };
    SignificanceColors = [
        0.3686 0.6314 0.7412;   % Condition 2: light blue / dark blue
        0.1882 0.4902 0.8118;   % Condition 3: light green / dark green
        0.0784 0.1569 0.5255;   % Condition 4: light purple / dark purple
        0       0       0      % Condition 5: black (or placeholder)
    ];
    
    fdr_results_1 = fdr_results(cc,:)';
    
    
    BROADNESS_plot_timeseries(data, time, 'ConditionXGroupColors',ConditionXGroupColors','SignificantWindows', fdr_results_1, 'SignificanceColors', SignificanceColors, 'SignificanceStyle', 'line','SignificanceLinePosition','below' ...
        ,  'SignificanceLineOffset', 0.05, 'SignificanceLinePadding', 0.05, 'YLimits', [-2, 2.8], 'XLimits', [-0.1, 3.4])
    %BROADNESS_plot_timeseries(data, time, 'ConditionXGroupColors', ConditionXGroupColors)
    if export_l ==1
    exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\ICA_timeseries_14components_network' num2str(cc) '.pdf'],'Resolution',300)
    end
end


%% ICA- TIME SERIES WITH SIGNIFICANCE INTERVAL - TWO NETWORKS, OUT OF 2 TOTAL COMPONENTS
% loading data
clear all
export_l =1;
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\best_ICA_components_2.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\fdr_results_bestICA_2components.mat');
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time
time = time_sel(1:1026);
data_temp = permute(best_ICA,[1 2 4 3]);
for cc = 1:size(data_temp,1)
    data = data_temp(cc,:,:,:);
    
    % Define color scheme for each condition × group
    ConditionXGroupColors = {
        [1 0 0];   % Condition 1: light red / dark red
        [0.3686 0.6314 0.7412];   % Condition 2: light blue / dark blue
        [0.1882 0.4902 0.8118];   % Condition 3: light green / dark green
        [ 0.0784 0.1569 0.5255];   % Condition 4: light purple / dark purple
        [0 0 0];   % Condition 5: light cyan / dark cyan
    };
    SignificanceColors = [
        0.3686 0.6314 0.7412;   % Condition 2: light blue / dark blue
        0.1882 0.4902 0.8118;   % Condition 3: light green / dark green
        0.0784 0.1569 0.5255;   % Condition 4: light purple / dark purple
        0       0       0      % Condition 5: black (or placeholder)
    ];
    
    fdr_results_1 = fdr_results(cc,:)';
    
    
    BROADNESS_plot_timeseries(data, time, 'ConditionXGroupColors',ConditionXGroupColors','SignificantWindows', fdr_results_1, 'SignificanceColors', SignificanceColors, 'SignificanceStyle', 'line','SignificanceLinePosition','below' ...
        ,  'SignificanceLineOffset', 0.05, 'SignificanceLinePadding', 0.05, 'YLimits', [-1.6, 1.5], 'XLimits', [-0.1, 3.4])
    %BROADNESS_plot_timeseries(data, time, 'ConditionXGroupColors', ConditionXGroupColors)
    % if export_l ==1
    % exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\ICA_timeseries_2components_network' num2str(cc) '.pdf'],'Resolution',300)
    % end
    if export_l ==1
    exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\ICA_timeseries_2components_network' num2str(cc) '_newlimit_.pdf'],'Resolution',300)
    end
end


%% TIME SERIES WITH SIGNIFICANCE INTERVAL - ALL NETWORKS, OUT OF 14 TOTAL COMPONENTS
% loading data
clear all
export_l =1;
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\all_ICA_components_14.mat');
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\fdr_results_allICA_14components.mat');
addpath('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara')
load('C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\TempChiara\TempChiara/time.mat'); %loading time
time = time_sel(1:1026);
data_temp = permute(icScores,[1 2 4 3]);
for cc = 1:size(data_temp,1)
    data = data_temp(cc,:,:,:);
    
    % Define color scheme for each condition × group
    ConditionXGroupColors = {
        [1 0 0];   % Condition 1: light red / dark red
        [0.3686 0.6314 0.7412];   % Condition 2: light blue / dark blue
        [0.1882 0.4902 0.8118];   % Condition 3: light green / dark green
        [ 0.0784 0.1569 0.5255];   % Condition 4: light purple / dark purple
        [0 0 0];   % Condition 5: light cyan / dark cyan
    };
    SignificanceColors = [
        0.3686 0.6314 0.7412;   % Condition 2: light blue / dark blue
        0.1882 0.4902 0.8118;   % Condition 3: light green / dark green
        0.0784 0.1569 0.5255;   % Condition 4: light purple / dark purple
        0       0       0      % Condition 5: black (or placeholder)
    ];
    
    fdr_results_1 = fdr_results(cc,:)';
    
    
    BROADNESS_plot_timeseries(data, time, 'ConditionXGroupColors',ConditionXGroupColors','SignificantWindows', fdr_results_1, 'SignificanceColors', SignificanceColors, 'SignificanceStyle', 'line','SignificanceLinePosition','below' ...
        ,  'SignificanceLineOffset', 0.05, 'SignificanceLinePadding', 0.05, 'YLimits', [-2, 2.8], 'XLimits', [-0.1, 3.4])
    %BROADNESS_plot_timeseries(data, time, 'ConditionXGroupColors', ConditionXGroupColors)
    if export_l ==1
    exportgraphics(gcf,['C:\Users\malva\OneDrive\Desktop\APPLIED_PHYSICS\Tesi\PCA\Broadness\reviews\AdvancedScience_Revision\ICA_analysis\all_ICA_timeseries_14components_network' num2str(cc) '.pdf'],'Resolution',300)
    end
end

