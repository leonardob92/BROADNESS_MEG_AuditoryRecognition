function BROADNESS_Visualizer(BROADNESS, Options)

% ========================================================================
%  BROADBAND BRAIN NETWORK ESTIMATION VIA SOURCE SEPARATION (BROADNESS) TOOLBOX
%  VISUALIZER
% ========================================================================
%
%  Please cite the first BROADNESS paper:
%  Bonetti, L., Fernandez-Rubio, G., Andersen, M. H., Malvaso, C., Carlomagno,
%  F., Testa, C., Vuust, P, Kringelbach, M.L., & Rosso, M. (2025). Advanced Science.
%  BROAD-NESS Uncovers Dual-Stream Mechanisms Underlying Predictive Coding in Auditory Memory Networks.
%  https://doi.org/10.1002/advs.202507878
%
% ========================================================================
%
%
%  This script visualizes:
%
%  - Dynamic brain activity map
%    Showing the activity of each brain voxel (or channel) over time using
%    the imagesc function. This simply provides a useful depiction of the
%    original data, but is not related to the actual BROADNESS computations.
%
%  - Network prominence
%    Quantified as % of variance explained based on the associated  
%    eigenvalues.
%
%  - Network time series
%    Time series of each brain network (and, if provided, experimental
%    condition or whatever the user stored in the 3rd dimension of the
%    original data matrix)
%
%  - Spatial activation patterns of the networks  
%    - All requested networks are plotted together in a 3D visualization of the brain.  
%    - NIFTI files are generated for each network, allowing further  
%      inspection in FSLeyes or similar software for visualization.
%    - Excel files are generated with thresholded activation patterns (mean plus 1SD)
%      for each brain network. 
%
% ------------------------------------------------------------------------
%  INPUT ARGUMENTS:
% ------------------------------------------------------------------------
%  - BROADNESS                                          : Structure outputted by BROADNESS_NetworkEstimation function.
%                                                         We recommend not to modify this structure before giving it
%                                                         as input to the current function (BROADNESS_Visualizer).
%      - BROADNESS.Variance_BrainNetworks               : Normalized eigenvalues (variance explained in % points).
%      - BROADNESS.Significant_BrainNetworks            : Significant brain networks according to Monte-Carlo simulations (MCS).   
%      - BROADNESS.ActivationPatterns_BrainNetworks     : Spatial activation patterns used to generate NIFTI files and
%                                                         the 3D brain plot.
%      - BROADNESS.TimeSeries_BrainNetworks             : Time series of the brain networks (independently for each conditions and participants
%                                                         if original data was provided for each condition and participant).
%                                                         Note that if data was provided for each participant,the time series plots show
%                                                         mean across participants and standard errors.
%      - BROADNESS.VariancePermutations                 : Variance explained by PCA on permuted data; only if MCS was computed.
%      - BROADNESS.Time                                 : Vector with time in seconds. 
%      - BROADNESS.OriginalData                         : Original data matrix. 
%
%
%  - Options                  : Structure containing information for plotting.
%      - Options.WhichPlots   : Binary vector to indicate which plots should be produced (e.g. Options.WhichPlots = [0 1 0 1 1])
%                               according to the following order:
%                               - 1)Dynamic brain activity map of the original data
%                               - 2)Variance explained by the networks
%                               - 3)Time series of the networks
%                               - 4)Activation patterns of the networks (3D)
%                               - 5)Activation patterns of the networks (nifti images)
%                               Default: All plots will be generated.
%      - Options.name_nii     : path plus name for nifti images to be produced (one for each brain network)
%                               (e.g. 'YOUR_OWN_PATH/'
%                               then: 'Brain_Network_#_PROGRESSIVE NUMBER' will be automatically
%                               added to the name of the nifti image)
%      - Options.MNI_coords   : MNI coordinates provided in the same order as your data
%                               (N x 3, where N is the brain voxel number)
%      - Options.ncomps       : Components (networks) indices to be plotted in all plots but Variance plot
%                               (e.g. [1:5] for first 5 components or [2 5] for components 2 and 5).
%                               If the field is not provided, the default is to plot:
%                               - the MCS significant networks (if MCS was computed)
%                               - the first 5 networks
%      - Options.ncomps_var   : Number (!) of components (networks) to plot in the Variance plot.
%                               Default: first 20 components.
%      - Options.Labels       : If the original data matrix is 3D, here you can provide the labels
%                               of the experimental conditions.
%                               Cell array containing characters, e.g. Options.Labels = {'Cond 1';'Cond 2';'Cond 3'}.
%                               Default: 'Condition X', where X is a progressive number.
%      - Options.color_PCs    : Array with RGB color for PCs (e.g. [1 0 1; 1 1 0; 0.5 0.6 0.2]).
%                               If not supplied, default colors will be provided. 
%      - Options.color_conds  : Array with RGB color for experimental conditions (e.g. [1 0 1; 1 1 0; 0.5 0.6 0.2]).
%                               If not supplied, default colors will be provided. 
%
%
%
%
%  NOTE 1: While the 3D plot produced by this function  
%  is a convenient way to quickly inspect the topographies of  
%  multiple networks at once, it is recommended to use the  
%  NIFTI files for an accurate depiction of the networks' topographies.  
%  NIFTI files can be visualized using FSLeyes or an equivalent 
%  software.  
%
%  NOTE 2: 3D plot is supported for any brain MNI space (e.g. 1,2,8mm, etc).
%          The nifti images are currently supported only for 8mm.
%          If you need nifti images in a different space, feel free to contact us. 
%

% ------------------------------------------------------------------------
%  AUTHORS:
%  Leonardo Bonetti & Mattia Rosso
%  leonardo.bonetti@clin.au.dk; leonardo.bonetti@psych.ox.ac.uk
%  mattia.rosso@clin.au.dk
%  Center for Music in the Brain, Aarhus University
%  Centre for Eudaimonia and Human Flourishing, Linacre College, University of Oxford
%  Aarhus (DK), Oxford (UK), Bologna (Italy), Updated version 06/07/2025
%
% ========================================================================


% NOTE: we acknowledge the NIFTI Toolbox, which is used by FREQNESS to generate nifti images.
% Jimmy Shen (2025). Tools for NIfTI and ANALYZE image
% (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image)
% MATLAB Central File Exchange. Retrieved July 05, 2025.







%% Controlling inputs and parsing some data/information

data = BROADNESS.OriginalData;
time = BROADNESS.Time;
ActPat = BROADNESS.ActivationPatterns_BrainNetworks;
TimeSeries = BROADNESS.TimeSeries_BrainNetworks;

% Compute mean and standard deviation if data is provided for single participants
sz = size(data);
non_singleton_dims = sum(sz > 1); %trick to get if the matrix is a vector
if non_singleton_dims == 4 %data provided for single participants
    TimeSeries_stde = std(TimeSeries,[],4) ./ sqrt(size(TimeSeries,4));    
    TimeSeries = mean(TimeSeries,4);
end

% Detect whether this BROADNESS struct is from PCA (has variance fields)
isPCA = isfield(BROADNESS,'Variance_BrainNetworks') && ...
    isnumeric(BROADNESS.Variance_BrainNetworks) && ...
    ~isempty(BROADNESS.Variance_BrainNetworks);

if ~isfield(Options,'WhichPlots') % if request of which plots should be prepared is not provided
    Options.WhichPlots = ones(1,5);  % Assigning default
end

% Checking if MNI coordinates are provided
if Options.WhichPlots(4) == 1 && ~isfield(Options,'MNI_coords')
    error('MNI coordinates must be provided for 3d plotting in brain template.. (Options.WhichPlots = [0 0 0 1 0])')
end

% Checking if path and name to nifti file are provided
if Options.WhichPlots(5) == 1 && ~isfield(Options,'name_nii')
    error('Path and name to nifti file must be provided for saving nifti images.. (Options.WhichPlots = [0 0 0 0 1])')
end

if isfield(Options,'Labels') % if user provided labels, they are extracted for later plotting purposes
    Labels = Options.Labels;
else % otherwise default is generated
    Labels = cell(1,size(data, 3));
    for condi = 1:size(data, 3)
        Labels(condi) = {['Cond ' num2str(condi)]};
    end
end

if isfield(Options,'ncomps_var') % if number of components for variance plot are provided, they are extracted here
    ncomps_var = Options.ncomps_var;
else %otherwise default is assigned
    ncomps_var = 20;
end


%assigning number of components to be plotted
if isfield(Options,'ncomps')
    ncomps = Options.ncomps;
elseif isPCA && isfield(BROADNESS,'Significant_BrainNetworks') && ~ischar(BROADNESS.Significant_BrainNetworks)
    ncomps = BROADNESS.Significant_BrainNetworks;   % PCA + MCS
else
    % Fallback for ICA or PCA without MCS: first up to 5 components
%     ncomps = 1:min(5, size(ActPat,2));
    ncomps = 1:5;%size(ActPat,2);
end
if length(ncomps) > size(TimeSeries,1)
    ncomps = 1:size(TimeSeries,1);
end
% if isfield(Options,'ncomps') %if components indices to be plotted are provided
%     ncomps = Options.ncomps; %extracting them
% elseif ~ischar(BROADNESS.Significant_BrainNetworks) %otherwise if MCS was previously computed
%     ncomps = BROADNESS.Significant_BrainNetworks; %extracting indices of significant PCs
% else %otherwise default is assigned
%     ncomps = 1:5;
% end

% for later plotting solutions
thresh_nsdt = 1; % how many std away from the mean, for thresholding the visualization




% Checking if colors are provided for PCs and experimental conditions, otherwise assigning default
if isfield(Options,'color_PCs') %PCs
    col_comp = Options.color_PCs;
    if size(col_comp,1) < length(ncomps)
        warning('You have more PCs than supplied colors.. thus assigning colors by default')
        col_comp = .8*cool(length(ncomps)); %color mapping for PCs
    end
else
    col_comp = .8*cool(length(ncomps)); %color mapping for PCs
end
if isfield(Options,'color_conds') %conditions
    col_cond = Options.color_conds;
    if size(col_cond,1) < size(data,3)
        warning('You have more experimental conditions than supplied colors.. thus assigning colors by default')
        col_cond = .8*lines(size(data,3)); %color mapping for PCs
    end
else
    col_cond = .8*lines(size(data,3)); %color mapping for conditions
end

%% 1)Dynamic brain activity map of the original data

if Options.WhichPlots(1) == 1
    
    % Plotting data using imagesc (sort of raster plots)
    disp('Generating dynamic brain activity plots...');

    data = mean(data,4); % Average across participants (if single-participant data was provided)
    % Computing global min and max
    MAX = max(data(:));
    MIN = min(data(:));
    % Computing symmetric scaling limit based on max absolute value
    maxAbs = max(abs([MIN, MAX]));
    cLim = [-maxAbs, maxAbs];
    
    % Plotting each experimental condition with consistent color scaling
    for condi = 1:size(data, 3) % Over experimental condition (or whatever the user placed in the 3rd dimension of the original data)
        figure
        imagesc(time, 1:size(data,1), data(:,:,condi));
        set(gcf, 'Color', 'w');
        colorbar;
        caxis(cLim);  % Applying symmetric color scale
        title(['Dynamic Brain Activity - ' Labels{condi}], 'FontWeight', 'bold', 'FontSize', 14);
        xlabel('Time (s)'); ylabel('Brain sources');
    end
end

    
%% 2)Variance explained by the networks

if Options.WhichPlots(2) == 1
    
    if ~isPCA
        warning('Variance field not found (likely because you provided ICA results).. Skipping variance plot');
    else
        
        disp('Generating variance explained plot...');
        
        figure
        % Plotting actual data with line + marker (e.g. star marker)
        plot(BROADNESS.Variance_BrainNetworks(1:ncomps_var), ...
            '-*', ...                         % line + star marker
            'DisplayName', 'Data', ...
            'LineWidth', 1.5, ...
            'MarkerSize', 6);
        hold on
        
        % If MCS permutations available, plotting them too
        if ~ischar(BROADNESS.Significant_BrainNetworks)
            plot(BROADNESS.VariancePermutations(1:ncomps_var), ...
                '-o', ...                     % line + circle marker
                'DisplayName', 'Random', ...
                'LineWidth', 1.5, ...
                'MarkerSize', 5);
        end
        
        % Additional setting for the plot
        grid minor
        legend('show', 'Location', 'northeast')
        set(gcf, 'Color', 'w')
        title('Variance Explained by Principal Components', 'FontWeight', 'bold', 'FontSize', 14);
        xlabel('Component #');
        ylabel('% Variance Explained');
    end

end


%% Time series of the networks

if Options.WhichPlots(3) == 1
    
    disp('Generating time series plots for brain networks...');
    
    if ~isPCA
        TimeSeries = permute(TimeSeries,[2 1 3]);
        if non_singleton_dims == 4 % data provided for single participants
            TimeSeries_stde = permute(TimeSeries_stde,[2 1 3]);
        end
    end
    
    for compi = 1:length(ncomps) %over selected PCs

        figure; hold on; set(gcf,'color','w');
        grid minor
        t = time; t = t(:);
        for condi = 1:size(data,3)
            mu = TimeSeries(:,ncomps(compi),condi);
            if non_singleton_dims == 4 % data provided for single participants
                se = TimeSeries_stde(:,ncomps(compi),condi);
                mu = mu(:); se = se(:);
                fill([t; flipud(t)], [mu+se; flipud(mu-se)], col_cond(condi,:), 'FaceAlpha',0.25, 'EdgeColor','none', 'HandleVisibility','off');
            end
            plot(t, mu, 'Color', col_cond(condi,:), 'LineWidth',2, 'DisplayName', Labels{condi});
        end
        xlim([t(1) t(end)]);
        legend('show')
        box on

        if isPCA
            title(['Time Series - Brain Networks # ' num2str(ncomps(compi)) ...
                ' - Var ' num2str(BROADNESS.Variance_BrainNetworks(ncomps(compi)))], ...
                'FontWeight','bold','FontSize',14);
        else
            title(['Time Series - Brain Networks # ' num2str(ncomps(compi))], ...
                'FontWeight','bold','FontSize',14);
        end
        %         title(['Time Series - Brain Networks # ' num2str(ncomps(compi)) ' - Var ' num2str(BROADNESS.Variance_BrainNetworks(ncomps(compi)))], 'FontWeight', 'bold', 'FontSize', 14);
        
        xlabel('Time (s)'); ylabel('Component Amplitude');
    end
end


%% 4)Activation patterns of the networks (3D) - in brain template

if Options.WhichPlots(4) == 1
    
    disp('Generating 3D topographic plots of brain networks...');
    
    % Some default settings
    skipper  = 1;          % parameter for donwsampling the visualization
    scale_size = 100;      % scaling factor for activation patterns in the brain
    
    % Opening figure
    openfig('BrainTemplate_GT.fig')
    hold on
    legend_handles = gobjects(1, length(ncomps)); % Preallocate legend handles
    
    for compi = 1:length(ncomps) %over selected PCs
        % Assigning temporary activation pattern to plot
        pat2plot = ActPat(:,ncomps(compi));
        pat2plot( pat2plot < mean(pat2plot)+thresh_nsdt*std(pat2plot) ) = nan;  % apply threshold
        pat2plot(isnan(Options.MNI_coords(:,1))) = nan;
        
        % Assigning temporary activation pattern to plot
        mni2plot = Options.MNI_coords;
        mni2plot(isnan(pat2plot(:,1)),:) = []; % clear from nans
        pat2plot(isnan(pat2plot)) = [];   % repeat for activation patterns
        % Re-scaling pattern to plot
        %         pat2plot = pat2plot/max(pat2plot);  % option #1
        pat2plot = (pat2plot-min(pat2plot))./(max(pat2plot)-min(pat2plot)).*(1-.01) + .01;  % option #2
        
        % Plotting in 3D
        for voxi = 1:skipper:length(pat2plot)
            plot3( mni2plot(voxi,1), mni2plot(voxi,2), mni2plot(voxi,3), '.', 'Color', col_comp(compi,:), 'MarkerSize', scale_size * pat2plot(voxi) );
            hold on
        end
        
        % Using NaN to avoid plotting actual points, just storing color info
        legend_handles(compi) = plot3(nan, nan, nan, '.', 'Color', col_comp(compi,:), 'MarkerSize', 12);    
    end
   
    % Adjusting legend properties
    legend(legend_handles, arrayfun(@(x) sprintf('Component %d', x), ncomps, 'UniformOutput', false), ...
        'FontSize', 14, ...  % Increasing font size
        'Location', 'northeastoutside'); % Moving legend outside the plot area
    set(legend_handles, 'MarkerSize', 20); % Increasing marker size in legend
    rotate3d on; axis off; axis vis3d; axis equal % Allowing rotation and removing axes
    set(gcf, 'Color', 'w'); % Setting figure background to white
    camlight; lighting gouraud;
    
    title('3D Spatial Patterns of Brain Networks', 'FontSize', 16, 'FontWeight', 'bold')
end


%% 5)Activation patterns of the networks (brain nifti images)

if Options.WhichPlots(5) == 1
    
    disp('Generating and saving NIFTI images of brain network activation patterns...');

    % Creating directory for storing BROADNESS NIFTI output files
    nifti_path = [Options.name_nii '/BROADNESS_Output/BROADNESS_nifti'];
    mkdir(nifti_path)
    
    % Loading template
    template_nii = load_nii('MNI152_8mm_brain_diy.nii.gz');
    
    % Geting template image data and initializing an empty volume
    nii_data = template_nii.img;
    nii_data(:) = 0;  % Setting all voxels to zero
    nii_data = double(nii_data);
    
    % Extracting affine transformation matrix from srow_x, srow_y, srow_z
    affine = [template_nii.hdr.hist.srow_x;
        template_nii.hdr.hist.srow_y;
        template_nii.hdr.hist.srow_z;
        0 0 0 1]; % Appending [0 0 0 1] to make it 4x4
    
    % Extracting MNI coordinates
    MNI_coords = Options.MNI_coords;
    
    for compi = 1:length(ncomps) % Over components given as input
        
        % Assigning temporary activation pattern to plot
        pat2plot = ActPat(:,ncomps(compi));
%         pat2plot( pat2plot < mean(pat2plot)+thresh_nsdt*std(pat2plot) ) = 0;  % apply threshold
        pat2plot( abs(pat2plot) < mean(abs(pat2plot))+thresh_nsdt*std(abs(pat2plot)) ) = 0;  % apply threshold
        
        % -----------------------
        % Create table of non-zero activations (index, MNI coords, activation)
        % -----------------------
        % Find voxels with non-zero activation after thresholding
        voxelIdxList = find(pat2plot ~= 0);
        
        if ~isempty(voxelIdxList)
            % Collect activations for those voxels
            NetworkActPats = pat2plot(voxelIdxList);
            
            % Choose MNI coordinates to use for the table:
            % Prefer the per-component MNI_coords (Options.MNI_coords) if available,
            % otherwise fall back to coordinates.MNI8 (global grid). The typical case
            % is that MNI_coords rows correspond to pat2plot entries, so we index into them.
            if exist('MNI_coords','var') && ~isempty(MNI_coords)
                MNIcoords_for_table = MNI_coords(voxelIdxList, :);
            else
                % If MNI_coords not provided, try to use coordinates.MNI8
                % NOTE: This assumes voxelIdxList indexes into coordinates.MNI8 appropriately.
                MNIcoords_for_table = coordinates.MNI8(voxelIdxList, :);
            end
            
            % Compose table data: running index within the list, coords, activations
            excel_data = [(1:size(NetworkActPats,1))', MNIcoords_for_table, NetworkActPats];
            
            % Create headers: Index, X, Y, Z, Activation
            tableHeaders = [{'Index','X','Y','Z','Activation'}];
            
            % Convert to table
            tbl = array2table(excel_data, 'VariableNames', tableHeaders);
            
            % Save the table to an Excel file in the nifti_path
            % Filename includes component number (ncomps(compi))
            table_filename = fullfile(nifti_path, sprintf('ActivationTable_BrainNetwork_%d.xlsx', ncomps(compi)));
            try
                writetable(tbl, table_filename);
                disp(['Saved activation table: ' table_filename]);
            catch ME
                warning('Could not write activation table to Excel: %s\nFalling back to .mat save. Error: %s', table_filename, ME.message);
                save(fullfile(nifti_path, sprintf('ActivationTable_Component_%d.mat', ncomps(compi))), 'tbl');
            end
        else
            disp(['No suprathreshold voxels for component ' num2str(ncomps(compi)) '; skipping table creation.']);
        end
        
        % -----------------------
        % Continue with voxel -> NIfTI mapping
        % -----------------------

        num_points = size(MNI_coords, 1);
        voxel_coords = zeros(num_points, 3);
        
        % Converting MNI Coordinates to Voxel Indices
        for ii = 1:num_points
            coord = [MNI_coords(ii, :) 1];  % Add homogeneous coordinate
            voxel = affine\coord';% inv(affine) * coord';  % Convert to voxel space
            voxel_coords(ii, :) = (voxel(1:3)); % Extract rounded voxel indices
            
            % Adjusting for voxel center vs edge (half voxel shift)
            voxel_coords(ii, :) = voxel_coords(ii, :) + 1;  % Subtract 1 voxel (adjust for 8mm shift)
            
        end
        
        % Assigning activation patterns to the corresponding voxels
        for ii = 1:num_points
            x = voxel_coords(ii, 1);
            y = voxel_coords(ii, 2);
            z = voxel_coords(ii, 3);
            
            % Ensuring indices are within image boundaries
            dims = size(nii_data);
            if x > 0 && x <= dims(1) && y > 0 && y <= dims(2) && z > 0 && z <= dims(3)
                %             if all([x, y, z] > 0) && all([x, y, z] <= size(nii_data))
                nii_data(x, y, z) = pat2plot(ii);%ActPat(ii, ncomps(compi));
            end
        end
        
        template_nii.img = nii_data; % Storing matrix within the image structure
        
        % Creating a NIFTI image from the 3D data matrix (8 mm resolution)
        nii = make_nii(nii_data, [8 8 8]);
        nii.img = nii_data;  % Storing matrix within image structure
        nii.hdr.hist = template_nii.hdr.hist;  % Copying header information from mask
        
        % Displaying saving progress
        disp(['Saving NIFTI images - PC ' num2str(ncomps(compi))])
        
        % Saving the NIFTI file
        if isPCA
            save_nii(nii, [nifti_path '/PCA_ActivationPattern_BrainNetwork_#' num2str(ncomps(compi)) '.nii']);
        else
            save_nii(nii, [nifti_path '/ICA_ActivationPattern_BrainNetwork_#' num2str(ncomps(compi)) '.nii']);
        end
    end
end

end
