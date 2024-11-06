%%

%% BROADNESS (USING PCA)

% 1) Please, read the function documentation for detailed, additional information.
% 2) The function is available at the following link: https://github.com/leonardob92/LBPD-1.0
% 3) Please, also note that if you want to plot the brain networks in nifti images
%    you should download the full LBPD ot the following toolboox: https://se.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

S = []; %opening a Matlab structure
S.H = dum(:,:,1:5); %providing the data (space x time x conditions); conditions here will be averaged in the functions and PCA computed on their average
S.permnum = 1; %number of permutations if you wish to use BROADNESS with the Monte-Carlo simulations approach
S.fig_l = 1; %if you want to plot some figures
S.sign_eig = '0'; %option to handle the sign of the eigenvector elements
S.namenii = 'YOUR_PATH/YOUR_NAME'; %path and name for the nifti files which the function provides to show the brain networks (PCA components)
S.time = time; %vector with time in seconds
S.rand_l = 1; %randomisation strategy (see more details in the function description)
S.onefig = 0; %option for additional figure

[ OUT ] = PCA_LBPD( S ) %actual function

%%