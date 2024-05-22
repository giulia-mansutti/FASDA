% Analysis of speckle frames acquired from a dysplastic nevus with crossed polarizers
% according to full-field amplitude speckle decorrelation angiography (FASDA) [1]
%
% Authors:  Giulia Mansutti, Néstor Uribe-Patarroyo
%
% GM:
% 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
% General Hospital, 40 Blossom Street, Boston, MA, USA;
% <gmansutti@mgh.harvard.edu>
%
% NUP: 
% 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
% General Hospital, 40 Blossom Street, Boston, MA, USA;
% <uribepatarroyo.nestor@mgh.harvard.edu>
%
% MGH Full-field amplitude speckle decorrelation angiography (FASDA) project
%
% References:
% [1] G. Mansutti, M. Villiger, B. E. Bouma, and N. Uribe-Patarroyo, 
% 'Full-field amplitude speckle decorrelation angiography', Biomedical
% Optics Express, May 2024 [submitted]
%
%
% Copyright and permission related to the usage of color maps:
%
% Copyright (c) 2014-2020 Peter Kovesi
% Centre for Exploration Targeting
% The University of Western Australia
% peter.kovesi at uwa edu au
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.


% Set variables related to directories of interest
addpath(genpath('Support Functions/'));     % Add path to support functions folder
dataPath = 'Data/Crossed Polarizers/';      % data folder's path
outputPath = 'Output/Crossed Polarizers/';  % output folder's path
[~, ~, ~] = mkdir(outputPath);              % create output folder
          
% Set variables related to speckle frames processing
filenameRoot = 'file'; % each frame is named as "file" followed by a number
nY = 1024;             % y-axis dimension of each speckle frame (i.e., image height, or number of rows)
nX = 1280;             % x-axis dimension of each speckle frame (i.e., image width, or number of columns)
tauVec = [0 1 15];     % vector containing the time lags of interest (expressed in number of frames) to compute the second-order auto-correlation function

% Set variables related to figures generation
savefigures = true;                 % Decide whether to save figures or not (figures are saved in the 'Output' folder) 
currFig = 0;                        % start generating figures from this index
colorMapFASDA = colorcet('L20');    % Colormap used to generate ASDI figures
colorMapBFI = colorcet('CBL3');     % Colormap used to generate BFI figures
speckleFrameIndex = 10;             % index of the sample speckle frame to visualize (can be any number between 1 and 500)
lowSNRpixelVal = 35;                % threshold value to define low SNR pixels
speckleFrameColorbarLim = [0 255];  % color bar limits for speckle frame figure
betaBFIColorbarLim = [25 65];       % color bar limits for beta_BFI figure
alphaASDIColorbarLim = [0.2 0.95];  % color bar limits for alpha_ASDI figure

% Set variables related to the spatial ensemble window
xEnsembleWindowHalfSize = 2;       % half of the window size along the x-dimension
xEnsembleWindowExp2Diameter = 2;   % diameter along the x-dimension
yEnsembleWindowHalfSize = 2;       % half of the window size along the y-dimension
yEnsembleWindowExp2Diameter = 2;   % diameter along the y-dimension

%% Load image data
fprintf('\t Loading speckle frames... ')

% get info about speckle frames in data folder
dirElements = dir(fullfile(dataPath, strcat(filenameRoot, '*.bmp'))); 
nFrames = size(dirElements, 1); % number of frames in data folder 

% create the matrix containing all acquired speckle frames
speckleMatrix = zeros(nY, nX, nFrames, 'single'); 
thisImage = 1;
for k = 1:nFrames
  speckleMatrix(:, :, thisImage) = imread(fullfile(dataPath, dirElements(k).name));
  thisImage = thisImage + 1;
end
fprintf('done.\n')

% Visualize a sample speckle frame 
fig1 = figure(currFig + 1); imagescnan(speckleMatrix(:, :, speckleFrameIndex), speckleFrameColorbarLim)
colorbar, colormap(gray(256)), axis image, title('sample speckle frame')
if savefigures
  filename = sprintf('sampleSpeckleFrame');  % file name to save sample speckle frame
  saveas(fig1,fullfile(outputPath, strcat(filename, '.png')));
end


%% Mask saturated pixels

% Define mask of saturated pixels as those with saturation for more than 1/8 of the frames
satPxMask = single(sum(speckleMatrix == 255, 3) <= nFrames / 8);
satPxMask(~satPxMask) = nan;

% Define low-signal pixels as those with a value lower than lowSNRpixelVal for more than 1/8 of the frames
lowSnrPxMask = single(sum(speckleMatrix <= lowSNRpixelVal, 3) <= nFrames / 8);
lowSnrPxMask(~lowSnrPxMask) = nan;

% Combine the two masks
imageMask = satPxMask .* lowSnrPxMask;

% Visualize the same sample speckle frame shown above, but masked
fig2 = figure(currFig + 2); 
imagescnan(speckleMatrix(:, :, speckleFrameIndex) .* imageMask, speckleFrameColorbarLim)
colorbar, colormap(gray(256)), axis image, title('sample speckle frame - masked')
if savefigures
  filename = sprintf('sampleSpeckleFrame_masked');
  saveas(fig2, fullfile(outputPath, strcat(filename, '.png')));
end

%% Create spatial ensemble to compute the 2nd order ACF g2

% Create spatial Gaussian ensemble window
ensembleWindow = NormAnisotropicGaussianExp2Diam(...
  [xEnsembleWindowHalfSize, yEnsembleWindowHalfSize],...
  [xEnsembleWindowExp2Diameter, yEnsembleWindowExp2Diameter]);

%% Calculate g2 and K, and normalize g2 by modulation

fprintf('\t Calculating g2(tau), g2ModNorm(tau), BFI and ASDI ... ')

% Transform matrices so indices are: 1=y, 2=repetitions, 3=x, 4=time
speckleMatrixG2 = reshape(speckleMatrix, nY, nX, nFrames, []); 
speckleMatrixG2 = permute(speckleMatrixG2, [1 3 2 4]);

% Calculate the 2nd order ACF g2
g2 = CalculateG2(speckleMatrixG2, tauVec, ensembleWindow);
% Calculate contrast from g2
speckleContrastG2 = CalculateContrastFromG2(g2);
% Normalize g2 by modulation and get the contrast K_g2
g2ModNorm = NormalizeByModulation(g2);

% Reorder the indices 
g2 = permute(g2, [1 3 2]);
speckleContrastG2 = permute(real(speckleContrastG2), [1 3 2]);
g2ModNorm = permute(g2ModNorm, [1 3 2]);

% Calculate the blood flow index (BFI): beta_BFI
beta_BFI = 1 ./ (speckleContrastG2.^2); % BFI

% Plot the blood flow index (BFI): beta_BFI
fig3 = figure(currFig + 3); imagescnan(beta_BFI .* imageMask, betaBFIColorbarLim), colorbar
  colormap(colorMapBFI), axis image, title('\beta_{BFI}')
if savefigures
  filename = sprintf('BFI');    % file name to save BFI image
  saveas(fig3, fullfile(outputPath, strcat(filename, '.png')));
end  

%% Calculate the Adaptive Speckle Decorrelation Index (ASDI) alpha_ASDI [1]

% Calculate weight factor w [1]
w_ASDI = (g2ModNorm(:, : ,1) - g2ModNorm(:, :, 2)) ./ (g2ModNorm(:, :, 1) - g2ModNorm(:, :, end));
% Calculate ASDI metric [1]
alpha_ASDI = (g2ModNorm(:, :, 1) - g2ModNorm(:, :, 2)) .* w_ASDI;

% Plot ASDI
fig6 = figure(currFig + 6); imagescnan(alpha_ASDI .* imageMask, alphaASDIColorbarLim), colorbar 
  colormap(colorMapFASDA), axis image, title('\alpha_{ASDI}')
if savefigures
  filename = sprintf('ASDI'); % file name to save ASDI image
  saveas(fig6, fullfile(outputPath, strcat(filename, '.png')));
end

fprintf('done. \n')


