function [ average ] = RunningArbitraryAndLateralAve(array, ensembleWindow, dim)
  % RunningArbitraryAndLateralAve calculates average of array
  %   Calculates the average of array along dimension dim using ensembleWindow
  %
  % Inputs:
  %   array:            matrix with the signal
  %   ensembleWindow:   window used for averaging
  %   dim:              dimension along which the mean is calculated
  %
  % Outputs:
  %   average:          average of array
  %
  % This script and its functions follow the coding style that can be
  % sumarized in:
  % * Variables have lower camel case
  % * Functions upper camel case
  % * Constants all upper case
  % * Spaces around operators
  % Authors:  Néstor Uribe-Patarroyo
  %
  % NUP: 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA;
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  %
  % MGH Full-field amplitude speckle decorrelation angiography (FASDA) project
  %
  % Copyright Néstor Uribe-Patarroyo (2021)
  
  
  % Use imfilter for coherent 2D averaging along 1st and all other spatial dimensions (indices >= 4)
  sumAx = imfilter(array, ensembleWindow, 'replicate');
  % Coherent averaging along dim (can be a vector with multiple dims)
  sumLat = sum(sumAx, dim);
  % Count number of elements to do averaging. prod is required in case dim is a vector
  average = sumLat / (prod(size(array, dim)) * sum(ensembleWindow(:)));
end

