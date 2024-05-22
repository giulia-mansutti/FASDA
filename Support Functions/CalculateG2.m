function [g2] = CalculateG2(speckleMatrix, tauVec, ensembleWindow, varargin)
  % CalculateG2 Calculates the 2nd order ACF along 2nd index
  %   Calculates the 2nd order autocorrelation function (ACF) along 2nd index,
  %   with optional inputs.
  %
  % Inputs:
  %   speckleMatrix:    matrix with the signal
  %   tauVec:           vector with desired delays
  %   ensembleWindow:   kernel for ensemble averaging in arbitrary dimensions.
  %   addEnsembleDims:  vector with list of dimensions for further ensemble averaging (dims are collapsed)
  %
  % Outputs:
  %   g2:      2nd order ACF
  %   
  % This script and its functions follow the coding style that can be
  % sumarized in:
  % * Variables have lower camel case
  % * Functions upper camel case
  % * Constants all upper case
  % * Spaces around operators
  %
  % Authors:  Néstor Uribe-Patarroyo
  %
  % NUP: 
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA;
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  %
  % MGH Full-field amplitude speckle decorrelation angiography (FASDA) project
  %
  % Changelog:
  %
  % V1.0 (2017-09-01): Initial version released
  % V2.0 (2021-01-09): Probabilistic angiography version
  % V3.0 (2024-05-01): FASDA version
  %
  % Copyright Néstor Uribe-Patarroyo (2021)
  
  if nargin < 4 || ~isempty(varargin{1})
    addEnsembleDims = [];
  else
    addEnsembleDims = varargin{1};
  end
  
  % Replace all the ensemble averages to include axial running averaging too
  if numel(ensembleWindow) > 1
    % Shift all dims >= 2 by two due to the similar shift done to the signal below
    if ~iscolumn(ensembleWindow)
      ensembleWindow = shiftdim(ensembleWindow, -1);
      ensembleWindow = permute(ensembleWindow, [2, 1, 3:ndims(ensembleWindow)]);
    end
    meanEnsemble = @(x, dim) RunningArbitraryAndLateralAve(x, ensembleWindow, [dim, addEnsembleDims]);
  else
    meanEnsemble = @(x, dim) mean(x, [dim, addEnsembleDims]);
  end
  
  % Get info on speckle matrix containing the signal
  nDims = ndims(speckleMatrix);
  colonOp = repmat({':'}, [1 nDims - 2]);
  nY = size(speckleMatrix, 1);
  corrWindow = size(speckleMatrix, 2);
  nOther = num2cell(size(speckleMatrix));
  nOther(1:2) = [];
  % Make addEnsembleDims singleton
  addEnsembleDimsCell = num2cell(addEnsembleDims - 2);
  nOther([addEnsembleDimsCell{:}]) = deal({1});
  % Fix tauVec size
  tauVec = unique(min(tauVec, corrWindow - 1));
  nTaus = numel(tauVec);

  % Now we create the complex-conjugate products maxtrices
  g2 = zeros(nY, nTaus, nOther{:}, 'like', speckleMatrix);

  m = 0;
  for thisTau = tauVec
    m = m + 1;
    
    % Get signals' indexes
    idx1 = 1:corrWindow - thisTau;
    idx2 = 1 + thisTau:corrWindow;
    
    % Get signals
    i1 = speckleMatrix(:, idx1, colonOp{:});
    i2 = speckleMatrix(:, idx2, colonOp{:});
    
    % Calculate products with proper normalization.
    % (Ensemble averages are in the 2nd index)    
    intProdIMean = bsxfun(@rdivide, meanEnsemble(i1 .* i2, 2), meanEnsemble(i1, 2) .* meanEnsemble(i2, 2));

    % Get meanEnsemble of intProd to obtain correlation coefficient of ensemble
    g2(:, m, colonOp{:}) = intProdIMean;
  end
  
end

