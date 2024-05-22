function [g2Norm] = NormalizeByModulation(g2)
  % NormalizeG2ByModulation normalizes g2 by modulation
  %
  % Inputs:
  %   g2:       matrix with the second-order ACF
  %
  % Outputs:
  %   g2Norm:   matrix with the second-order ACF g2 normalizede by modulation
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
  % Copyright Néstor Uribe-Patarroyo (2021)
  
    nDims = ndims(g2);
    colonOp = repmat({':'}, [1 nDims - 2]);
    
    idealContrast = 1; % Defined as std/mean
    g2Norm = (g2 - 1) ./ (g2(:, 1, colonOp{:}) - 1) * (idealContrast .^ 2) + 1;
end

