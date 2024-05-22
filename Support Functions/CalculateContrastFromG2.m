function [Contrast_g2] = CalculateContrastFromG2(g2)
%CalcContrast2 Calculates the contrast K from the 2nd order ACF
  %
  % Inputs:
  %   g2:            4-D matrix containing the 2nd order ACF
  %   
  % Outputs:
  %   Contrast_g2:   contrast K
  %   
  % This script and its functions follow the coding style that can be
  % sumarized in:
  % * Variables have lower camel case
  % * Functions upper camel case
  % * Constants all upper case
  % * Spaces around operators
  %
  % Authors:  Néstor Uribe-Patarroyo, Giulia Mansutti
  %
  % NUP:
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA;
  % <uribepatarroyo.nestor@mgh.harvard.edu>
  %
  % GM:
  % 1. Wellman Center for Photomedicine, Harvard Medical School, Massachusetts
  % General Hospital, 40 Blossom Street, Boston, MA, USA;
  % <gmansutti@mgh.harvard.edu>
  
    
    nDims = ndims(g2);
    colonOp = repmat({':'}, [1 nDims - 2]);
    
    Contrast_g2 = sqrt(g2(:, 1, colonOp{:}) - 1);
end