function kernelXY = NormAnisotropicGaussianExp2Diam(halfSize, exp2Diams)
  % NormAnisotropicGaussianExp2Diam generates a normalized Gaussian window
  % with input params
  %   Generates an anisotropic gaussian window of diamension diameters,
  %   and with diameters exp2DiamX and exp2Diams(2)
  %
  % Inputs:
  %   halfSize:     along [x, y] axes, final dimensions are [x, y] * 2 + 1
  %   exp2Diams:    exp2 diameters along [x-, y-] axes
  %
  % Outputs:
  %   h:      Gaussian window
  %
  % This function follows the coding style that can be sumarized in:
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
  % Changelog:
  %
  % V1.0 (2024-05-21): Initial version released
  
  [xMat, yMat] = meshgrid(-halfSize(1):halfSize(1), -halfSize(2):halfSize(2));
  % We assume a diameter of 0 really means an infinite diameter
  if (exp2Diams(1) ~= 0) && (exp2Diams(2) ~= 0)
    % 2D
    expArgument = -(8 .* xMat .^ 2 / (exp2Diams(1) ^ 2) + 8 .* yMat .^ 2 / (exp2Diams(2) ^ 2));
    kernelXY = exp(expArgument);
  elseif (exp2Diams(1) == 0) && (exp2Diams(2) ~= 0)
    % Along Y
    kernelX = xMat == 0;
    expArgument = -(8 .* yMat .^ 2 / (exp2Diams(2) ^ 2));
    kernelY = exp(expArgument);
    kernelXY = kernelY .* kernelX;
  elseif (exp2Diams(1) ~= 0) && (exp2Diams(2) == 0)
    % Along X
    kernelY = yMat == 0;
    expArgument = -(8 .* xMat .^ 2 / (exp2Diams(1) ^ 2));
    kernelX = exp(expArgument);
    kernelXY = kernelX .* kernelY;
  else
    kernelXY = 1;
  end
  % Make small numbers really zero
  kernelXY(kernelXY < eps * max(kernelXY(:))) = 0;
  % Normalize
  sumKernel = sum(kernelXY(:));
  if sumKernel ~= 0
    kernelXY  = kernelXY / sumKernel;
  end
end
