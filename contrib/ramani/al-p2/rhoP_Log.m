function [drho, tE] = rhoP_Log(x, params)
%
% Evaluate the derivative of the Logarithmic 'rho' in the paper 
% "Highly Undersampled MRI Reconstruction" by Tryasko et al, IEEE TMI
%
% x      ->      Point at which rho needs to be evaluated
% sigma  ->      Parameter that decides the concavity of rho
%
% Normalize rho such that r(1,sigma) = 1

sigma = params.Prior.sigma;
den = params.Prior.den;

tS = tic;
drho = 1./((x + sigma)*den);
tE = toc(tS);