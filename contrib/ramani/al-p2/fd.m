function [h, v, tE] = fd(x)
% Compute horizontal and vertical finite differences on the 
% current estimate x and then compute the absolute of the 
% finite-difference-gradient
% 
% Use periodic boundaries


%% Horizontal finite difference
tS = tic;
h = imfilter(x, [-1 1], 'circular');

%% Vertical finite difference
v = imfilter(x, [-1 1]', 'circular');
tE = toc(tS);
