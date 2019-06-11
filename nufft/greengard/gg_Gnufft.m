function ob = gg_Gnufft(om,M,Msp,R)
%function ob = gg_Gnufft(om, M, Msp, R)
%
% Object constructor for Leslie Greengard's 2-dimensional NUFFT using separable
% gaussian interpolator.
%
% in
%	om [Min,2]	"digital" frequencies in radians
%	M [1]		image dimensions (N1,N2,...,Nd)
%	Msp [1]		# of neighbors used
%	R [1]		oversampling ratio
% out
%	ob              the NUFFT object
%
% A type 1 NUFFT is evaluated using x = G' * b
% A type 2 NUFFT is evaluated using b = G * x
%
% Copyright 2008-6-6	Will Grissom and Jeff Fessler	The University of Michigan


% use greengard-supplied tau
tau = 1/M/M*pi/R/(R-0.5)*Msp; % width of gaussian

% initialize nufft
st = gg_Gnufft_init(om,M,Msp,R,tau);

% build an object with it.
ob = Fatrix([size(om,1) M*M], st, 'forw', @gg_Gnufft_forw, 'back', @gg_Gnufft_back, 'caller', mfilename);