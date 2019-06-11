 function st = gg_Gnufft_init(om, M, Msp, R, tau)
%function st = gg_Gnufft_init(om, M, Msp, R, tau)
%
% Initialize structure for 2-dimensional NUFFT using gaussian interpolator,
% particularly the interpolation matrix in sparse format.
% caution: this routine can require a lot of memory!
% in
%	om [Min,2]	"digital" frequencies in radians
%	M [1]		image dimensions (N1,N2,...,Nd)
%	Msp [1]		# of neighbors used
%	R [1]		oversampling ratio
%       tau [1]         Width of Guassian interpolator
% out
%	st.E1		[*M]
%       st.E2x          [*M]
%       st.E2y          [*M]
%       st.E3           [Msp]
%       st.E4           [*M]   Deconvolution kernel
%	st.Nd,Jd,Kd,om	copies of inputs
%
% *Nd is shorthand for prod(Nd).
% (Nd) is shorthand for (N1,N2,...,Nd)
%
% Like fft(), the NUFFT expects the signals to be x(0,0)?
%
% Copyright 2008-6-6	Will Grissom and Jeff Fessler	The University of Michigan

if nargin ~= 5, help(mfilename), error args, end

st.tau	= tau;

st.M	= M;
st.Msp	= Msp;
st.om	= om;
st.R    = R;
st.Mr   = R*M;

% precomputations
l = [-Msp+1:Msp]';

st.E3 = circshift(exp(-((pi*l/st.Mr).^2)/st.tau),Msp+1); % part of gaussian interpolation
                                                         % kernel
st.k = [-M/2:M/2-1]';
[st.kdc1,st.kdc2] = ndgrid(st.k);
st.E4 = exp(st.tau*(st.kdc1(:).^2+st.kdc2(:).^2)); % deconvolution kernel

% find the nearest grid points
st.m1 = floor(om(:,1)/2/pi*st.Mr);
st.m2 = floor(om(:,2)/2/pi*st.Mr);
st.xi1 = 2*pi*st.m1/st.Mr;
st.xi2 = 2*pi*st.m2/st.Mr;

% calculate irregular-point-dependent terms
st.E1 = exp(-((om(:,1)-st.xi1).^2 + (om(:,2)-st.xi2).^2)/4/st.tau);
st.E2x = exp(pi*(om(:,1)-st.xi1)./st.Mr./st.tau);
st.E2y = exp(pi*(om(:,2)-st.xi2)./st.Mr./st.tau);

