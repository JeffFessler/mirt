function swc = myswt2(x, n, lo_d, hi_d)
%SWT2 Discrete stationary wavelet transform 2-D.
% swc = SWT2(X,N,Lo_d,Hi_d)
% compute the stationary wavelet decomposition as above,
% given these filters as input:
% lo_d is the decomposition low-pass filter and
% hi_d is the decomposition high-pass filter.
% lo_d and hi_d must be the same length.
% swc = [H(:,:,1:N) ; V(:,:,1:N); D(:,:,1:N); A(:,:,N)].
%
%  Filters lo_d and hi_d must be normalized by dividing by sqrt(2)
% (necessary for implmentation of Stationary Discrete WT)
%
% See also DWT2, WAVEDEC2.

% based on matlab's swt2.m
% Sathish Ramani, October 11, 2009

% Preserve initial size.
s = size(x);

% Set DWT_Mode to 'per'.
old_modeDWT = dwtmode('status','nodisp');
modeDWT = 'per';
dwtmode(modeDWT,'nodisp');

% Compute stationary wavelet coefficients.
evenoddVal = 0;
evenLEN	= 1;
a = zeros(s(1),s(2),n);
h = zeros(s(1),s(2),n);
v = zeros(s(1),s(2),n);
d = zeros(s(1),s(2),n);

for k = 1:n

	% Extension
	lf = length(lo_d);
	x  = wextend('2D', modeDWT, x, [lf/2,lf/2]);

	% Decomposition
	y = wconv2('row', x, lo_d);
	a(:,:,k) = wkeep2(wconv2('col',y,lo_d), s, [lf+1,lf+1]);
	h(:,:,k) = wkeep2(wconv2('col',y,hi_d), s, [lf+1,lf+1]);
	y = wconv2('row',x,hi_d);
	v(:,:,k) = wkeep2(wconv2('col',y,lo_d), s, [lf+1,lf+1]);
	d(:,:,k) = wkeep2(wconv2('col',y,hi_d), s, [lf+1,lf+1]);

	% upsample filters
	lo_d = dyadup(lo_d, evenoddVal, evenLEN);
	hi_d = dyadup(hi_d, evenoddVal, evenLEN);

	% New value of x
	x = a(:,:,k);

end
swc = cat(3,h,v,d,a(:,:,n));

% Restore DWT_Mode.
dwtmode(old_modeDWT,'nodisp');
