function SWC = myswt2(x, n, lo_D, hi_D)
%SWT2 Discrete stationary wavelet transform 2-D.
%   SWC = SWT2(X,N,Lo_D,Hi_D)
%   compute the stationary wavelet decomposition as above,
%   given these filters as input: 
%     Lo_D is the decomposition low-pass filter and
%     Hi_D is the decomposition high-pass filter.
%     Lo_D and Hi_D must be the same length.
%   SWC = [H(:,:,1:N) ; V(:,:,1:N); D(:,:,1:N); A(:,:,N)].
%
%   Filters lo_D and hi_D must be normalized by dividing by sqrt(2) (necessary for implmentation of Stationary Discrete WT)
%
%   See also DWT2, WAVEDEC2.

%   Sathish Ramani, October 11, 2009

% Preserve initial size.
s = size(x);

% Set DWT_Mode to 'per'.
old_modeDWT = dwtmode('status','nodisp');
modeDWT = 'per';
dwtmode(modeDWT,'nodisp');

% Compute stationary wavelet coefficients.
evenoddVal = 0;
evenLEN    = 1;
a = zeros(s(1),s(2),n);
h = zeros(s(1),s(2),n);
v = zeros(s(1),s(2),n);
d = zeros(s(1),s(2),n);

for k = 1:n

    % Extension.
    lf = length(lo_D);
    x  = wextend('2D',modeDWT,x,[lf/2,lf/2]);

    % Decomposition.
    y = wconv2('row',x,lo_D);    
    a(:,:,k) = wkeep2(wconv2('col',y,lo_D),s,[lf+1,lf+1]);
    h(:,:,k) = wkeep2(wconv2('col',y,hi_D),s,[lf+1,lf+1]);
    y = wconv2('row',x,hi_D);
    v(:,:,k) = wkeep2(wconv2('col',y,lo_D),s,[lf+1,lf+1]);
    d(:,:,k) = wkeep2(wconv2('col',y,hi_D),s,[lf+1,lf+1]);

    % upsample filters.
    lo_D = dyadup(lo_D,evenoddVal,evenLEN);
    hi_D = dyadup(hi_D,evenoddVal,evenLEN);

    % New value of x.
    x = a(:,:,k);

end
SWC = cat(3,h,v,d,a(:,:,n));

% Restore DWT_Mode.
dwtmode(old_modeDWT,'nodisp');
