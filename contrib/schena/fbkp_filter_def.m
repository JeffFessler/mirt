function filt =fbkp_filter_def(filter, len, d)
%%%
%%%  Sub-Function:  designFilter
%%%
% Returns the Fourier Transform of the filter which will be 
% used to filter the projections
%
% INPUT ARGS:   filter - either the string specifying the filter 
%               len    - the length of the projections
%               d      - the fraction of frequencies below the nyquist
%                        which we want to pass
%
% OUTPUT ARGS:  filt   - the filter to use on the projections


order = max(64,2^nextpow2(2*len));

% First create a ramp filter - go up to the next highest
% power of 2.

filt = 2*( 0:(order/2) )./order;
w = 2*pi*(0:size(filt,2)-1)/order;   % frequency axis up to Nyquist 

switch filter
case 'ram-lak'
   % Do nothing
case 'shepp-logan'
   % be careful not to divide by 0:
   filt(2:end) = filt(2:end) .* (sin(w(2:end)/(2*d))./(w(2:end)/(2*d)));
case 'cosine'
   filt(2:end) = filt(2:end) .* cos(w(2:end)/(2*d));
case 'hamming'  
   filt(2:end) = filt(2:end) .* (.54 + .46 * cos(w(2:end)/d));
case 'hann'
   filt(2:end) = filt(2:end) .*(1+cos(w(2:end)./d)) / 2;
otherwise
   error('Invalid filter selected.');
end

filt(w>pi*d) = 0;                      % Crop the frequency response
filt = [filt' ; filt(end-1:-1:2)'];    % Symmetry of the filter
