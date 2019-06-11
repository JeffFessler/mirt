  function window = fbp2_window(n, window)
%|function window = fbp2_window(n, window)
%| compute an apodizing window of length n and fft shift it

if nargin == 1 && streq(n, 'test'), fbp2_window_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

if ischar(window)
	if isempty(window) || streq(window, 'boxcar') || streq(window, 'ramp')
		window = ones(n,1);

	elseif streq(window, 'boxcar,', 7)
		cut = sscanf(window, 'boxcar,%g');
		window = my_boxcar(n, cut);

	elseif streq(window, 'hamming,', 8)
		cut = sscanf(window, 'hamming,%g');
		window = my_hamming(n, cut);

	elseif streq(window, 'hanning,', 8)
		cut = sscanf(window, 'hanning,%g');
		window = my_hann(n, cut);

	elseif streq(window, 'hann')
		window = my_hann(n, 1.0);
	elseif streq(window, 'hann50')
		window = my_hann(n, 0.5);
	elseif streq(window, 'hann75')
		window = my_hann(n, 0.75);
	elseif streq(window, 'hann80')
		window = my_hann(n, 0.80);

	else
		fail('unknown window %s', window)
	end

elseif length(window) ~= n
	error 'bad window length'
end

window = fftshift(window);


% my_boxcar()
function window = my_boxcar(n, cutoff)
w = round(cutoff * n);
ii = [0:n-1]'-n/2;
window = (abs(ii) < w/2);

% my_hann()
function window = my_hann(n, cutoff)
w = round(cutoff * n);
ii = [0:n-1]'-n/2;
window = 0.5 * (1 + cos(2*pi*ii/w)) .* (abs(ii) < w/2);


% my_hamming()
function window = my_hamming(n, cutoff)
w = round(cutoff * n);
ii = [0:n-1]'-n/2;
window = (0.54 + 0.46 * cos(2*pi*ii/w)) .* (abs(ii) < w/2);


% fbp2_window_test
function fbp2_window_test
n = 128;
types = {'boxcar', 'boxcar,0.6', 'hann', 'hamming,0.6'};
win = zeros(n, numel(types));
for it=1:numel(types)
	win(:,it) = fbp2_window(n, types{it});
end
plot(0:n-1, win, '-o'), xtick([0 n/2 n-1])
legend(types)
