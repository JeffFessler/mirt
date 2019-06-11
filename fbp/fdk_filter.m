 function proj = fdk_filter(proj, window, dsd, dfs, ds)
%function proj = fdk_filter(proj, window, dsd, dfs, ds)
%|
%| step 2 of FDK cone-beam CT reconstruction:
%| filter the (zero padded) projections
%|
%| in
%|	proj	[ns nt na]
%|	window	[npad]		or 'ramp' 'hann'
%|				or see fbp2_window()
%|
%| out
%|	proj	[ns nt na]
%|
%| todo: should allow view-dependent ds for mojette
%|
%| Jeff Fessler and others, University of Michigan

if nargin < 5, help(mfilename), error(mfilename), end

[ns nt na] = size(proj);
npadh = 2^ceil(log2(2*ns-1));
%printm('ns=%d npadh=%d', ns, npadh)

if isinf(dsd) % parallel-beam
	H = fdk_fan_filter('flat', npadh, ds, [], window); % [npad 1]
elseif isinf(dfs) % flat
	H = fdk_fan_filter('flat', npadh, ds, [], window); % [npad 1]
elseif dfs == 0 % arc
	H = fdk_fan_filter('arc', npadh, ds, dsd, window);
end
H = ds * H; % differential for discrete-space convolution vs integral

proj = ifft_sym( fft(proj, npadh, 1) .* repmat(H, [1 nt na]) );
proj = proj(1:ns,:,:); % back to original unpadded size
end % fdk_filter()


% fdk_fan_filter()
% apodized filter frequency response
function H = fdk_fan_filter(type, n, ds, dsd, window)

if streq(type, 'flat')
	h = fbp_ramp('flat', n, ds);
else
	h = fbp_ramp('arc', n, ds, dsd);
end
H = reale(fft(fftshift(h)));

if ischar(window)
	if streq(window, 'ramp')
		window = ones(n,1);
	elseif streq(window, 'hann')
		window = ir_hann_periodic(n); % cf hann(n, 'periodic');
	else
		window = fftshift(fbp2_window(n, window));
	end
elseif length(window) ~= n
	error 'bad window length'
end

H = H .* fftshift(window);
end % fdk_fan_filter()
