 function ws = unwrapping_sps_manysets(w, y, delta, R, varargin)
%function ws = unwrapping_sps_manysets(w, y, delta, R, [options])
%|
%| Phase unwrapping of multiple data sets (w denotes frequency)
%| using a separable quadratic surrogates algorithm.
%| can be 1d, 2d, etc., depending on R.
%|
%| cost(w) = sum_{i=1 to n} sum_{j=1 to n}
%|		|yi*yj| (1 - cos(w * (dj-di) + \angle(yi) - \angle(yj))) + R(w)
%|
%| in
%|	w	[np 1]		initial estimate
%|	y	[np M]		M sets of measurements
%|	delta	[M]		vector of M offsets
%|				note - the first delta is normally just 0.
%|	R			penalty object (see Reg1.m)
%|
%| options
%|	niter			# of iterations (default 1)
%|	chat
%|
%| out
%|	ws	[np niter]	iterates
%|
%| Use "mri_field_map_reg test" to see how it works.
%|
%| Copyright 2007-12-15, Amanda Funai, University of Michigan

if nargin == 1 && streq(w, 'test')
	run_mfile_local('mri_field_map_reg test')
return
end
if nargin < 3, help(mfilename), error(mfilename), end

arg.niter = 1;
arg.chat = false;
arg = vararg_pair(arg, varargin);

[np M] = size(y);

if M ~= numel(delta), fail 'need delta to be [M]', end

ang = angle(y);

% calculate the magnitude here to avoid recomputing each iteration
mag = zeros(np, M*M);
wj = zeros(np, M*M);

% todo: there is a bug here compared to the paper!

set = 1;
for i=1:M
	for j=1:M
		wj(:,set) = abs(y(:,i)).^2 .* abs(y(:,j)).^2;
		mag(:,set) = abs( conj(y(:,i)) .* y(:,j) );
		set = set + 1;
	end
end
wjtotal = sum(abs(y).^2, 2);

%% todo: comment what is going on here!
for i=1:size(wj,2)
	wj(:,i) = div0(wj(:,i), wjtotal .* abs(y(:,1)).^2);
end


% loop over iterations
ws = zeros(length(w(:)), arg.niter);
ws(:,1) = w;

for iter = 2:arg.niter
	if arg.chat, printm('unwrap iteration %d', iter-1), end

	% num & denom contribution for each surrogate function
	grad = 0;
	denom = 0;
	set = 1;
	for i=1:M
		for j=1:M
			s = w .* (delta(j) - delta(i)) + ang(:,i) - ang(:,j);
%			if any(isnan(s)), 's', keyboard, end
			grad = grad + wj(:,set) .* mag(:,set) ...
				.* (delta(j) - delta(i)) .* sin(s);
%			if any(isnan(grad)), 'grad', keyboard, end
			sr = mod(s + pi, 2*pi) - pi;
%			if any(isnan(sr)), 'sr', keyboard, end
			denom = denom + wj(:,set) .* mag(:,set) ...
				.* (delta(j) - delta(i))^2 .* nufft_sinc(sr / pi);
%			if any(isnan(denom)), 'denom', keyboard, end
			set = set + 1;
		end
	end

	if ~isempty(R)
		num = grad + R.cgrad(R, w);
		den = denom + R.denom(R, w);
	end

	w = w - div0(num, den);

	if arg.chat, printm('Range %g %g', min(w), max(w)), end

	ws(:,iter) = w;

	if any(isnan(ws(:))), warn 'bug: nan', keyboard, end
end
