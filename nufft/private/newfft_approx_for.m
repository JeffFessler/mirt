 function X = newfft_approx_for(st, x, om)
%function X = newfft_approx_for(st, x, om)
%|
%| (approximate) forward NUFFT

Nd = st.Nd;
Kd = st.Kd;

dims = size(x);
dd = length(Nd);
if ndims(x) < dd, fail 'input signal has too few dimensions', end
if any(dims(1:dd) ~= Nd), fail 'input signal has wrong size', end


% the usual case is where L=1, i.e., there is just one input signal.
if ndims1(x) == dd
	x = x .* st.sn; % apply scaling factors
	Xk = col(fftn_fast(x, Kd)); % [*Kd] oversampled FFT, padded at end

% otherwise, collapse all excess dimensions into just one
else
	xx = reshape(x, [prod(Nd) prod(dims((dd+1):end))]); % [*Nd *L]
	L = size(xx, 2);
	Xk = zeros(prod(Kd),L); % [*Kd *L]
	for ll=1:L
		xl = reshape(xx(:,ll), [Nd 1]); % l'th signal
		xl = xl .* st.sn; % scaling factors
		Xk(:,ll) = col(fftn_fast(xl, [Kd 1]));
	end
end

if ~isempty(st.phase_before)
	Xk = Xk .* repmat(st.phase_before(:), ncol(Xk));
end

% interpolate using precomputed sparse matrix or table
if isfield(st, 'interp_for') % table or other functional method
	if ~isvar('om') || isempty(om)
		om = st.om;
	end
	if isempty(om), fail 'om or st.om required', end

	X = st.interp_for(Xk, om);
else
	if isvar('om') && ~isempty(om) && ~isequal(om, st.om)
		fail 'om given that does not match st.om'
	end
	X = st.p * Xk; % [M *L]
end

if ndims1(x) > dd
	M = size(om,1);
	X = reshape(X, [M dims((dd+1):end)]); % [M (L)]
end

if isfield(st, 'phase_after') && ~isempty(st.phase_after)
	if isnumeric(st.phase_after)
		X = X .* repmat(st.phase_after, ncol(X));
	else
		phase = st.phase_after(om);
		X = X .* repmat(phase, ncol(X));
	end
end


% ndims1()
% variant of ndims() that returns "1" if x has size [N 1]
function nd = ndims1(x)
nd = ndims(x);
if nd == 2 && size(x,2) == 1
	nd = 1;
end
