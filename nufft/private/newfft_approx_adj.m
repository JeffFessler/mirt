 function x = newfft_approx_adj(st, X, om)
%function x = newfft_approx_adj(st, X, om)
%|
%| Adjoint of (approximate) NUFFT
%|
%| in
%|	st			structure precomputed by nufft_init()
%|	X	[M (L)]
%|	om	[M 1]		optional frequency vector (default: st.om)
%|
%| out
%|	x	[(Nd) (L)]	signal(s)/image(s)
%|
%| Copyright 2003-6-1, Jeff Fessler, University of Michigan

Nd = st.Nd;
Kd = st.Kd;

dims = size(X);
dd = length(Nd);

if length(dims) > 2 || dims(2) > 1
	Lprod = prod(dims(2:end));
	X = reshape(X, [st.M Lprod]);	% [M *L]
else
	Lprod = 1; % usual case
end


% adjoint via conjugate
if isfield(st, 'phase_after') && ~isempty(st.phase_after)
	if isnumeric(st.phase_after)
		X = X .* repmat(conj(st.phase_after), ncol(X));
	else
		phase = st.phase_after(om);
		X = X .* repmat(conj(phase), ncol(X));
	end
end

% adjoint of interpolator using precomputed sparse matrix or table
if isfield(st, 'interp_adj') % table or other functional method
	if ~isvar('om') || isempty(om)
		om = st.om;
	end
	if isempty(om), fail 'om or st.om required', end
	if length(om) ~= dims(1), fail 'om vs X size mismatch', end
	Xk_all = st.interp_adj(X, om);
else
	if isvar('om') && ~isempty(om) && ~isequal(om, st.om)
		fail 'om given that does not match st.om'
	end
	if length(st.om) ~= dims(1), fail 'om vs X size mismatch', end
	Xk_all = full(st.p' * X);		% [*Kd *L]
end

if ~isempty(st.phase_before) % adjoint via conjugate
	Xk_all = Xk_all .* repmat(conj(st.phase_before(:)), ncol(Xk_all));
end

x = zeros(prod(Kd), Lprod);			% [*Kd *L]
for ll=1:Lprod
	Xk = reshape(Xk_all(:,ll), [Kd 1]);		% [(Kd)]
	x(:,ll) = prod(Kd) * col(ifftn_fast(Xk));	% scale factor!
end
x = reshape(x, [Kd Lprod]);			% [(Kd) *L]

% eliminate zero padding from ends
if length(Nd) == 1
	x = x(1:Nd(1),:);			% [N1 *L]
elseif length(Nd) == 2
	x = x(1:Nd(1),1:Nd(2),:);		% [N1 N2 *L]
elseif length(Nd) == 3
	x = x(1:Nd(1),1:Nd(2),1:Nd(3),:);	% [N1 N2 N3 *L]
else % general
	arg = 'x(1:Nd(1)';
	for id=2:numel(Nd)
		arg = [arg sprintf(',1:Nd(%d)', id)];
	end
	arg = [arg ');'];
	x = eval(arg);
	warn 'fix: untested'
end

% scaling factors
x = reshape(x, [prod(Nd) Lprod]);		% [*Nd *L]
snc = conj(col(st.sn));				% [*Nd 1]
x = x .* repmat(snc, [1 Lprod]);		% scaling factors
x = reshape(x, [Nd dims(2:end)]);		% [(Nd) (L)]
