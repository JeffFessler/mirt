 function x = nufft_adj(X, st)
%function x = nufft_adj(X, st)
%| Apply adjoint of d-dimensional NUFFT to spectrum vector(s) X
%| in
%|	X	[M,(L)]
%|	st			structure precomputed by nufft_init()
%| out
%|	x	[(Nd),(L)]	signal(s)/image(s)
%|
%| Copyright 2003-6-1, Jeff Fessler, University of Michigan

% if no arguments, then run a simple test
if nargin < 1
	help(mfilename)

	N1 = 4;
	N2 = 8;
	n_shift = [2.7 3.1]; % random shifts to stress it
	o1 = 2 * pi * [0.0 0.1 0.3 0.4 0.7 0.9]';
	o2 = flipud(o1);
	om = [o1 o2];
	st = nufft_init(om, [N1 N2], [4 6], 2*[N1 N2], n_shift, 'minmax:tuned');
	X = [1:length(o1)].^2'; % test spectrum
	xd = dtft2_adj(X, om, N1, N2, n_shift);
	xn = nufft_adj([X X X], st);
	printm('nufft vs dtft max%%diff = %g', max_percent_diff(xd,xn(:,:,2)))
%	s2 = nufft2_init([o1 o2], N1, N2, 4, 6, 2*N1, 2*N2, n_shift, 0, 'best');
%	printm('nufft vs nufft2 max%%diff = %g', max_percent_diff(s2.p,st.p))
return
end

% extract attributes from structure
Nd = st.Nd;
Kd = st.Kd;
dims = size(X);
if dims(1) ~= st.M, fail 'size', end

% adjoint of interpolator using precomputed sparse matrix
if length(dims) > 2 || dims(2) > 1
	Lprod = prod(dims(2:end));
	X = reshape(X, [st.M Lprod]); % [M *L]
else
	Lprod = 1; % the usual case
end

if ~isfield(st, 'interp_table_adj')
	Xk_all = full(st.p' * X);		% [*Kd *L]
else
	Xk_all = feval(st.interp_table_adj, st, X);
end

for ll=1:Lprod
	Xk = reshape(Xk_all(:,ll), [Kd 1]);	% [(Kd)]
	tmp = prod(Kd) * col(ifftn_fast(Xk));	% scale factor!
	if ll == 1
		x = zeros(prod(Kd), Lprod, class(tmp)); % [*Kd *L]
	end
	x(:,ll) = tmp;
end
x = reshape(x, [Kd Lprod]);			% [(Kd) *L]

% eliminate zero padding from ends  fix: need generic method
switch length(Nd)
case 1
	x = x(1:Nd(1),:); % [N1 *L]
case 2
	x = x(1:Nd(1),1:Nd(2),:); % [N1 N2 *L]
case 3
	x = x(1:Nd(1),1:Nd(2),1:Nd(3),:); % [N1 N2 N3 *L]
case 4
	x = x(1:Nd(1),1:Nd(2),1:Nd(3),1:Nd(4),:); % [N1 N2 N3 N4 *L]
otherwise
	error 'only up to 4D implemented'
end

% scaling factors
x = reshape(x, [prod(Nd) Lprod]);		% [*Nd *L]
snc = conj(col(st.sn));				% [*Nd 1]
x = x .* snc(:,ones(1,Lprod));			% scaling factors
x = reshape(x, [Nd dims(2:end)]);		% [(Nd) (L)]
