 function X = dtft(x, omega, varargin)
%function X = dtft(x, omega [,options])
%|
%| Compute d-dimensional DTFT of signal x at frequency locations omega
%|
%| in
%|	x	[(Nd) L]	signal values
%|	omega	[M dd]		frequency locations (radians), dd = numel(Nd)
%|
%| option
%|	'n_shift' [dd 1]	use [0:N-1]-n_shift (default [0 .. 0])
%|	'how'			'outer' (default) big outer product
%|				'loop' reduce memory use (slower)
%|				'arrayfun' uses arrayfun()
%|
%| out
%|	X	[M L]		DTFT values
%|
%| Requires enough memory to store M * prod(Nd) size matrices (for testing)
%|
%| Copyright 2001-9-17, Jeff Fessler, University of Michigan
%| 2013-03-22, Daniel Weller added arrayfun version and other improvements
%| 2013-03-27, JF converted to vararg

if nargin == 1 && streq(x, 'test'), dtft_test(0), return, end
if nargin == 1 && streq(x, 'time'), dtft_test(1), return, end
if nargin < 2, ir_usage(), end

arg.n_shift = 0;
arg.how = 'outer';
arg = vararg_pair(arg, varargin);

dd = size(omega, 2);
Nd = size(x);

if numel(arg.n_shift) == 1
	arg.n_shift = repmat(arg.n_shift, dd);
end
n_shift = arg.n_shift;
if numel(n_shift) ~= dd
	fail 'n_shift size bad')
end

if dd == 1 && numel(Nd) == 2 && Nd(2) == 1 % 1D
	Nd = Nd(1);
end

if length(Nd) == dd % just one image
	x = x(:);
elseif length(Nd) == dd+1 % multiple images
	Nd = Nd(1:(end-1));
	x = reshapee(x, prod(Nd), []); % [*Nd L]
else
	error 'bad input signal size'
end

% dsw alternative to the loop:
% nn = arrayfun(@(nd,nshift) (0:(nd-1))-nshift, Nd, n_shift, 'UniformOutput', false);
for id=1:dd
	nn{id} = [0:(Nd(id)-1)] - n_shift(id);
end

% nn = ndgrid_jf('cell', nn);
if dd > 1
	[nn{:}] = ndgrid(nn{:});
end

switch arg.how
case 'outer'
	X = dtft_outer(x, omega, nn);
case 'loop'
	X = dtft_loop(x, omega, Nd, nn);
case 'arrayfun'
	X = dtft_arrayfun(x, omega, nn);
otherwise
	fail('unknown how "%s"', arg.how)
end


% dtft_outer()
function X = dtft_outer(x, omega, nn);
X = 0;
dd = ncol(omega);
for id=1:dd % add up phases
	X = X + omega(:,id) * col(nn{id})'; % [M *Nd]
end
X = exp(-1i*X) * x;


% dtft_loop()
% loop way: slower but less memory
function X = dtft_loop(x, omega, Nd, nn);
M = nrow(omega);
X = zeros(numel(x)/prod(Nd),M); % [L M]
%t1 = col(nn{1})';
if ncol(omega) > 3
	fail 'only up to 3d done'
end
if ncol(omega) < 3
	for dd = (ncol(omega)+1) : 3
		nn{dd} = 0; % dummy 0's
	end
	omega(1,3) = 0; % trick: make '3d'
end
t1 = nn{1}(:)';
t2 = col(nn{2})';
t3 = col(nn{3})';
for mm=1:M
	tmp = omega(mm,1)*t1 + omega(mm,2)*t2 + omega(mm,3)*t3;
	X(:,mm) = exp(-1i * tmp) * x;
end
X = X.'; % [M L]


% dtft_arrayfun()
% by Dan Weller, 2013-03-27
function X = dtft_arrayfun(x, omega, nn);
nn = cellfun(@(x) col(x).',nn,'UniformOutput',false); % make row vectors
nn = cat(1,nn{:}); % [dd *Nd]

M = nrow(omega);
X = arrayfun(@(m) exp((-1i*omega(m,:)) * nn) * x, 1:M, 'UniformOutput', false); % each cell is [1 L]
X = cat(1,X{:}); % [M L]

% X = exp(-1i*(omega * nn)) * x; % [M L]


% dtft_test
% simple test
function dtft_test(do_time)
Nd = [4 6 5] * 2^(1+do_time);
n_shift = [1 3 2];
rng(0), x = randn(Nd); % test signal
o1 = 2*pi*[0:(Nd(1)-1)]'/Nd(1); % test with uniform frequency locations
o2 = 2*pi*[0:(Nd(2)-1)]'/Nd(2);
o3 = 2*pi*[0:(Nd(3)-1)]'/Nd(3);
[o1 o2 o3] = ndgrid(o1, o2, o3);
om = [o1(:) o2(:) o3(:)];
cpu etic
Xd = dtft(x, om, 'n_shift', n_shift);
cpu etoc outer
cpu etic
Xl = dtft(x, om, 'n_shift', n_shift, 'how', 'loop');
cpu etoc loop
cpu etic
Xa = dtft(x, om, 'n_shift', n_shift, 'how', 'arrayfun');
cpu etoc arrayfun
printm('loop max %% difference = %g', max_percent_diff(Xl,Xd))
printm('afun max %% difference = %g', max_percent_diff(Xa,Xd))
Xf = fftn(x);
Xf = Xf(:) .* exp(1i * (om * n_shift(:))); % phase shift
printm('fftn max %% difference = %g', max_percent_diff(Xf,Xd))
