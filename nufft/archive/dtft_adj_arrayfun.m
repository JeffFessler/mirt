 function x = dtft_adj_arrayfun(X, omega, Nd, n_shift, usearrayfun)
%function x = dtft_adj(X, omega, Nd, n_shift, usecellfun)
%|
%| Compute adjoint of d-dim DTFT for spectrum X at frequency locations omega
%|
%| in
%|	X	[M L]		dD DTFT values
%|	omega	[M d]		frequency locations (radians)
%|	n_shift [d 1]		use [0:N-1]-n_shift (default [0 ... 0])
%|	usearrayfun			1 to reduce memory use (slower)
%| out
%|	x	[(Nd) L]	signal values
%|
%| Requires enough memory to store M * (*Nd) size matrices. (For testing only.)
%|
%| Copyright 2003-4-13, Jeff Fessler, University of Michigan
%| Revised 2013-3-22, Daniel Weller, University of Michigan

% if no arguments, then run a simple test
if nargin < 2
	help(mfilename)
	Nd = [4 6 5]; Nd = [50 20 10];
	n_shift = [2 1 3]; 
	n_shift = 0*[2 1 3]; n_shift = [20 10 4];
	% test with uniform frequency locations:
	o1 = 2*pi*[0:(Nd(1)-1)]'/Nd(1);
	o2 = 2*pi*[0:(Nd(2)-1)]'/Nd(2);
	o3 = 2*pi*[0:(Nd(3)-1)]'/Nd(3);
	[o1 o2 o3] = ndgrid(o1, o2, o3);
	X = o1 + o2 - o3; % test spectrum
	om = [o1(:) o2(:) o3(:)];
	td = tic(); xd = dtft_adj(X(:), om, Nd, n_shift); td = toc(td);
	tl = tic(); xl = dtft_adj(X(:), om, Nd, n_shift, 1); tl = toc(tl);
	printm('loop max %% difference = %g (te = %g/%g)', max_percent_diff(xl,xd), tl, td)
	tc = tic(); xc = dtft_adj_arrayfun(X(:), om, Nd, n_shift, 1); tc = toc(tc);
	printm('arrayfun max %% difference = %g (te = %g/%g)', max_percent_diff(xc,xd), tc, td)
	Xp = X .* reshape(exp(-1i * om * n_shift(:)), size(X));
	xf = ifftn(Xp) * prod(Nd);
	printm('ifft max %% difference = %g', max_percent_diff(xf,xd))
return
end

if ~isvar('n_shift') || isempty(n_shift), n_shift = zeros(size(Nd)); end
if ~isvar('usearrayfun') || isempty(usearrayfun), usearrayfun = 0; end

% if length(Nd) == 1
% 	nn{1} = [0:(Nd(1)-1)] - n_shift(1);
% elseif length(Nd) == 2
% 	nn{1} = [0:(Nd(1)-1)] - n_shift(1);
% 	nn{2} = [0:(Nd(2)-1)] - n_shift(2);
% 	[nn{1} nn{2}] = ndgrid(nn{1}, nn{2});
% elseif length(Nd) == 3
% 	nn{1} = [0:(Nd(1)-1)] - n_shift(1);
% 	nn{2} = [0:(Nd(2)-1)] - n_shift(2);
% 	nn{3} = [0:(Nd(3)-1)] - n_shift(3);
% 	[nn{1} nn{2} nn{3}] = ndgrid(nn{1}, nn{2}, nn{3});
% else
% 	'only 1D-3D done'
% end
dd = size(omega,2);
Nd(end+1:dd) = 1;
n_shift(end+1:dd) = 0;

nn = cell(1,dd);
for id=1:dd % fixed: dd
	nn{id} = (0:(Nd(id)-1))-n_shift(id);
end
[nn{:}] = ndgrid(nn{:}); % fixed: dd
nn = cellfun(@(x) col(x),nn,'UniformOutput',false);
nn = cat(2,nn{:}); % [*Nd dd]
omega = omega.'; % [dd M]

if usearrayfun
    x = arrayfun(@(n) exp(1i.*(nn(n,:)*omega)) * X,1:prod(Nd),'UniformOutput',false);
    x = cat(1,x{:});
else
    x = exp(1i.*(nn*omega)) * X;
end
x = reshape(x, [Nd ncol(X)]); % [Nd L]
