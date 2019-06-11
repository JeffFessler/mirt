 function [C, ugibb] = Csparse(arg0, arg1, arg2, arg3)
%function C = Csparse('b2info', kappa, mask, chat)
%function C = Csparse('center', Cscale, mask, chat)
%function C = Csparse('maskleak', mask, chat)
%function C = Csparse('maskleak3', [beta_xy beta_z], mask, chat)
%function C = Csparse('maskbar', mask, chat)
%function C = Csparse('test')
%function Csparse('plot', C, mask)
% create 2D (or 1D or 3D) penalty matrix C
% R = beta * C'C is the penalty Hessian (in quadratic case)
% in
%	kappa	[nx,ny]		as in :srp paper, i.e.
%				sqrt( (A.^2)'*nder2(:) ./ sum(A.^2)' )
%	mask	[nx,ny[,nz]]	support
% out
%	C	[nc,np]		np = # of pixels in mask
%
% maskleak: penalize across mask boundary, like in ASPIRE (C*1 != 0)
%	(matches "-" penalty in ASPIRE)
% maskbar: barrier across mask boundary (C*1 = 0)
% b2info: as in 'b2info' penalty in ASPIRE
% center: as in 'center' penalty in ASPIRE

if nargin == 1 && streq(arg0, 'test'), Csparse_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

warning 'Csparse is obsolete.  C2sparse is recommended instead!'

betas = 1;

if strcmp(arg0, 'maskleak3')
	Cscale = 1;
	mask = arg2 ~= 0;
	betas = arg1;
	chat = 0;		if 1==exist('arg3'), chat = arg3; end

elseif strcmp(arg0, 'maskleak') || strcmp(arg0, 'maskbar')
	Cscale = 1;
	mask = arg1 ~= 0;
	chat = 0;		if 1==exist('arg2'), chat = arg2; end

elseif strcmp(arg0, 'center')
	Cscale = arg1;
	mask = arg2 ~= 0;
	chat = 0;		if 1==exist('arg3'), chat = arg3; end
	arg0 = 'maskleak';	% trick!

elseif strcmp(arg0, 'b2info')
	kappa = arg1;
	mask = kappa~=0;	if 1==exist('arg2'), mask = arg2; end
	chat = 0;		if 1==exist('arg3'), chat = arg3; end
	if ~islogical(mask), error 'mask must be logical', end
	if any(mask(:,1)) || any(mask(:,end)) || any(mask(1,:)) || any(mask(end,:))
		error 'mask must not include outer row/col'
	end

%
% test by calculating penalty for an image here and in aspire
%
elseif streq(arg0, 'test')
	nx = 20; ny = 12;	nb = nx+1; na = ny;
	mask = zeros(nx,ny);	mask(4:14,4:8) = ones(11,5);
	mask = mask ~= 0;

	A = [spdiag(mask(:), 'nowarn'); zeros(ny,nx*ny)];
	rng(0)
	yy = rand(nb,na);
	xx = mask .* rand(nx,ny);
	dd = 1 + 1 * rand(nb,na);
	kappa = mask+0;
	jb = mask(:);
	kappa(jb) = sqrt( (A(:,jb).^2)'*dd(:) ./ sum(A(:,jb).^2)' );

whos
error('hmm')
	delete('t,A.wtf')
	wtf_write('t,A.wtf', A, nx, ny, nb, na)
	save('t,d', 'dd')
	save('t,y', 'yy')
	save('t,x', 'xx')

	err = yy(:)-A*xx(:);
%	C = Csparse('maskleak', mask);		method = '@1@alg@0,quad,1,-'
	C = Csparse('b2info', kappa, mask, 1);	method = '@1@alg@0,quad,1,b2info'

	printm('like=%.11e penal=%.11e', ...
		err' * spdiag(dd, 'nowarn') * err / 2, norm(C * xx(jb)).^2 / 2)
	eval(['!i best,pwls2 - t,x.mat t,y.mat - t,d.mat t,A.wtf ' method])
	C = 0;
	return

elseif strcmp(arg0, 'plot')
	C = arg1;
	mask = arg2 ~= 0;
%	im(221, kappa, 'kappa')
	if im
		subplot(338), spy(C), title 'C'
	end
	if ncol(C) < 500
		im(339, C'*C, 'C''C')
	else
		im(339, embed(diag(C'*C), mask), 'Diag(C''C)'), cbar
	end
	return

else
	help Csparse
	error([arg0 ' unknown'])
end

	[nx,ny,nz] = size(mask);

%
% C for 1st-order differences
%
if 1
	dx	= sparse(filtmat('1d', [-1 1 0], nx));
	if ~strcmp(arg0, 'maskleak')
		dx(1,1) = 0;
	end
%	dx	= dx(2:end,:);
	dx	= kron(speye(ny*nz), dx);
	C	= dx;	clear dx
	if ny > 1
		dy	= sparse(filtmat('1d', [-1 1 0], ny));
		dy(1,1) = 0;
%		dy	= dy(2:end,:);
		if nz > 1
			dy	= kron(speye(nz), kron(dy, speye(nx)));
		else
			dy	= kron(dy, speye(nx));
		end
		b.xy = betas(1);
		C	= sqrt(b.xy) * [C; dy];	clear dy
	end
	if nz > 1
		dz	= sparse(filtmat('1d', [-1 1 0], nz));
%		dz	= dz(2:end,:);
		dz(1,1) = 0;
		dz	= kron(dz, speye(nx*ny));
		b.z = betas(2);
		C	= [C; sqrt(b.z)*dz];	clear dz
	end
end


if strcmp(arg0, 'maskleak') || strcmp(arg0, 'maskleak3')
	C = Cscale * C(:, mask(:));
	if chat
		Csparse('plot', C, mask);
	end
	return
end

% only rows (and cols) of C where all elements are in mask
if strcmp(arg0, 'maskbar')
	t = abs(C) * mask(:);
	t = full(t == sum(abs(C)')');
	C = C(t,:);
	C = C(:,mask(:));
	return
end


if ~strcmp(arg0, 'b2info'),	error bug, end

	if any(kappa(~mask(:))), error('kappa outside of mask'), end

	%
	% 'fix edges' of kappa, as in ASPIRE
	%
	offset = [-1 1 -nx nx nx-1 nx+1 1-nx -1-nx];
	for jj = find(mask(:))'
		bad = offset(mask(jj + offset) == 0);
		if bad % ~= []
			kappa(jj + bad) = kappa(jj) * ones(size(bad));
		end
	end

	t = abs(C) * spdiag(1:nx*ny, 'nowarn');
	t2 = max(t')';
	t1 = sum(t')' - t2;	% trick works with 1st order only!
	if any(t1 == t2 | t1 <= 0 | t2 <= 0)
		if sum(t1 == 0) ~= nx + ny, error bug_t1, end
		if sum(t2 == 0) ~= nx + ny, error bug_t2, end
		if any(t1(t1 ~= 0) == t2(t2 ~= 0)), error bug_t12, end
		tdiag = zeros(nx*ny,1);
		tdiag(t1 ~= 0) = sqrt(kappa(t1(t1 ~= 0)) .* kappa(t2(t1 ~= 0)));
		% all the above because now we keep zero rows in C...
	else
		tdiag = sqrt(kappa(t1) .* kappa(t2));
	end
	C = spdiag(tdiag, 'nowarn') * C;

	C = C(:,mask(:));

	if chat
		Csparse('plot', C, mask);
	end

	if nargout > 1
		ugibb.x = zeros(nx,ny);
		ugibb.y = zeros(nx,ny);
		ugibb.p = zeros(nx,ny);
		ugibb.n = zeros(nx,ny);

		t = (t2 == t1 + 1) & (tdiag ~= 0);
		sum(t)
		ugibb.x(t2(t)) = tdiag(t);

		t = (t2 == t1 + nx) & (tdiag ~= 0);
		sum(t)
		ugibb.y(t2(t)) = tdiag(t);
		if chat
			im([ugibb.x ugibb.y])
		end
	end

%
% self test
%
function Csparse_test
% run test with 'kappa' for b2info
x = shepplogan(64,64,1); x = x / max(x(:));
mask = true(size(x));
mask([1 end],:) = false; mask(:,[1 end]) = false;
im pl 2 2
im(1, mask+x, 'mask + x'), cbar

kappa = 2 - double(x > 0.5);
kappa = mask .* kappa;
im(2, kappa, 'kappa'), cbar

C = Csparse('b2info', kappa, mask, 0);
if im, im subplot 3, spy(C), end

tmp = embed(sum(abs(C))', mask);
im(4, tmp, '|C|''1'), cbar
