function [xs, precon] = eml_osem(x, Gb, yi, ci, ri, varargin)
%function xs = eml_osem(x, Gb, yi, ci, ri, [options])
% E-ML-OSEM algorithm for image reconstruction from Poisson emission data
% (ordered subsets expectation maximization)
% model: Y_i ~ Poisson(c_i [G x]_i + r_i)
% in
%	x	[np,1]		initial estimate
%	Gb	[nd,np]		Gblock object (see eml_osem_test.m)
%	yi,ci,ri [nb,na]	see em_fbp.m (for model too)
% option
%	'niter			# iterations (default: 1+1)
%	'isave			which iterations to archive (default: all)
%	'pixmax			upper constraint for pixel values
%	'precon [np,{1|nblock}]	preconditioners.  recommend: 'classic'
%	'relax0	[1] or [2]	relax0 or (relax0, relax_rate)
% out
%	xs [np,length(isave)]	updated image vectors each iteration
%
% Modification of the file eml_osem_sse.m of Jeff Fessler

if nargin < 3, help(mfilename), error(mfilename), end

% backward compatability with old usage
% eml_osem(x, Gb, yi, ci, ri, niter, pixmax, precon, relax0)
if length(varargin) && isnumeric(varargin{1})
	[xs, precon] = eml_osem_orig(x, Gb, yi, ci, ri, varargin{:});
return
end

% defaults
arg.niter = 1;
arg.isave = [];
arg.pixmax = inf;
arg.chat = false;
arg.relax0 = 1;
arg.relax_rate = 0;
arg.precon = 'classic';
arg = vararg_pair(arg, varargin);

if isempty(arg.isave), arg.isave = 0:arg.niter; end

Gb = block_ob(Gb, 'ensure'); % make it a block object (if not already)
nblock = block_ob(Gb, 'n');
starts = subset_start(nblock);

if ~isvar('ci') | isempty(ci)
	ci = ones(size(yi));
end
if ~isvar('ri') | isempty(ri)
	ri = zeros(size(yi));
end

if length(arg.relax0) == 1
	arg.relax_rate = 0;
elseif length(arg.relax0) == 2
	arg.relax_rate = arg.relax0(2);
	arg.relax0 = arg.relax0(1);
else
	error relax
end

eml_check(yi, ci, ri, 'os', nblock);
[nb na] = size(yi);

%
% Test of SSE projectors
%
% Asume the following - (information to be contained in Gb)
%
nthread = 1;
offset_r = 0.0;
fov = 500;
nx = Gb.arg.idim(1);
ny = Gb.arg.idim(2);
dx = fov/nx;
dy = -dx;
offset_x = 0.0;
offset_y = 0.0;
orbit_low = 0.0;
orbit_high = pi;
dr = 528 / (nx+2);
strip_width = dr;

% new initial estimate
x = ones(nx,ny);
%
% precompute the preconditioners for classic OSEM
%
if streq(arg.precon, 'classic')
 %  precon = zeros(size(Gb,2), nblock);
  precon = zeros(nx,ny, nblock);
  for iset=1:nblock
    ticker([mfilename ' : precon'], iset, nblock)
    istart = starts(iset);
    ia = istart:nblock:na;
    
    Asum = jmh_sse_mex('jmh,sse,back',int32(1), int32(nb), single(offset_r), ...
                        single(dr),single(dx), single(dy), single(offset_x), ...
                        single(offset_y), uint8(Gb.arg.mask), int32(na), ...
                        single(orbit_low), single(orbit_high), int32(nblock), ...
                        int32(iset-1), reshape(col(ci(:,ia)),[nb,na/nblock]));
    
 %    Asum = Gb{istart}' * col(ci(:,ia));
    Asum(Asum == 0) = Inf; % avoid divide by 0
    precon(:,:, iset) = 1 ./ Asum;
  end, clear Asum
  
%
% This 'fast' approach uses the same preconditioner for each subset
% which helps convergence (if diminishing relaxation is used).
% Of course, convergence is not really essential in the unregularized case,
% but this saves a bit of memory by using just one preconditioner, which
% ought to work when the subsets are reasonably balanced.
% But, if diminishing relaxation is not used, then this can cause problems
% because pixels can get stuck at zero.  So I no longer recommend this.
%
elseif streq(arg.precon, 'fast')
	warning 'Using fast preconditioner rather than classic OSEM'
	precon = Gb' * ci(:); % complete backprojection
	precon(precon == 0) = Inf; % avoid divide by 0
	precon = nblock ./ precon;

else % user-supplied precon
	precon = arg.precon;
	if ncol(precon) ~= nblock & ncol(precon) ~= 1
		error 'precon columns must be 1 or nblock'
	end
end


%
% loop over iterations
%
xs = zeros(nx,ny, length(arg.isave));
if any(x <= 0), error 'need x > 0', end
x = min(x, arg.pixmax);
if any(arg.isave == 0)
	xs(:,:,find(arg.isave == 0)) = x;
end

for iter = 1:arg.niter
	relax = arg.relax0 / (1 + arg.relax_rate * (iter-2));

	%
	% loop over subsets
	%
	for iset=1:nblock
		ticker(mfilename, [iter iset], [arg.niter nblock])

		istart = starts(iset);
		ia = istart:nblock:na;

		% predicted measurements
                yp = jmh_sse_mex('jmh,sse,proj',int32(1), int32(nb), single(offset_r), ...
                        single(dr),single(dx), single(dy), single(offset_x), ...
                        single(offset_y), uint8(Gb.arg.mask), int32(na), ...
                        single(orbit_low), single(orbit_high), int32(nblock), ...
                        int32(iset-1), single(reshape(x,[nx,ny])));
		
                %yp = reshape(Gb{istart} * x, nb, length(ia));
		
                yp = ci(:,ia) .* yp + ri(:, ia);
		yp(yp == 0) = inf;	% avoids /0 error

		pre = precon(:,:,min(iset,ncol(precon)));	% 1 or iset

		dhi = ci(:,ia) .* (yi(:,ia) ./ yp - 1);
                
                grad = jmh_sse_mex('jmh,sse,back',int32(1), int32(nb), single(offset_r), ...
                                   single(dr),single(dx), single(dy), single(offset_x), ...
                                   single(offset_y), uint8(Gb.arg.mask), int32(na), ...
                                   single(orbit_low), single(orbit_high), int32(nblock), ...
                                   int32(iset-1), reshape(dhi,[nb,na/nblock]));

                %grad = Gb{istart}' * dhi(:);

                % Filter outside mask
                x = x.*Gb.arg.mask;
		x = x + relax * (x .* pre) .* grad;

		% fix: implement ahn's shift trick?  no, just use classic.
		nneg = sum(x < 0);
		if nneg, printm('oh no! %d negatives', nneg), end

		x = max(x,0);	% caution: if x becomes zero, it is always zero!
		x = min(x,arg.pixmax);
	end

	if arg.chat, printm('Range %g %g', min(x), max(x)), end
	if any(arg.isave == iter)
		xs(:,:,find(arg.isave == iter)) = x;
	end
end


%
% eml_osem_orig()
% arguments: pixmax, precon, relax0
%
function [xs, precon] = eml_osem_orig(x, Gb, yi, ci, ri, niter, varargin)
arg = {'niter', niter};
opt = {'pixmax', 'precon', 'relax0'};
for ii=1:length(varargin)
	arg = {arg{:}, opt{ii}, varargin{ii}};
end
[xs, precon] = eml_osem(x, Gb, yi, ci, ri, arg{:});
