 function ob = Gnufft(varargin)
%function ob = Gnufft([mask,] args)
%|
%| Construct Gnufft object that computes nonunform samples of the FT
%| of signals with dimensions [(Nd)] approximately via the NUFFT.
%| For exact (but slow) computation, use Gdsft instead.
%|
%| The arguments (args) are simply a cell array of the all the arguments
%| that will be passed to "nufft_init()" in the appropriate order.
%| See Gnufft_test.m for example usage.
%|
%| Alternatively, one can input a struct of the type returned by nufft_init().
%|
%| Basically, you create a system matrix object by calling:
%|	A = Gnufft( ... )
%| and then you can use it thereafter by typing commands like
%|	y = A * x;
%| which will auto-magically evaluate the DSFT samples.
%| This is useful for iterative image reconstruction in MRI.
%|
%| Optional arguments
%|	mask		logical support array
%|
%| out
%|	ob		fatrix2 (default) or Fatrix if Gnufft('Fatrix', ...)
%|
%| Copyright 2003-6-1, Jeff Fessler, University of Michigan

if nargin == 1 && streq(varargin{1}, 'test')
	run_mfile_local('Gnufft_test0')
	run_mfile_local('Gnufft_test')
return
end
if nargin < 1, help(mfilename), error(mfilename), end

arg.class = 'fatrix2';
arg.mask = [];
if ischar(varargin{1})
	arg.class = varargin{1};
	varargin = varargin(2:end);
end
if islogical(varargin{1})
	arg.mask = varargin{1};
	varargin = varargin(2:end);
end

if length(varargin) ~= 1
	help(mfilename), error 'bad number of arguments'
end

% cell case with nufft_init args
if iscell(varargin{1})
	arg.arg = varargin{1};
	arg.st = nufft_init(arg.arg{:}); % initialize NUFFT structure

% struct case with nufft_init structure given
elseif isstruct(varargin{1})
	arg.st = varargin{1};

else
	error 'bad argument type'
end

switch arg.class
case 'Fatrix'
	if isempty(arg.mask)
		arg.mask = true([arg.st.Nd 1]); % [(Nd)]
	end
%	arg.new_mask = @Gnufft_new_mask;
	arg.dim = [nrow(arg.st.om) sum(arg.mask(:))]; % M x np
	ob = Fatrix(arg.dim, arg, ...
		'forw', @Gnufft_forw, 'back', @Gnufft_back, ...
		'gram', @Gnufft_gram);
case 'fatrix2'
	odim = nrow(arg.st.om); % M
	forw = @(arg, x) nufft(x, arg.st);
	if isempty(arg.mask)
		arg = rmfield(arg, 'mask'); % avoid potential confusion
		back = @(arg, y) nufft_adj(y, arg.st); % simple
	else
		back = @(arg, y) fatrix2_maskit(arg.mask, nufft_adj(y, arg.st));
	end
	ob = fatrix2('idim', arg.st.Nd, 'odim', odim, 'arg', arg, ...
		'does_many', 1, ...
		'gram', @Gnufft_gram, ...
		'forw', forw, 'back', back);
otherwise
	fail('class')
end


%{
% Gnufft_new_mask()
% this is simply too dangerous
function ob = Gnufft_new_mask(ob, mask)
if ~isequal(size(mask), size(ob.arg.mask))
	error 'new mask size does not match old size'
end
ob.arg.mask = mask;
if isa(ob, 'Fatrix')
	ob.arg.dim = [nrow(ob.arg.st.om) sum(ob.arg.mask(:))]; % M x np
	ob.dim = ob.arg.dim;
else % fatrix2
	ob.mask = mask;
	ob.np = sum(mask(:));
end
%}


% Gnufft_forw(): y = A * x
function y = Gnufft_forw(arg, x)

if size(x,1) == arg.dim(2) % [np (L)]
	x = embed(x, arg.mask); % [(Nd) (L)]
end
y = nufft(x, arg.st); % [M (L)]


% Gnufft_back(): x = A' * y
% in
%	y	[M L]
% out
%	x	[np L]
function x = Gnufft_back(arg, y)

x = nufft_adj(y, arg.st); % [(Nd) L]

%Ld = size(y); Ld = Ld(2:end);
Ns = prod(arg.st.Nd);
x = reshapee(x, Ns, []); % [*Nd *L]
x = x(arg.mask,:); % [np *L]
%x = reshape(x, [Ns Ld]); % [np (L)] % not needed if L can be scalar only!
