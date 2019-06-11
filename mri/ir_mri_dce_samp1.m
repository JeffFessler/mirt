 function [samp1 samp2] = ir_mri_dce_samp1(dims, varargin)
%function [samp1 samp2] = ir_mri_dce_samp1(dims, varargin)
%|
%| Generate dynamic kspace data for DCE MRI simulations.
%|
%| in
%|	dims			dynamic object dims, e.g. [nx ny Nt]
%|				Nt is # of "phase encode groups"
%|
%| option
%|	'pattern'	(char)	(default: 'v01')
%|	'n_tr_merge'	(int)	see ir_mri_dce_obj1 (default: 1)
%|	'n_low'		[2]	size of fully sampled k-space center
%|				(default: [16 16])
%|	'Nframe'	[1]	# of times to fully sample k-space center,
%|				i.e., # of "frames" (default: Nt)
%|				must be a divisor of Nt; affects trade-off
%|				between compute time and samples per frame, etc.
%|
%| out
%|	samp1	[(N) Nframe]	k-space sampling pattern (binary)
%|	samp2	[(N) Nt]	k-space sampling pattern (binary)
%|				reshaped version of [(N) Nt/Nframe Nframe]
%|
%| 2014-08-25 Jeff Fessler and Mai Le, University of Michigan

if nargin == 1 && streq(dims, 'test'), ir_mri_dce_samp1_test, return, end
if nargin < 1, ir_usage, end

arg.pattern = '';
arg.n_low = [16 16];
arg.Nframe = dims(end); % Nt
arg.n_tr_merge = 1; % silly default
arg.seed = 0;
arg.chat = 0;
arg = vararg_pair(arg, varargin);

rng(arg.seed)

Nt = dims(end);
Nd = Nt * arg.n_tr_merge;
Nframe = arg.Nframe;
nsplit = Nt / Nframe;
if nsplit ~= round(nsplit)
	fail('Nframe=%d not divisor of Nt=%d', arg.Nframe, Nt)
end
Np = dims(1:(end-1)); % spatial dimensions (pixels)

switch arg.pattern
case {'', 'v01'}
	samp1 = ir_mri_dce_samp1_v01(Np, ...
		arg.n_low, Nframe, Nd, arg.chat);
otherwise
	fail('only sampling pattern "%s"', arg.pattern)
end

samp1_col = reshape(samp1, [prod(Np) Nframe]); % [*N Nframe]
samp2 = false(prod(Np), Nt/Nframe, Nframe); % [*N Nt/Nframe Nframe]
for ir=1:Nframe
	tmp = find(samp1_col(:,ir)); % length should be multiple of nsplit
	index = randperm(numel(tmp));
	index = reshape(index, [], nsplit);
	for is=1:nsplit
%		samp2(tmp(index(:,is)), (ir-1)*nsplit + is) = true;
		samp2(tmp(index(:,is)), is, ir) = true;
	end
end
samp2 = reshape(samp2, [Np Nt]); % [(N) Nt]


% ir_mri_dce_samp1_v01()
function samp = ir_mri_dce_samp1_v01(dims, n_low, Nframe, Nd, chat)

if numel(dims) ~= numel(n_low)
	fail 'n_low'
end
n_pe_per_frame = Nd / Nframe;
if n_pe_per_frame ~= round(n_pe_per_frame)
	fail('Nframe=%d not divisor of Nd=%d', Nframe, Nd)
end
if chat
	printm('under-sampling = np/n_pe_per_frame = %d/%d = %.1f', ...
		prod(dims), n_pe_per_frame, prod(dims) / n_pe_per_frame)
end
n_high = n_pe_per_frame - prod(n_low);
if n_high < 0
	fail('n_pe_per_frame=%d < nlow=%d', n_pe_per_frame, prod(n_low))
end
if chat
	printm('n_low=%d n_high=%d np=%d', prod(n_low), n_high, prod(dims))
end

Np = prod(dims);
samp = false(Np, Nframe);
low = true(n_low);
samp_lo = col(padn(low, dims));
tmp = [1:prod(dims)]';
index_hi = tmp(~samp_lo);
for ir = 1:Nframe
	tmp = randperm(numel(index_hi), n_high); % todo: poisson disk preferable
	keep_hi = index_hi(tmp);
	samp_hi = false(prod(dims), 1);
	samp_hi(keep_hi) = true;
	if any(samp_lo .* samp_hi), fail 'bug', end
	samp(:,ir) = samp_lo | samp_hi;
end
samp = reshape(samp, [dims Nframe]); % [(N) Nframe]


% ir_mri_dce_samp1_test()
function ir_mri_dce_samp1_test

[dyn_obj dce] = ir_mri_dce_obj1('chat', 0);
Nframe = 48;
[samp1 samp2] = ir_mri_dce_samp1(size(dyn_obj), ...
	'n_tr_merge', dce.n_tr_merge, 'Nframe', Nframe, 'chat', 1);
pr size(samp1)
pr size(samp2)
if 1 % check that samp2 collapse to samp1
	tmp = reshape(samp2, size(dyn_obj,1), size(dyn_obj,2), [], Nframe);
	tmp = squeeze(sum(tmp, 3)) > 0;
	jf_equal(tmp, samp1)
end
im(samp1(:,:,1:8:end))
