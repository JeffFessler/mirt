 function kspace = ir_mri_dce_kspace1(dyn_obj, samp, varargin)
%function kspace = ir_mri_dce_kspace1(dyn_obj, samp, varargin)
%|
%| Generate dynamic kspace data for DCE MRI simulations.
%|
%| in
%|	dyn_obj	[(N) Nt]	complex dynamic object
%|	samp	[(N) Nt]	sampling pattern
%|				(reshaped from [(N) Nt/Nframe Nframe])
%|
%| option
%|	Nframe			divisor of Nt for "# of frames" (default: Nt)
%|	'smap'	[(N) Ncoil]	coil sensitivity maps (default: [] single coil)
%|
%| out
%|	kspace	[(N) Nframe Ncoil]	k-space data, with DC in "middle"
%|
%| 2014-08-25 Jeff Fessler and Mai Le, University of Michigan

if nargin == 1 && streq(dyn_obj, 'test'), ir_mri_dce_kspace1_test, return, end
if nargin < 2, ir_usage,, end

arg.Nframe = []; % default to Nt below
arg.smap = []; % default to single coil
arg.chat = 0;
arg = vararg_pair(arg, varargin);

if isempty(arg.Nframe)
	arg.Nframe = size(dyn_obj, ndims(dyn_obj));
end

kspace = ir_mri_dce_kspace1_do(dyn_obj, samp, arg.Nframe, arg.smap, arg.chat);


function kspace = ir_mri_dce_kspace1_do(dyn_obj, samp, Nframe, smap, chat)

dims = size(dyn_obj);
Np = dims(1:(end-1));
Nt = dims(end);

jf_equal(size(samp), size(dyn_obj))
if ~islogical(samp), fail 'samp must be logical', end

if ~isempty(smap)
	tmp = size(smap);
	Ncoil = tmp(end);
	jf_equal(tmp(1:end-1), Np)
else
	Ncoil = 1;
end

% check that Nframe is divisor of Nt and unique sampling within rep
if Nframe ~= Nt
	nsplit = Nt / Nframe;
	if nsplit ~= round(nsplit)
		fail('Nframe=%d not divisor of Nt=%d', Nframe, Nt)
	end
	samp1 = reshape(samp, prod(Np), [], Nframe);
	samp1 = squeeze(sum(samp1, 2));
	samp1 = reshape(samp1, [Np Nframe]);
%	im('row', 6, tmp)
	if any(samp1(:) > 1), fail 'sampling not unique per rep', end
	samp1 = samp1 == 1;
end

kspace = complex(zeros(prod(Np), Nframe, Ncoil, 'single'));

ticker reset
for ic=1:Ncoil
	ticker(mfilename, ic, Ncoil)
	if isempty(smap)
		x = dyn_obj;
	else
		tmp = stackpick(smap, ic);
		x = repmat(tmp, [1 1 Nt]) .* dyn_obj;
	end

	switch ndims(dyn_obj)
	case 3 % 2D+T
		kspace_x = ir_fftshift2(fft2(ir_fftshift2(x))); % [nx ny Nt]
	otherwise
		fail 'only 2d done'
	end

	kspace_x = kspace_x .* samp; % apply sampling pattern

	kspace_x = reshape(kspace_x, prod(Np), [], Nframe); % [*N (Nt/Nframe) Nframe]
	for ir=1:Nframe
		% here is where realistic "model mismatch" occurs
		% trick: this "sum" needs unique samples within each rep
		kspace(:,ir,ic) = sum(kspace_x(:,:,ir), 2);
	end
end

kspace = reshape(kspace, [Np Nframe Ncoil]);


% ir_mri_dce_kspace1_test()
function ir_mri_dce_kspace1_test

[dyn_obj, dce] = ir_mri_dce_obj1('chat', 0);
Nframe = 12;
[samp1 samp2] = ir_mri_dce_samp1(size(dyn_obj), ...
	'n_tr_merge', dce.n_tr_merge, 'Nframe', Nframe, 'chat', 1);
nx = size(dyn_obj,1);
ny = size(dyn_obj,2);
smap = ir_mri_sensemap_sim('nx', nx, 'ny', ny, 'ncoil', 2, 'coil_distance', 2);
% im(smap)
kspace = ir_mri_dce_kspace1(dyn_obj, samp2, 'Nframe', Nframe, 'smap', smap, 'chat', 0);
% pr size(kspace)

ir=1:4:Nframe;
tmp = log(abs(kspace(:,:,ir,2)));
tmp(~samp1(:,:,ir)) = -8; % avoid warning about infty
im('row', 1, tmp, [-8 0] + max(tmp(:)), 'log(|kspace|)')
