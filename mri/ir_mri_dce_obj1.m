 function [dyn_obj dce] = ir_mri_dce_obj1(varargin)
%function [dyn_obj dce] = ir_mri_dce_obj1(varargin)
%|
%| Generate dynamic object for DCE MRI simulations of GRE sequence.
%|
%| option
%|	'TR'			[s] repetition time (default: 5e-3)
%|	'TE'			[s] echo time (default: 0)
%|	'flip'			[rad] flip angle (default: pi/6)
%|	'duration_s'		[s] (default: 4*60=240)
%|	'n_tr_round'		make total #TR a multiple of this (default: 24)
%|				for ease dividing up scan into various "frames"
%|	'n_tr_merge'		# of TR to "merge" into one "time point"
%|				to reduce memory (default: 100)
%|	'labels' [(N)]		brainweb labels image
%|				(default: brainweb slice 93)
%|
%| out
%|	dyn_obj	[(N) Nt]	[au] dynamic object magnetization
%|				Nt = (duration_s / TR) / n_tr_merge
%|	dce	struct		.TR .TE etc.
%|
%| 2014-08-21 Jeff Fessler and Mai Le, University of Michigan

if nargin == 1 && streq(varargin{1}, 'test'), ir_mri_dce_obj1_test, return, end
if ~nargout && nargin < 1, ir_usage, end

arg.TR = 5e-3; % 5ms TR
arg.TE = 0;
arg.flip = pi/6; % 30 degree flip
arg.duration_s = 4*60; % [s]
arg.n_tr_round = 24; % round duration down to be a multiple of this # of TRs
arg.n_tr_merge = 100; % # of TR to "merge" into one "time point"
arg.labels = []; % [(N)] brainweb values are 0 to 10
arg.Ktrans = [0.2 0.6 3.0];
arg.kep = [1.3 2.0 6.0];
arg.lesion_legend = {'slow', 'moderate', 'rapid'};
arg.phase = 'bilin'; % spatially smooth image phase to stress test
arg.brainweb_dir = '';
arg.brainweb_slice = 93;
arg.chat = 0;
arg = vararg_pair(arg, varargin);

dce = arg;

n_tr_init = arg.duration_s / arg.TR; % initial # of TR
n_tr_mult = arg.n_tr_round * arg.n_tr_merge;
n_tr = n_tr_mult * floor(n_tr_init / n_tr_mult); % actual # of TR
if arg.chat && n_tr ~= n_tr_init
	printm('#TR rounded down from %d to %d to be multiple of %d', ...
		n_tr_init, n_tr, n_tr_mult);
end

if isempty(arg.labels)
	if isempty(arg.brainweb_dir)
		arg.brainweb_dir = '~/l/dat/phantom,digital/brainweb/ms/';
	end
	img_path = [arg.brainweb_dir 'phantom_1.0mm_msles2_crisp.fld'];
	if exist(img_path, 'file')
		labels = fld_read(img_path, 'slice', arg.brainweb_slice);
		labels = flipdim(labels,2); % patient facing up
	else
		labels = ellipse_im([181 217], 'shepplogan-brainweb');
	%	fail('cannot find brainweb phantom "%s"', img_path)
	end
%	unique(labels(:))
	labels = labels(2:end,2:end); % even size
	labels = padn(labels, [2*2*3*16 216]); % more fft friendly
	[nx, ny] = size(labels);
%	pr factor(nx)
%	pr factor(ny)
	im clf, im(labels), cbar

	if 1 % add lesions as additional labels
		ell = [-32 30 9 9 0 1;
			32 15 7 7 0 1;
			-27 -50 4 4 0 1 ];
		nlesion = nrow(ell);
		assert(nlesion == numel(arg.Ktrans), 'nlesion')
		jf_equal(size(arg.kep), size(arg.Ktrans))
		nbw = max(labels(:)); % # of brainweb tissue types
		% recall brainweb starts at 0 for background
		if arg.chat
			printm('full frame sampling time %g [s]', nx*ny*arg.TR)
		end
		ig = image_geom('nx', nx, 'ny', ny, 'dx', 1);
		for it = 1:nlesion
			tmp = ellipse_im(ig, ell(it,:));
		%	unique(tmp)
			labels(tmp ~= 0) = nbw + it;
		end
		im(ig.x, ig.y, labels)
	else
		nlesion = 0;
	end

	% tissue class mr properties
	nclass = nbw + nlesion;
	tcmr = struct('Ktrans',{}, 'kep',{}, ...
		'name',{}, 'pd',{}, 'r1',{}, 'r2',{}, 'r2s',{});
	for ic = 1:nclass
		if ic <= nbw
			tmp = mri_brainweb_params(ic);
			tcmr(ic).Ktrans = 0;
			tcmr(ic).kep = 0;
		else
			tcmr(ic).Ktrans = arg.Ktrans(ic-nbw);
			tcmr(ic).kep = arg.kep(ic-nbw);
			tmp = mri_brainweb_params('white-matter');
		end
		tcmr(ic).name = tmp.name;
		tcmr(ic).pd = tmp.pd;
		tcmr(ic).r1 = div0(1000, tmp.t1);
		tcmr(ic).r2 = div0(1000, tmp.t2);
		tcmr(ic).r2s = div0(1000, tmp.t2s);
	end

	[dyn_obj, dce] = ir_mri_dce_obj1_do(labels, tcmr, dce, ...
		arg.TR, arg.TE, arg.flip, n_tr);
	dce.labels = labels; % save for later ROI use

else
	fail 'only default brainweb done'
end

if nx*ny > n_tr
	warn('nx*ny=%d > n_tr=%d', nx*ny, n_tr)
end

% add (optional) smooth phase to make simulation more realistic
switch arg.phase
case {'', 'none'}
case 'bilin'
	x = [1:nx]/nx;
	y = [1:ny]/ny;
	[xx yy] = ndgrid(x, y);
	phase =	exp((2i * pi) * (xx.*yy));
	dyn_obj = dyn_obj .* repmat(single(phase), [1 1 size(dyn_obj,3)]);

otherwise
	fail('unknown phase pattern: %s', arg.phase)
end



% ir_mri_dce_obj1_do()
function [dyn_obj, dce] = ir_mri_dce_obj1_do(labels, tcmr, dce, ...
		TR, TE, flip, n_tr)

dce.Ktrans = [tcmr(:).Ktrans]';
dce.kep = [tcmr(:).kep]';
% pr n_tr
dce.ti = (([0:((n_tr/dce.n_tr_merge)-1)] + 0.5) * dce.n_tr_merge * TR) / 60; % [min]
[dce.Ct dce.Cp] = ir_mri_dce_aif1(dce.ti, dce.Ktrans, dce.kep);
% plot(dce.ti, dce.Ct)

dce.r10 = [tcmr(:).r1]';
dce.r1d = ir_mri_dce_r1d(dce.r10, dce.Ct);
% plot(dce.ti, dce.r1d)

dce.r2s = [tcmr(:).r2s]'; % R2* in [1/s]
dce.sig = ir_mri_dce_gre(dce.r1d, 'TR', TR, 'TE', TE, 'r2s', dce.r2s);
% plot(dce.ti, dce.sig)

Nt = numel(dce.ti); % # of merged TR ("fine" time sampling)
dyn_obj = zeros([prod(size(labels)) Nt], 'single'); % [(*N) Nt]

% time curve for each tissue class
labs = unique(labels(:));
labs(labs == 0) = []; % ignore background
for il=1:numel(labs)
	ic = labs(il);
	tmp = labels == ic;
%	im(tmp), drawnow
%	dyn_obj = dyn_obj + tmp(:) * dce.sig(ic,:);
	dyn_obj(tmp,:) = repmat(dce.sig(ic,:), [sum(tmp(:)) 1]);
end

dyn_obj = reshape(dyn_obj, [size(labels) Nt]); % [(N) Nt]


% ir_mri_dce_obj1_test()
function ir_mri_dce_obj1_test
[dyn_obj, dce] = ir_mri_dce_obj1('chat', 1);
if im, pr dce, end
im(abs(dyn_obj(:,:,1:25:end)))
