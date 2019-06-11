 function X = nufft_table_interp(st, Xk, order, flips, om)
%function X = nufft_table_interp(st, Xk, order, flips, [om])
%|
%| table-based nufft
%|
%| in
%|	st	structure	formed by nufft_init (through nufft_init_table)
%|	Xk	[*Kd nc]	over-sampled DFT coefficients
%|	order	0|1		0th or 1st-order interpolation
%|				default is 0 for backward compatability
%|	flips	0|1		sign flips? (for real interpolator, even N)
%|
%| options
%|	om	[M,1]		frequency locations, overriding st.om
%|
%| out
%|	X	[M,nc]		NUFFT values
%|
%| Copyright 2004-3-30, Jeff Fessler and Yingying Zhang, University of Michigan

if nargin < 2, ir_usage, end
if ~isvar('order') || isempty(order)
	order = 0; % default 0th order for backward compatability
end
if ~isvar('flips') || isempty(flips)
	flips = zeros(size(st.Nd)); % default no flip for backward compatability
end

if nargin < 5
	om = st.om;
end

dd = numel(st.Kd);

% t = omega / gamma
tm = zeros(size(om)); % must be double for mex file
for id=1:dd
	gam = 2*pi / st.Kd(id);
	tm(:,id) = om(:,id) / gam;
end

if size(Xk,1) ~= prod(st.Kd), fail 'Xk size problem', end

class_Xk = class(Xk);
% force Xk to be complex, as needed for pointers in the mex files.
if ~isa(Xk, 'double'), Xk = double(Xk); end % double also needed by mex

nc = size(Xk, 2);
%X = zeros(size(om,1), nc);

arg = {int32(st.Jd), int32(st.Ld), double(tm), int32(order), int32(flips)};

Xk = complexify(Xk);

switch dd
case 1
	X = interp1_table_mex(Xk, st.h{1}, arg{:});

case 2
	Xk = reshape(Xk, [st.Kd nc]);
	X = interp2_table_mex(Xk, st.h{1}, st.h{2}, arg{:});

case 3
	Xk = reshape(Xk, [st.Kd nc]);
	X = interp3_table_mex(Xk, st.h{1}, st.h{2}, st.h{3}, arg{:});

case 4
	Xk = reshape(Xk, [st.Kd nc]);
	X = interp4_table_mex(Xk, st.h{1}, st.h{2}, st.h{3}, st.h{4}, arg{:});

otherwise
	fail '> 4d not done'
end

% apply phase shift
if isfield(st, 'phase_shift') && ~isempty(st.phase_shift)
	X = X .* repmat(st.phase_shift, [1 nc]);
end

X = cast(X, class_Xk); % if input is single, then make output match
