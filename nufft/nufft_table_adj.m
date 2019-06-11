 function Xk = nufft_table_adj(st, X, order, flips, om)
%function Xk = nufft_table_adj(st, X, order, flips, om)
%| adjoint of table-based nufft interpolation.
%|
%| in
%|	st		structure from nufft_init
%|	X [M,nc]	DTFT values (usually nc=1)
%|	order	0|1	0th or 1st-order interpolation
%|			default is 0 for backward compatability
%|	flips	0|1	sign flips? (for real table with even N)
%|	om [M,1]	optional (default st.om)
%|
%| out
%|	Xk [*Kd,nc]	DFT coefficients
%|
%| Copyright 2004-3-30, Jeff Fessler and Yingying Zhang, University of Michigan

if nargin < 2, ir_usage, end

if ~isvar('order') || isempty(order)
	order = 0; % default 0th order for backward compatability
end

if ~isvar('flips') || isempty(flips)
	flips = zeros(size(st.Nd)); % default no flips for backward compatability
end

if nargin < 4
	om = st.om;
end

dd = length(st.Kd);

tm = zeros(size(om)); % must be double for mex file!
for id=1:dd
	gam = 2*pi / st.Kd(id);
	tm(:,id) = om(:,id) / gam; % t = omega / gamma
end

if size(X,1) ~= size(om,1)
	fail 'X size problem'
end

nc = size(X,2);

% adjoint of phase shift
class_X = class(X);
if isfield(st, 'phase_shift') && ~isempty(st.phase_shift)
	X = X .* repmat(conj(st.phase_shift), [1 nc]);
end

% convert X to complex double for mex file
if ~isa(X, 'double'), X = double(X); end

arg = {int32(st.Jd), int32(st.Ld), double(tm), int32(st.Kd(1:dd)), ...
	int32(order), int32(flips)};

X = complexify(X);

switch dd
case 1
	Xk = interp1_table_adj_mex(X, st.h{1}, arg{:});

case 2
	Xk = interp2_table_adj_mex(X, st.h{1}, st.h{2}, arg{:});

case 3
	Xk = interp3_table_adj_mex(X, st.h{1}, st.h{2}, st.h{3}, arg{:});

case 4
	Xk = interp4_table_adj_mex(X, st.h{1}, st.h{2}, st.h{3}, st.h{4}, arg{:});

otherwise
	fail '> 4d not done'
end
Xk = cast(Xk, class_X);
