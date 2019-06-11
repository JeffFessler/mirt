 function ob = fatrix2_mtimes2(ob1, ob2, varargin)
%function ob = fatrix2_mtimes2(ob1, ob2, options)
%|
%| Construct fatrix2 object that is the product of two objects:
%|	ob = ob1 * ob2.
%| Requires size(ob1,2) == size(ob2,1) as in matrix multiplication.
%|
%| in
%|	ob1	*atrix		any object that can do "mtimes" and "size"
%|	ob2	*atrix		''
%|
%| out
%|	ob	fatrix2		ob1 * ob2
%|
%| caution: this will have limited (if any) "block" capabilities
%| that are inherited from ob1 (if it has any):
%|	ob{iblock} = ob1{iblock} * ob2		todo: under construction
%|
%| Copyright 2010-12-04, Jeff Fessler, University of Michigan

if size(ob1, 2) ~= size(ob2, 1)
	pr 'size(ob1,2)'
	pr 'size(ob2,1)'
	fail('size mismatch')
end

if isfield(ob1, 'idim') && isfield(ob2, 'odim') ...
	&& ~isequal(ob1.idim, ob2.odim)
	if numel(ob1.idim) > 1 && ob1.idim(end) == 1
		ob1.idim = ob1.idim(1:end-1); % remove trailing "1"
	end
	if numel(ob2.odim) > 1 && ob2.odim(end) == 1
		ob2.odim = ob2.odim(1:end-1); % remove trailing "1"
	end
	if ~isequal(ob1.idim, ob2.odim)
		fail('idim=[%s] odim=[%s] mismatch', ...
			num2str(ob1.idim), num2str(ob2.odim))
	end
end

if isa(ob2, 'fatrix2')
	idim = ob2.idim;
	imask = ob2.imask;
else
	idim = size(ob2,2);
	imask = [];
end

if isa(ob1, 'fatrix2') % && (numel(idim) > 1 || ~isempty(imask)) % caution
	odim = ob1.odim;
	omask = ob1.omask;
else
	odim = size(ob1,1);
	omask = [];
end

%{

2016-10-24 commented ok; out to have block objects but the cascade will not,
at least for now

if isfield(ob2, 'nblock') && ~isempty(ob2.nblock)
	fail 'block object on rhs unsupported'
end

if isfield(ob1, 'nblock') && ~isempty(ob1.nblock)
	fail 'block object on lhs could be added - ask jf'
%{
	forw_block = @(arg, x, iblock, nblock)
		arg.ob1{iblock} * (arg.ob2 * x)
	arg.ob1.handle_forw_block(arg.ob1.arg, arg.ob2 * x, iblock, nblock);
	back_block = @(arg, y, iblock, nblock)
	arg.ob2' * arg.ob1.handle_back_block(arg.ob1.arg, arg.ob2 * x, iblock, nblock);
	argblock = {'blockify_data', ob1.blockify_data, ...
		'block_setup', @fatrix2_block_setup, ...
		'mtimes_block', @fatrix2_mtimes_block};
%}
end

%block_args = {};

%}

if isa(ob1, 'fatrix2') && isa(ob2, 'fatrix2') % fatrix * fatrix
%{
	forw = @(arg, x) (arg.ob1.scale * arg.ob2.scale) * ...
		arg.ob1.handle_forw(arg.ob1.arg, ...
		arg.ob2.handle_forw(arg.ob2.arg, x));
	back = @(arg, y) conj(arg.ob1.scale * arg.ob2.scale) * ...
		arg.ob2.handle_back(arg.ob2.arg, ...
		arg.ob1.handle_back(arg.ob1.arg, y));
%}
	forw = @(arg, x) fatrix2_do_forw(arg.ob1, fatrix2_do_forw(arg.ob2, x));
	back = @(arg, y) fatrix2_do_back(arg.ob2, fatrix2_do_back(arg.ob1, y));

%	block_args = fatrix2_mtimes2_block_args(ob1, ob2);

else % product of something else with fatrix
	forw = @(arg, x) arg.ob1 * (arg.ob2 * x);
	back = @(arg, y) arg.ob2' * (arg.ob1' * y);
	if isa(ob1, 'fatrix2') || isa(ob2, 'fatrix2')
		warn('multiplying %s * %s may fail', class(ob1), class(ob2))
	end
end

arg.ob1 = ob1;
arg.ob2 = ob2;

ob = fatrix2('arg', arg, ...
	'idim', idim, 'imask', imask, ...
	'odim', odim, 'omask', omask, ...
	'forw', forw, 'back', back, ...
	'power', @fatrix2_mtimes2_power, ...
	'sparse', @fatrix2_mtimes2_sparse);
%	block_args{:}, ...
%	arg.block{:},


%{

% fatrix2_mtimes2_block_args()
% kludge to handle Gdiag * fatrix
function block_args = fatrix2_mtimes2_block_args(ob1, ob2)
block_args = {};
if ~streq(ob1.caller, 'Gdiag') || ~isequal(ob1.idim, ob2.odim) ...
	return
end
if ~isempty(ob2.handle_forw_block)
	block_args = {'forw', @fatrix2_mtimes2_block_Gdiag_forw};
end
if ~isempty(ob2.handle_back_block)
	block_args = {block_args{:}, 'back', @fatrix2_mtimes2_block_Gdiag_back};
end


function y = fatrix2_mtimes2_block_Gdiag_forw(arg, x, iblock, nblock)
y = arg.ob1.scale * arg.ob2.scale * ...
	arg.ob2.handle_forw_block(arg.ob2.arg, x, iblock, nblock);
ia = istart:nblock:arg.ob2.odim(end);

	back = @(arg, y) conj(arg.ob1.scale * arg.ob2.scale) * ...
		arg.ob2.handle_back(arg.ob2.arg, ...
		arg.ob1.handle_back(arg.ob1.arg, y));
todo

function x = fatrix2_mtimes2_block_Gdiag_back(arg, x, iblock, nblock)
todo

%}



% fatrix2_mtimes2_power(): A.^p
function ob = fatrix2_mtimes2_power(ob, pow)
ob1 = ob.arg.ob1;
if isfield(ob1, 'caller') && streq(ob1.caller, 'diag', 4)
	ob.arg.ob1 = ob.arg.ob1 .^ pow;
	ob.arg.ob2 = ob.arg.ob2 .^ pow;
else
	error 'power defined only for diag * object'
end


% fatrix2_mtimes2_sparse(): sparse(A)
function ob = fatrix2_mtimes2_sparse(ob)
ob = sparse(ob.arg.ob1) * sparse(ob.arg.ob2);
