 function [Af, Aa, dif] = test_adjoint(A, varargin)
%function [Af, Aa, dif] = test_adjoint(A, varargin)
%|
%| test that the adjoint of a 'matrix object' A really is its transpose
%|
%| in
%|	A		Fatrix or fatrix2 typically
%|
%| option
%|	'tolre'	dbl	tolerance (relative) for equality (default: 0)
%|	'big'	0|1	for big objects, use random data
%|	'tol'	dbl	tolerance for adjoint for 'big' option (mpd)
%|	'nrep'	int	multiple realizations for 'big' option
%|	'complex' 0|1	if 1, test with complex data (default: 0)
%|	'warn'	0|1	if 1, just print warning. 0 (default), fail if >tol
%|	'chat'	0|1	default 0
%|
%| Copyright 2005-8-2, Jeff Fessler, University of Michigan

if nargin < 1, ir_usage, end

arg.big = 0;
arg.tolre = 0;
arg.tol = 1e-5;
arg.nrep = 1;
arg.chat = 0;
arg.warn = false;
arg.complex = 0;
arg = vararg_pair(arg, varargin);

test_adjoint_big(A, arg.nrep, arg.tol, arg.complex, arg.warn, arg.chat);
if arg.big, return, end % stop here if 'big'

[nd np] = size(A);

try
	% must be very careful here because ob(:,:) calls full(ob) by default
	% which in turns calls fatrix2_subsref_colon() which is "too" clever
	% because it uses the smaller dimension, so transpose has no effect.
	% 2013-07-16 added 'col' option to 'full' to force using the column
	% dimension for the loop (for both ob and ob') thereby truly testing.

	if arg.chat, printm 'A(:,:)', end
	if isa(A, 'Fatrix')
		Af = A(:,:); % not for fatrix2! calls full()
	else
		Af = full(A, 'col'); % trick
	end
	if arg.chat, printm 'done', end

	if arg.chat, printm 'A''(:,:)', end
	Aa = A';
	if isa(Aa, 'Fatrix')
		Aa = Aa(:,:); % not for fatrix2! calls full()
	else
		Aa = full(Aa, 'col'); % trick
	end
	if arg.chat, printm 'done', end

catch
	warn 'A(:,:) or its adjoint failed!?'

	Af = zeros(nd,np);
	Aa = A';
	% forward projection
	for jj=1:np
		x = zeros(np,1);
		x(jj) = 1;
		Af(:,jj) = A * x(:);
	end

	% back projection
	for ii=1:nd
		y = zeros(nd, 1);
		y(ii) = 1;
		Aa(:,ii) = A' * y(:);
	end
end

dif = Af - Aa';
relerr = max(abs(dif(:))) / max(abs(Af(:)));
if relerr > arg.tolre
	printm('adjoint of %s is imperfect (relerr=%g):', inputname(1), relerr)
	printm('adjoint real minmax: %g %g', minmax(real(dif)).')
	printm('adjoint imag minmax: %g %g', minmax(imag(dif)).')
	printm('adjoint real max diff = %g%%', max_percent_diff(real(Af), real(Aa')))
	printm('adjoint imag max diff = %g%%', max_percent_diff(imag(Af), imag(Aa')))
	printm('adjoint class=%s range: %g %g', class(Aa), minmax(Aa(:)).')
else
	tmp = class(A);
	if streq(tmp, 'Fatrix') || streq(tmp, 'fatrix2')
		tmp = sprintf('%s:%s', tmp, A.caller);
	end
	if arg.tolre
		printm('adjoint of %s matches within %g < %g, %s', ...
			inputname(1), relerr, arg.tolre, tmp)
	else
		printm('adjoint of %s appears numerically exact, %s', ...
			inputname(1), tmp)
	end
end


% test_adjoint_big()
% test for big operators using random vectors
function test_adjoint_big(A, nrep, tol, do_complex, do_warn, chat)
rng(0)
for ii=1:nrep
	x = rand(size(A,2),1) - 0.5;
	y = rand(size(A,1),1) - 0.5 ;
	if do_complex
		x = x + 1i * rand(size(A,2),1) - 0.5;
		y = y + 1i * rand(size(A,1),1) - 0.5 ;
	end
	Ax = A * x;
	if ~isreal(Ax) && ~do_complex
		fail 'must test complex systems with ''complex'' option'
	end
%	v1 = y' * (A * x);
	v1 = ir_dot_double(conj(y), A * x);
%	v2 = (x' * (A' * y))';
	v2 = ir_dot_double(conj(x), A' * y)';

	mpd = max_percent_diff(v1, v2);
	if mpd/100 > tol
		pr v1
		pr v2
		pr '[mpd/100 tol]'
		if do_warn
			warn 'adjoint mismatch'
		else
			fail 'adjoint mismatch'
		end

	elseif chat
		pr v1
		pr v2
	end
end
