 function n = nrms(x, xtrue, arg)
%function n = nrms(x, xtrue, arg)
%|
%| normalized rms error
%|
%| in
%|	x	[[Nd] Nrep]	multiple x's for one xtrue
%|	xtrue	[[Nd]]
%|	arg			use '2' to do along 2nd dim for same size
%|
%| out
%|	scalar
%|	vector if x and xtrue have same size and arg = 2

if nargin == 1 && streq(x, 'test'), nrms_test, return, end

if nargin < 2, ir_usage, end

if nargin == 3 && arg == 2
	if any(size(x) ~= size(xtrue))
		error 'need same size'
	end
	n = sqrt(mean(abs(xtrue - x).^2));
return
end

xtrue = double(xtrue(:));
np = numel(xtrue);
nrep = numel(x) / np;
x = reshape(x, np, nrep);
x = double(x);
n = abs(x - repmat(xtrue, [1 nrep]));

xnorm = norm(xtrue);
enorm = sqrt(sum(n.^2)');
if xnorm
	n = enorm / xnorm;
elseif ~enorm
	n = 0; % 0/0 = 0
else
	n = inf;
end

doprint = nargout < 1;
base = '';
if doprint
	fprintf('%snrms(%s,%s) =', base, inputname(1), inputname(2))
	for ii=1:length(n)
		fprintf(' %g%%', n(ii) * 100)
	end
	fprintf('\n')
	if ~nargout
		clear n
	end
end

function nrms_test
x = ones(5,1);
y = 1.1*ones(5,1);
z = 1.2*ones(5,1);
nrms(y,x)
nrms([y z], x)
