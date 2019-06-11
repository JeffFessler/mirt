 function yi = interp1_jump(xj, yj, xi, varargin)
%function yi = interp1_jump(xj, yj, xi, {arguments for interp1})
%|
%| Generalization of matlab's "interp1" to allow xj with repeated values,
%| for interpolation of a function that has "jumps" (discontinuities),
%| such as is caused by k-edges for mass attenuation coefficients.
%| If the first option is 'monodown' then the function is expected to
%| be monotone decreasing except for certain jumps that may not correspond
%| to xj's with equal values.
%| The (remaining) options are passed to 'interp1'.
%|
%| Copyright 2004-5-2, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(xj, 'test'), interp1_jump_test, return, end
if nargin < 3, ir_usage, end

if length(xj) ~= length(yj), error 'xj and yj have different lengths', end

jjump = find(diff(xj) == 0);
if length(varargin) && streq(varargin{1}, 'monodown')
	varargin = {varargin{2:end}};
	yjump = find(diff(yj) > 0);
	jjump = unique([jjump; yjump]);
end
npiece = 1 + length(jjump);
if npiece == 1
	yi = interp1(xj, yj, xi, varargin{:});
	return
end

yi = zeros(size(xi));
done = zeros(size(xi));
for ip=1:npiece
	if ip == 1
		jlist = [1:jjump(1)];
	elseif ip == npiece
		jlist = [(1+jjump(npiece-1)):length(xj)];
	else
		jlist = [(1+jjump(ip-1)):jjump(ip)];
	end
	x = xj(jlist);
	y = yj(jlist);
	if ip == 1
		ilist = find(xi <= max(x));
	elseif ip == npiece
		ilist = find(min(x) <= xi);
	else
		ilist = find(min(x) <= xi & xi <= max(x));
	end
	if isempty(ilist), continue, end

	if length(x) > 1
		yi(ilist) = interp1(x, y, xi(ilist), varargin{:});
		done(ilist) = 1;
	elseif length(x) == 1
		if any(x == xi(ilist))
			yi(ilist) = y(x == xi(ilist));
		else
			warning 'bug?'
			keyboard
		end
	end
end

% for anything left over, use linear interpolation
if any(~done)
	[xj jj] = unique(xj);
	yi(~done) = interp1(xj, yj(jj), xi(~done), 'linear');
end


function interp1_jump_test
x = [0 0.5 1 1 2 3 3 4 5];
%x = [0 0.5 1 1.1 2 3 3.1 4 5];
y = [1 0 0 1 1 2 1 2 2];
t = linspace(-0.5,0.5+max(x),1001);
f = interp1_jump(x, y, t, 'pchip', 'extrap');
if im
	clf, subplot(211)
	plot(x, y, 'o', t, f, '-')
end

x = [0:5];
y = [4 2 1 2 1 0];
f = interp1_jump(x, y, t, 'monodown', 'pchip', 'extrap');
if im
	subplot(212)
	plot(x, y, 'o', t, f, '-')
end
