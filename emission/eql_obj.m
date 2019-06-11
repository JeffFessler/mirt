  function obj = eql_obj(xs, G, yi, ci, ri, R, mask)
%|function obj = eql_obj(xs, G, yi, ci, ri, R, mask)
%| compute quadratically penalized emission Poisson likelihood
%| for each "column" of x (but x can be multi-dim if mask provided)
%|
%| Copyright May 1999, Jeff Fessler

if nargin < 3, help(mfilename), error(mfilename), end
if ~isvar('ci') || isempty(ci)
	ci = ones(size(yi(:)));
end
if ~isvar('ri') || isempty(ri)
	ri = zeros(size(yi(:)));
end
if ~isvar('R'), R = []; end
if nargin == 7
	dims = size(xs);
	xs = reshape(xs, prod(dims(1:(end-1))), dims(end));
	xs = xs(mask(:), :);
end

like = zeros(ncol(xs),1);
penal = zeros(ncol(xs),1);
for kk=1:ncol(xs)
	ticker(mfilename, kk, ncol(xs))
	x = xs(:,kk);
	yp = ci .* (G * x) + ri; % predicted measurements
	if any(yp(:) == 0 & yi(:) ~= 0)
		error 'bug'
	end
	good = yi > 0;
%	like(kk) = sum(yi .* log(yp) - yp); % old basic way
	like(kk) = sum(yi(good) .* log(yp(good))) - sum(yp); % new better way
	if ~isempty(R)
		if isstruct(R)
			penal(kk) = R.penal(R, x);
		else
			penal(kk) = norm(col(R*x))^2/2; % norm(C * x).^2 / 2;
		end
	end
end
obj = like - penal;
