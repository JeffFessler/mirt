function y = jf_prctile(pn, x, p, dim)
% mimic matlab's prctile which is in the "stats" toolbox
if streq(x, 'test'), jf_prctile_test, y = []; return, end

if nargin < 3, ir_usage, end

% work with permuted dimensions if needed
if nargin > 3 && ~isempty(dim)
	order = [dim:ndims(x) 1:dim-1];
	x = permute(x, order);
	y = jf_prctile(pn, x, p);
	y = ipermute(y, order);
return
end

if ndims(x) > 2
	dim = size(x);
	x = reshapee(x, [], prod(dim(2:end)));
	y = jf_prctile(pn, x, p);
	y = reshapee(y, [], dim(2:end));
return
end

if any(isnan(x(:)))
	fail 'nan values unsupported'
end
if ~isreal(x)
	fail 'complex values unsupported'
end

% at this point x is 2d and we work along 1st dimension

n = size(x, 1);
q = [0 100*(0.5:(n-0.5))./n 100]';
x = sort(x,1);
x = [x(1,:); x(1:n,:); x(n,:)];
y = interp1q(q, x, p(:));


function jf_prctile_test
%x = zeros(1,51,1);
%x(1,:,1) = 200:-4:0;
%x(1,:,1) = 0:2:100;
x = rand(10,40,20);
p = [20 50 60];
y1 = jf_prctile([], x, p, 2);
if exist('prctile') == 2
	y2 = prctile(x, p, 2);
	jf_equal(y1, y2)
end
