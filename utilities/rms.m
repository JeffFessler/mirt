 function [rh, sh] = rms(x)
%function [rh, sh] = rms(x)
%|
%| Root-mean-squared values of vector x (typically an error vector)
%|
%| in
%|	x: [N M]	error vectors
%| out
%|	rh: [1 M]	estimated rms error
%|	sh: [1 M]	estimated std. dev. of rms estimate

if nargin < 1, ir_usage, end
if streq(x, 'test'), rms_test, return, end

N = nrow(x);	if (N==1), fail('use column vector'), end
bs = mean(x, 1, 'double');
st = std(x);
rh = sqrt(mean(abs(x).^2, 1, 'double'));

var_mse = (2*st.^4 + 4*st.^2 .* bs.^2)/N;
if rh > 0
	sh = sqrt( var_mse / 4 ./ rh.^2 );
else
	sh = 0;
	if nargout > 1
		warning 'sh meaningless'
	end
end

if ~nargout
	base = '';
	fprintf('%srms(%s) =', base, inputname(1))
	fprintf(' %g', rh)
	fprintf('\n')
	clear rh
end


function rms_test
rng(0)
err = (rand(10^4,3)-0.5) * sqrt(12);
rms(err)
