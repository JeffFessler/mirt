  function y = penalty_mex_call(diff_str, x, offsets, ndim)
%|function y = penalty_mex_call(diff_str, x, offsets, ndim)
%| convenience routine for calling penalty_mex
%| to handle both real and complex case

if nargin < 1, help(mfilename), error(mfilename), end

offsets = int32(offsets);
ndim = int32(ndim);

if ~isreal(x)
	yr = penalty_mex(diff_str, single(real(x)), offsets, ndim);
	yi = penalty_mex(diff_str, single(imag(x)), offsets, ndim);
	y = complex(yr, yi);
else
	y = penalty_mex(diff_str, single(x), offsets, ndim);
end
