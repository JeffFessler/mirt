 function scale = ir_best_scale(x, y)
%function scale = ir_best_scale(x, y)
%| scale factor for x so it best matches y
%| min_scale || y - x * scale ||
x = x(:);
y = y(:);
scale = sum(col(conj(x) .* y), 'double') / sum(col(abs(x).^2), 'double');
