 function trl_check(yi, bi, ri)
%function trl_check(yi, bi, ri)
% check the validity of transmission data from Poisson transmission model

if any(yi(:) < 0), warning 'negative yi values!?', end

if ~isempty(ri)
	if ~dims_same(yi, ri, 'scalar_ok', 1), error 'yi vs ri size', end
	if any(ri(:) < 0), error 'negative ri value?', end
else
	ri = 0;
end

if ~isempty(bi)
	if ~dims_same(yi, bi, 'scalar_ok', 1), error 'bi vs ri size', end
	if any(bi(:) < 0), error 'negative bi value?', end
else
	bi = 1;
end

if any(yi(:) & ~(bi(:) + ri(:))), error 'Poisson model mismatch', end
