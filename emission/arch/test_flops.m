%	test_flops.m
%	compare ML-EM vs QPL-EMDP flops

if ~isvar('G')
	!wt -chat 0 dsc 2 na 80 >! t.dsc
	!wt gen t
	[G, n.x, n.y, n.b, n.a] = wtfmex('load', 't.wtf');

	mask = reshape(sum(G) ~= 0, n.x, n.y);
	im(mask)

	G = G(:, mask(:)~=0);
return
end


if ~isvar('C')
	C = sqrt(2^(-4)) * Csparse('maskleak', mask, 1);
return
end

	x = ones(sum(mask(:)), 1);
	yi = G * x;
	ci = ones(size(yi));
	ri = 0.001 + zeros(size(yi));


	flops(0)
	x = eml_em(x, G, yi(:), ci(:), ri(:), []);
	flops

	flops(0)
	x = eql_emdp(x, G, yi(:), ci(:), ri(:), [], C);
	flops
