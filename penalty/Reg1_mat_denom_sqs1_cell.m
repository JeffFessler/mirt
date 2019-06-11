  function denom = Reg1_mat_denom_sqs1_cell(C1s, mask, pots, ws, x)
%|function denom = Reg1_mat_denom_sqs1_cell(C1s, mask, pots, ws, x)
%|
%| helper function for SQS denominator for (e.g., edge-preserving) regularizer
%| in
%|	C1s	{M}		cell of finite differencing matrices
%|	mask	[(N)]
%|	pots	{M}		cell of potential function strum
%|	ws	strum		weighting for each m
%|	x	[(N)] or [np]
%| out
%|	denom	[(N)] or [np]

[x ei] = embed_in(x, mask);
denom = 0;
for mm=1:numel(C1s)
	Cm = C1s{mm};
	d = Cm * x;
	Cm = abs(Cm);
	ck = Cm * ones(size(x)); % |C|*1
	pot = pots{mm};
	wt = pot.wpot(d); % potential function Huber curvatures
	wt = wt .* reshape(ws.col(mm), size(wt));
	tmp = Cm' * (wt .* ck);
	denom = denom + tmp;
end
if ei.column
	denom = denom(mask);
else
	denom = denom .* mask; % apply mask if needed
end
