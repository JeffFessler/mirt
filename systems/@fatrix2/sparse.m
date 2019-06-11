 function sp = sparse(ob)
%function sp = sparse(ob)
% return sparse version of object
% this uses less memory than full(ob)

if ~isempty(ob.handle_sparse)
	sp = ob.handle_sparse(ob);
return
end

% extract one column at time, determine number of nonzeros, then combine
% todo: use blocks of 10 or 100 at a time instead of 1?

np = size(ob,2);
val_cell = cell(np,1);
nnz_all = zeros(np,1);

z = zeros(size(ob,2), 1);
for jj=1:np
	x = z; x(jj) = 1;
	tmp = ob * x;
%	tmp = ob(:,jj); % error!?
	tmp = sparse(double(tmp));
	val_cell{jj} = tmp;
	nnz_all(jj) = nnz(tmp);
end

nnz_sum = sum(nnz_all);

rows = zeros(nnz_sum,1);
cols = zeros(nnz_sum,1);
vals = zeros(nnz_sum,1);

start = 0;
for jj=1:np
	[row col val] = find(val_cell{jj});
	ii = start + [1:nnz_all(jj)];
	rows(ii) = row;
	cols(ii) = jj;
	vals(ii) = val;
	start = start + nnz_all(jj);
end

sp = sparse(rows, cols, vals, size(ob,1), size(ob,2), nnz_sum);
