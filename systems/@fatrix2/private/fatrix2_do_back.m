 function x = fatrix2_do_back(ob, y)
%function x = fatrix2_do_back(A, y)
%|
%| implement 'back' function: A' * x
%| in array mode, accounting for idiag, odiag, scale
%|
%| trick: conj() of idiag odiag scale already done in ctranspose.m

y = fatrix2_apply_diag(y, ob.odiag); % [odim *L]
x = ob.handle_back(ob.arg, y); % [idim *L]
x = fatrix2_apply_diag(x, ob.idiag); % [idim *L]

if ~isequal(ob.scale, 1)
	x = ob.scale * x;
end

end % fatrix2_do_back()


% fatrix2_apply_diag()
function z = fatrix2_apply_diag(z, adiag)
if ~isempty(adiag)
%	z = adiag .* z; % relies on singleton expansion of 2016b
	z = bsxfun(@times, adiag, z);
end

end % fatrix2_apply_diag()
