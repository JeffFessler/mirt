 function y = fatrix2_do_forw(ob, x)
%function y = fatrix2_do_forw(A, x)
%|
%| implement 'forw' function: A * x
%| in array mode, accounting for idiag, odiag, scale

x = fatrix2_apply_diag(x, ob.idiag); % [idim *L]
y = ob.handle_forw(ob.arg, x); % [odim *L]
y = fatrix2_apply_diag(y, ob.odiag); % [odim *L]

if ~isequal(ob.scale, 1)
	y = ob.scale * y;
end

end % fatrix2_do_forw()


% fatrix2_apply_diag()
function z = fatrix2_apply_diag(z, adiag)
if ~isempty(adiag)
%	z = adiag .* z; % relies on singleton expansion of 2016b
	z = bsxfun(@times, adiag, z);
end

end % fatrix2_apply_diag()
