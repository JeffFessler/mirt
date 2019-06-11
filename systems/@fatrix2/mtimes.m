 function y = mtimes(ob, x)
%function y = mtimes(ob, x)
% A * x
% A' * y
% y' * A
% x' * A'
% B * A		(calls private/fatrix2_mtimes2)
% etc

if isnumeric(ob)
	if isscalar(ob) % scalar * object
		y = fatrix2_scalar_times(ob, x);
	else % row_vector(s) * object
		y = (x' * ob')'; % trick: x <-> ob
	end
return
end

if ~isnumeric(x) % object * object
	if islogical(x)
		warn 'fatrix2 * logical is illogical; reconsider!'
		x = single(x);
%	elseif iscell(x) % untested feature for mai le
		% do nothing here; continue to "ordinary" multiplication below
	else
		y = fatrix2_mtimes2(ob, x);
	return
	end
end

y = fatrix2_mtimes_vector(ob, x); % "ordinary" fatrix * vector multiplication

end % mtimes()


% fatrix2_mtimes_vector()
% ordinary matrix-vector multiplication, but generalized to support arrays
% full multiplication
% y = scale * diag(odiag) * ob * diag(idiag) * x;
function y = fatrix2_mtimes_vector(ob, x)

fun_forw = @(ob, x) fatrix2_do_forw(ob, x);

[nd np] = size(ob);
dimx = size(x);
if dimx(1) == np % [np (L)] vector mode
	diml = [dimx(2:end) 1];
	LL = prod(diml);
	if LL == 1 % usual single column vector case
		x = iembed(ob, x); % [idim]
		y = fun_forw(ob, x); % [odim]
		y = y(:); % [*odim 1]
		y = oselect(ob, y); % [nd 1]
	elseif ob.does_many % multiple column vector without loop
		x = reshape(x, [np LL]); % [np *L]
		x = iembed(ob, x); % [idim *L]
		y = fun_forw(ob, x); % [odim *L]
		y = reshape(y, [prod(ob.odim) LL]); % [*odim *L]
		y = oselect(ob, y); % [nd *L]
		y = reshape(y, [nd diml]); % [nd (L)]
	else % multiple column vector case with loop
		x = reshape(x, [np LL]); % [np *L]
		tmp = iembed(ob, x(:,1)); % [idim]
		tmp = fun_forw(ob, tmp); % [odim]
		y = zeros([prod(ob.odim) LL], class(tmp)); % [*odim *L]
		y(:,1) = tmp(:);
		for ll=2:LL
			tmp = iembed(ob, x(:,ll)); % [idim]
			tmp = fun_forw(ob, tmp); % [odim]
			y(:,ll) = tmp(:);
		end
		y = oselect(ob, y); % [nd *L]
		y = reshape(y, [nd diml]); % [nd (L)]
	end

elseif first_dim_match(dimx, ob.idim)

	diml = [dimx((numel(ob.idim)+1):end) 1];
	LL = prod(diml);
	if LL == 1 % a single array
		y = fun_forw(ob, x); % [odim]
	elseif ob.does_many
		x = reshape(x, [ob.idim LL]); % [idim *L]
		y = fun_forw(ob, x); % [odim *L]
		y = reshape(y, [ob.odim diml]); % [odim (L)]
	else % loop
		x = reshape(x, [prod(ob.idim) LL]); % [*idim *L]
		tmp = reshape(x(:,1), [ob.idim 1]); % [idim]
		tmp = fun_forw(ob, tmp); % [odim]
		y = zeros([prod(ob.odim) LL], class(tmp)); % [*odim *L]
		y(:,1) = tmp(:);
		for ll=2:LL
			tmp = reshape(x(:,ll), [ob.idim 1]); % [idim]
			tmp = fun_forw(ob, tmp); % [odim]
			y(:,ll) = tmp(:);
		end
		y = reshape(y, [ob.odim diml]); % [odim (L)]
	end

else
	pr size(x)
	pr ob.idim
	fail 'input dimension size mismatch'
end

end % fatrix2_mtimes_vector()


% first_dim_match()
% see if first dimensions of dimx match idim properly
% basically, dimx = [idim (L)]
function yn = first_dim_match(dimx, idim)
yn = numel(dimx) >= numel(idim) && isequal(dimx(1:numel(idim)), idim);
if ~yn && idim(end) == 1 % handle case [... 1]
	yn = first_dim_match(dimx, idim(1:end-1));
end
end % first_dim_match()
