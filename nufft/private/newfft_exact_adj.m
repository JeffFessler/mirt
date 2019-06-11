 function x = newfft_exact_adj(st, X, om)
%function x = newfft_exact_adj(st, X, om)
%|
%| adjoint of exact NUFFT

nthread = 1; % todo
useloop = false; % todo

if ~isvar('om') || isempty(om)
	om = st.om;
end
if isempty(om), fail 'om or st.om required', end

if exist('dtft_mex') == 3
	x = jf_mex('dtft,adjoint', double(om'), double(X), ...
		int32(st.Nd), int32(nthread));
	if any(st.n_shift)
		fail 'n_shift not done'
	end
else
	x = dtft_adj(X, om, st.Nd, st.n_shift, useloop);
end
