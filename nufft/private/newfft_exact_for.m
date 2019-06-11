 function X = newfft_exact_for(st, x, om)
%function X = newfft_exact_for(st, x, om)
%|
%| exact forward NUFFT

nthread = 1; % todo
how = 'outer'; % todo

if ~isvar('om') || isempty(om)
	om = st.om;
end
if isempty(om), fail 'om or st.om required', end

if exist('dtft_mex') == 3
	X = jf_mex('dtft,forward', double(om'), double(x), int32(nthread));
	if any(st.n_shift)
		fail 'n_shift not done'
	end
else
	X = dtft(x, om, 'n_shift', st.n_shift, 'how', how);
end
