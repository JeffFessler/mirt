 function st = nufft_table_change(st, om)
% change to a new set of frequencies

if ~isvar(st, 'phase_shift')
	error 'only done for "table" version'
end

st.om = om;
st.M = nrow(om);
st.phase_shift = exp(1i * (st.om * st.n_shift(:)));
