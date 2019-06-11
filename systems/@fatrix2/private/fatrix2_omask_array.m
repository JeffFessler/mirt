function out = fatrix2_omask_array(st)
if isempty(st.omask)
	out = true([st.odim 1]);
else
	out = st.omask;
end
