function out = fatrix2_imask_array(st)
if isempty(st.imask)
	out = true([st.idim 1]);
else
	out = st.imask;
end
