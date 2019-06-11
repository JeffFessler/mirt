 function ob = ctranspose(ob)
%function ob = ctranspose(ob)

if ob.is.empty
	 error empty
end
	if ob.is.subref
		ob.is.transpose_after_sub = ~ob.is.transpose_after_sub;
	else
		ob.is.transpose = ~ob.is.transpose;
	end
