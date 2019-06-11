 function ob = ctranspose(ob)

if ob.is.empty
	error empty
end
ob.is.transpose = ~ob.is.transpose;
