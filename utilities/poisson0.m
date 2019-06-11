 function data = poisson0(xm)
%function data = poisson0(xm)
%	generate poisson random vector with mean xm
%	brute force by summing exponentials - dumb loops

	for ii=1:length(xm)
		g = exp(-xm);
		data(ii) = 0;
		t = rand(1,1);
		while t > g
			t = t * rand(1,1);
			data(ii) = data(ii) + 1;
		end
	end
