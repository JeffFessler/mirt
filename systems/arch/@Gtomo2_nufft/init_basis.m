function ob = init_basis(ob)
% initialization of the image basis function (blob/Gauss/...)
% by S. Matej

if streq(ob.basis.type, 'KB')
	dim1 = floor(ob.basis.diam/2);
	dim2 = 1+ 2*dim1;
	row = [-dim1:dim1] .* [-dim1:dim1];
	onev(1:dim2) = 1;
	kernel = row'*onev + onev'*row;
	%normft = ...
	%kaiser_bessel_ft(0,ob.basis.diam ,ob.basis.shape,ob.basis.m,ob.basis.dim)
	if ob.basis.dim==2
		kernel = sqrt(kernel);
		ob.basis.kernel = kaiser_bessel( ...
			kernel, ob.basis.diam, ob.basis.shape, ob.basis.m);

	elseif ob.basis.dim==3
		for i = 1:dim2
			kernel3(:,:,i) = kernel+row(i);
		end
		kernel3 = sqrt(kernel3);
		ob.basis.kernel = kaiser_bessel(...
			kernel3, ob.basis.diam, ob.basis.shape, ob.basis.m);

	else
		error(sprintf('basis function dimension %g not supported', ob.basis.dim))
	end

	norm=sum(sum(sum(ob.basis.kernel)));
	ob.basis.kernel = ob.basis.kernel/norm;

elseif streq(ob.basis.type,'Gauss')
	error('Gaussian basis not yet implemented - kernel not calculated')

elseif isempty(ob.basis.type) | streq(ob.basis.type,'pixel') | ...
	streq(ob.basis.type,'no')

else
	error(sprintf('basis function %s not implemented', ob.basis.type))

end
