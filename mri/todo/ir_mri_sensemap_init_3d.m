 function sinit = ir_mri_sensemap_init_3D(init, ykj, bodycoil, sizeI, thresh, maskObj)

%no complicated initialization like fit 1st order

nx = sizeI(1);
ny = sizeI(2);
nz = sizeI(3);

if ~isempty(maskObj)
	warning('using the maskObj instead of thresholding to determine good pixels');
end

if isempty(init) || ischar(init)
	sinit = zeros(nx, ny, nz);
	tmp = ykj ./ bodycoil; % usual ratio

	ind = isnan(tmp);
	tmp(ind) = 1;

	if isempty(maskObj)
		good = abs(bodycoil) > thresh * max(abs(bodycoil(:)));
	else
		good = maskObj;
	end

	if isempty(init) || strcmp(init, 'ratio')
		% set all uncertain map values to median of good ones
		disp('Using zero background.');
		tmp(~good) = 0; % median(abs(tmp(good)));;
		tmp = reshape(tmp, [nx,ny,nz]);
	elseif isempty(init) || strcmp(init, 'median')
		% set all uncertain map values to median of good ones
		disp('Using median background.');
		tmp(~good) = median(tmp(good));
		tmp = reshape(tmp, [nx,ny,nz]);
	elseif isempty(init) || strcmp(init, 'avg')
		% set all uncertain map values to median of good ones
		disp('Using mean background.');
		tmp(~good) = mean(abs(tmp(good))) .* exp(1i*mean(angle(tmp(good))));
		tmp = reshape(tmp, [nx,ny,nz]);
	elseif strcmp(init, 'zeros')
		tmp = zeros([nx ny,nz]);

	else
		error('Do not recognize initialization type.')
	end

	sinit(:,:,:) = reshape(tmp, [nx ny nz]);

else
	sinit = init;
end

end %
