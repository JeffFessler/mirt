 function full_data = ir_mri_dyn_data_share(undersamp_data, sampling_pattern)
%function full_data = ir_mri_dyn_data_share(undersamp_data, sampling_pattern)
%|
%| A data-sharing approach to imputing missing k-space data in dynamic MRI.
%| Takes undersampled dynamic data and uses data-sharing (0th order interp)
%| to fill in gaps.
%| Output k-space data is "fully sampled"
%| Data-sharing fills in unsampled k-space locations with values
%| from nearest (in time) sampling of that location, using mean for ties.
%| Leaves zeros if k-space location not ever sampled.
%| Currently works only for 2D+T.
%|
%| in
%|	sampling_pattern [Nx Ny Nf] logical array
%|	undersamp_data [Ns Nc] where Ns = a*Nx*Ny*Nf and a is undersamp factor
%|			or [Nx Ny Nf Nc] with zeros where data is not sampled
%|
%| out
%|	full_data [Nx Ny Nf Nc] complex array
%|
%| 02/05/14 Mai Le, University of Michigan
%| 2014-09-07 tweaks by Jeff Fessler

if nargin == 1 && streq(undersamp_data, 'test')
	ir_mri_dyn_data_share_test
return
end

if nargin < 2, help(mfilename), error(mfilename), end

if ~islogical(sampling_pattern)
	fail 'sampling_pattern must be logical'
end

Nx = size(sampling_pattern,1);
Ny = size(sampling_pattern,2);
Nf = size(sampling_pattern,3);

if (size(undersamp_data,1) == Nx ...
&& size(undersamp_data,2) == Ny ...
&& size(undersamp_data,3) == Nf)
	undersamp_data = reshape(undersamp_data, Nx*Ny*Nf, []);
	undersamp_data = undersamp_data(sampling_pattern(:), :);
elseif ndims(undersamp_data) > 2
	fail('bad data size')
end

[Ns Nc] = size(undersamp_data);

flip_samp = flipdim(sampling_pattern,3);

samps = {sampling_pattern, flip_samp};

zfill = reshape(embed(undersamp_data, sampling_pattern), Nx,Ny,Nf,Nc);

for ii = 1:2 % forward and backward passes
	samp = samps{ii};
	interp_ind = cumsum(int32(samp),3);
	interp_ind_pos = interp_ind + int32(interp_ind == 0);
	siip(:,:,:,ii) = interp_ind_pos; % keep!!

	full_interp_ind_pos = Nx*Ny*(interp_ind_pos(:)-1) + int32(repmat((0:Nx*Ny-1)', [Nf 1])) + 1;

	if Nc > 1
		full_interp_ind_pos = repmat(full_interp_ind_pos,[Nc 1]) + int32(kron((0:Nc-1)',Nx*Ny*Nf*ones(Nx*Ny*Nf,1)));	
	end

	[~, inds] = sort(samp, 3, 'descend');
	sinds(:,:,:,ii) = inds; % keep!!

	full_inds = Nx*Ny*(inds(:)-1) + repmat((1:Nx*Ny)',[Nf 1]);
	if Nc > 1
		full_inds = repmat(full_inds,[Nc 1]) ...
			+ kron((0:Nc-1)', Nx*Ny*Nf*ones(Nx*Ny*Nf,1));
	end
	
	if ii == 2
		curr_zfill = col(flipdim(zfill,3));
	else
		curr_zfill = col(zfill);
	end
	nonzero_elms(:,ii) = curr_zfill(full_inds);
	
	zero_less_fill = reshape(nonzero_elms(:,ii),Nx,Ny,Nf,Nc);
	curr_full_data = reshape(zero_less_fill(full_interp_ind_pos),Nx,Ny,Nf,Nc);

	if ii == 2 
		curr_full_data = flipdim(curr_full_data, 3);
	end

	sfull_data(:,:,:,:,ii) = curr_full_data;
end

siip1 = double(siip(:,:,:,1));
c1_ndx = (siip1-1)*Nx*Ny+repmat(reshape(0:Nx*Ny-1,Nx,Ny),[1 1 Nf])+1;
frames1 = col(sinds(:,:,:,1));
c1 = col(reshape(frames1(c1_ndx),Nx,Ny,Nf));

siip2 = double(flipdim(siip(:,:,:,2),3));
siip2 = repmat(siip2(:,:,1) + 1,[1 1 Nf]) - siip2;
c2_ndx = (siip2-1)*Nx*Ny+repmat(reshape(0:Nx*Ny-1,Nx,Ny),[1 1 Nf])+1;
frames2 = col(sinds(:,:,:,1)); % need to have 1!
c2 = col(reshape(frames2(c2_ndx),Nx,Ny,Nf));

c3 = col(reshape(kron((1:Nf),ones(Nx*Ny,1)),Nx,Ny,Nf));
%[c1 c2 c3 col(1:Nx*Ny*Nf)]

choose_flip = repmat(reshape(abs(c2 - c3) < abs(c1 - c3),Nx,Ny,Nf),[1 1 1 Nc]);

tie = repmat(reshape((abs(c2 - c3) == abs(c1 - c3)) & (c1 ~= c2), Nx,Ny,Nf), [1 1 1 Nc]); 

flip_weights = choose_flip + 0.5*tie;
full_data = sfull_data(:,:,:,:,1) .* (1-flip_weights) ...
		+ sfull_data(:,:,:,:,2) .* flip_weights;


function ir_mri_dyn_data_share_test
down = 2;
% n_tr_merge = 100;
n_tr_merge = 112; % 14 * 16 / 2
nx = 192 / down;
ny = 216 / down;
Nf = 12;
n_low = [14 16];
% nx * ny + prod(n_low) * Nf
% Nt = round(480 / down^0);
Nt = 20 * prod(n_low) * Nf / n_tr_merge;
[samp1 samp2] = ir_mri_dce_samp1([nx ny Nt], 'n_tr_merge', n_tr_merge, ...
	'n_low', n_low, 'Nframe', Nf, 'chat', 0);

im plc 2 2
im(1, 'row', 3, samp1)
xtrue = repmat(1:Nf, [nx*ny 1]);
xtrue = reshape(xtrue, [nx ny Nf]);

im(2, 'row', 3, xtrue)

xsamp = xtrue .* samp1;
im(3, 'row', 3, xsamp)
drawnow

xshare = ir_mri_dyn_data_share(xsamp, samp1);
im(4, 'row', 3, xshare)
