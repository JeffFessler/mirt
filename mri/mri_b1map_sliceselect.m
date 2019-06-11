 function [zmaps omaps cost] = mri_b1map_sliceselect(yy, alpha, varargin)
%function [zmaps omaps cost] = mri_b1map_sliceselect(yy, alpha, [options])
%|
%| Estimate "B1+ map" for each of ncoil coils
%| from sequence of reconstructed images with ntip different nominal tip angles.
%| Model:
%|	todo
%| in
%|	yy	[nx ny nmeasure]	noisy complex images
%|	alpha	[nmeasure ncoil]	indicates relative tip angles as well as
%|					which coils are operating
%|
%| options
%|	profile		'real_rect1', 'gauss', etc (see load_tables)...
%|	l2b_zmap	log_2(beta_zmap), regularization parameter for b1map
%|	l2b_zmap_1pass	log_2(beta_zmap) for first pass
%|	order		regularization order (default: 2)
%|	niter		# of iterations
%|	niter_1pass	# of iterations for first pass
%|	init_zmap	initial b1map for iterations
%|	init_omap	initial object magnetization map (should be real)
%|	updateo_length	how often to update object - 1: every time, 2: every other time, etc.
%|	isave		which iterations to save. (default: last)
%|	userfun		user defined function handle (see default below)
%|	userarg		arguments to user defined function
%|	kappa		1 to use kappa and modified penalty (for zmap), (default: 1)
%|	scale		adjust by median magnitude of object, (default: 1)
%| out
%|	zmaps	[nx ny ncoil niter]	regularized b1map estimate(s)
%|	omaps	[nx ny ncoil niter]	unregularized object estimate(s)
%|	cost	[niter]			value of cost function
%|
%|
%| Because this method smoothly interpolates over regions with signal voids,
%| it is an improvement over the usual "double angle" approach to B1 mapping.
%| Described in ISBI 2007 (p616), 2008 ISMRM (p3415), 2009 ISMRM (p2609).
%|
%| Copyright 2008-7, Amanda Funai, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if streq(yy, 'test'), mri_b1map_test, return, end
if streq(yy, 'spatial'), mri_b1map_test_spatial, return, end
if streq(yy, 'spatial2'), mri_b1map_test_spatial2, return, end
if streq(yy, 'test_less'), mri_b1map_test_less, return, end
if streq(yy, 'function'), mri_b1map_function_test, return, end

% defaults
arg.profile = 'real_rect1';
%arg.profile = 'short_trunc_sinc';
arg.l2b_zmap = -3;
arg.l2b_zmap_1pass = -10;
arg.order = 2;
arg.niter = 50;
arg.niter_1pass = 5;
arg.init_zmap = [];
arg.init_omap = [];
arg.updateo_length = 10;
arg.isave = [];
arg.userfun = @userfun_default;
arg.userarg = {};
arg.kappa = 1;
arg.scale = 1;
arg.rotate = 0;
arg.rotatedir = 1;

arg = vararg_pair(arg, varargin);
if isempty(arg.isave), arg.isave = arg.niter; end

[nx ny nmeasure] = size(yy);
yy = reshape(yy, [nx*ny nmeasure]); % could generalize for 3D later

[nmeasure2 ncoil] = size(alpha);
if (nmeasure ~= nmeasure2), error 'size of yy and alpha not consistent', end

%% run through first pass

[z1pass o1pass cost1pass] = mri_b1map_iterate(yy, alpha, arg, nx, ny, ...
	ncoil, 1, ones(nx,ny,ncoil), ones(nx,ny), ...
	arg.l2b_zmap_1pass, arg.niter_1pass, 1);


%% NOTE - CHANGED THE SCALING!!!!!
if(arg.scale)
	%scale_factor = median(o1pass(o1pass(:) > .01 * max(o1pass(:))))
	scale_factor = median(o1pass(o1pass(:) > .1 * max(o1pass(:))));
%	o_filt = medfilt2(o1pass,[2 2]);
%	scale_factor = median(o_filt(o_filt(:) > .01 * max(o_filt(:))))
%	arg.l2b_zmap = arg.l2b_zmap + 2*log2(scale_factor);
%scale_factor = median(o1pass(o1pass(:) > .001 * max(o1pass(:))))
	yy = yy / scale_factor;
end

%keyboard

%% run through main algorithm
[zmaps omaps cost] = mri_b1map_iterate(yy, alpha, arg, nx, ny, ncoil, 0, ...
	z1pass, o1pass, arg.l2b_zmap, arg.niter,scale_factor);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mri_b1map_kappa_calc
%
function [mask1 mask2] = mri_b1map_kappa_calc(alpha, arg, nx, ny, ncoil, first_pass, z1pass, o1pass)

z1pass = reshape(z1pass,[nx ny ncoil]);
o1pass = reshape(o1pass,[nx ny]);

if(first_pass)
	mask1 = logical(ones(nx,ny,ncoil));
	mask2 = logical(ones(nx,ny,ncoil));
else
	[h h_deriv h_deriv2] = load_tables(arg.profile);

	kmask_1 = o1pass > .1 * max(o1pass(:));
	%kmask_1 = ones(size(o1pass));
	for i=1:ncoil
		kmask(:,:,i) = kmask_1;
	end

	if(arg.kappa)
		z1pass = reshape(z1pass,[nx*ny ncoil]);
		mask1 = zeros(nx,ny,ncoil);
		mask2 = zeros(nx,ny,ncoil);
	%	cross_terms = zeros(nx,ny,ncoil);
		sum_factors = zeros(ncoil,1);

		for icoil=1:size(alpha,2)
		for im=1:size(alpha,1)
			%[F_prime] = F_2_prime( alpha(im,:) * z1pass', h_deriv)';
			%F_prime = reshape(F_prime,[nx ny]);
			%F2 = F_2( alpha(im,:) * z1pass', h)';
			%F2 = reshape(F2,[nx ny]);
			%F2R = real(F2);
			%F2I = imag(F2);

			%z = alpha(im,:) * z1pass';
			%z = reshape(z, [nx ny]);

			[da_F2R db_F2R da_F2I db_F2I] = ...
				F_2RI_prime2(alpha(im,:) * real(z1pass)', ...
       			        alpha(im,:) * imag(z1pass)', h_deriv, h);

			da_F2R = reshape(da_F2R,[nx ny]);
			db_F2R = reshape(db_F2R,[nx ny]);
			da_F2I = reshape(da_F2I,[nx ny]);
			db_F2I = reshape(db_F2I,[nx ny]);

			mask1(:,:,icoil) = mask1(:,:,icoil) + alpha(im,icoil).^2 * ((da_F2R).^2 + (da_F2I).^2);
			mask2(:,:,icoil) = mask2(:,:,icoil) + alpha(im,icoil).^2 * ((db_F2R).^2 + (db_F2I).^2);
		%	cross_terms(:,:,icoil) = cross_terms(:,:,icoil) ...
		%	+ alpha(im,icoil).^2 * (da_F2R .* db_F2R + da_F2I .* db_F2I);
		%	temp = mask1;
		%	temp2 = mask2;
		%	sum_factors(icoil) = sum_factors(icoil) + alpha(im,icoil).^2;
		end

		%mask1(:,:,icoil) = mask1(:,:,icoil) ./ sum_factors(icoil);
		%mask2(:,:,icoil) = mask2(:,:,icoil) ./ sum_factors(icoil);
		end

		mask1 = mask1 .* kmask;
		mask1(mask1 == 0) = mean(mask1(:) ~= 0);
		mask1 = sqrt(mask1);

		mask2 = mask2 .* kmask;
		mask2(mask2 == 0) = mean(mask2(:) ~= 0);
		mask2 = sqrt(mask2);

	else
		mask1 = logical(ones(nx,ny,ncoil));
		mask2 = logical(ones(nx,ny,ncoil));
	end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% mri_b1map_init()
%
function [zmap, omap] = mri_b1map_init(yy, alpha, arg, h, nx, ny);

if isempty(arg.init_zmap)
    isinitial_zmap = 0;

[nmeasure ncoil] = size(alpha);
[np nmeasure] = size(yy);

found = 0;
for p = 1:nmeasure
  for q = 1:nmeasure
  if(sum(alpha(q,:) == 2 * alpha(p,:)) == ncoil)
    found = found+1;
    found_matrix(found,:) = [p q];
  end
  end
end

if(found > 0)
    % Quick fix if we found too many double angle formulas
    % e.g. if we use tips of [1 2 4] for instance
    if(found > ncoil)
       found = ncoil;
    end
    for u=1:found
    num = yy(:,found_matrix(u,2));
    den = yy(:,found_matrix(u,1));

    alpha_sum(u,:) = alpha(found_matrix(u,1),:);

    %% testing - find an init for the abs and the phase
    absz_estimate = acos(num ./ (2*den));
    absz_estimate = real(absz_estimate);

    % trying something else
    %absz_estimate = acos(abs(num ./ (2*den)));
    %absz_estimate = abs(absz_estimate);

    absz_estimate(isnan(absz_estimate(:))) = 0;

    F_1_absz = absz_estimate .* h_bloch(absz_estimate,h);
    F_1_2absz = 2 * absz_estimate .* h_bloch(absz_estimate*2,h);

    anglez_estimate_temp1 = angle(den) - angle(F_1_absz);
    anglez_estimate_temp2 = angle(num) - angle(F_1_2absz);
    [junk max_loc] = max([F_1_absz F_1_2absz],[],2);

    anglez_estimate = size(anglez_estimate_temp1);
    anglez_estimate(max_loc==1) = anglez_estimate_temp1(max_loc == 1);
    anglez_estimate(max_loc==2) = anglez_estimate_temp2(max_loc == 2);

    anglez_sum(u,:) = anglez_estimate;
    absz_estimate_sum(u,:) = absz_estimate;

    end

if (found < ncoil && 0)
finding = found;
for u=1:ncoil
        if(isempty(find(found_matrix(:) == u)))
        finding = finding+1;
        alpha_sum(finding,:) = alpha(u,:);
        if(arg.rotate)
            rotation_angle = (u-1)/ncoil * 360 * arg.rotatedir;
            a = reshape(mean(anglez_sum(1:found,:),1),[nx ny]);
            b = reshape(mean(absz_estimate_sum(1:found,:),1),[nx ny]);
            temp1 = imrotate(a,rotation_angle);
            anglez_sum(finding,:) = reshape(imresize(temp1,[nx ny]),[nx*ny 1]);
            temp2 = imrotate(b,rotation_angle);
            absz_estimate_sum(finding,:) = reshape(imresize(temp2,[nx ny]),[nx*ny 1]);
        else
            anglez_sum(finding,:) = mean(anglez_sum(1:found,:),1);
            absz_estimate_sum(finding,:) = mean(absz_estimate_sum(1:found,:),1);
        end
    end
end
end


if (found < ncoil && 0)
% estimate omap based on rows with DAM
omap_temp = mri_b1map_omap(yy(:,found_matrix(:)),...
                [eye(size(found_matrix,1)); 2*eye(size(found_matrix,1))], ...
                real(absz_estimate_sum.* exp(1i * anglez_sum))', ...
                imag(absz_estimate_sum.*exp(1i*anglez_sum))',h);
omap_temp_mask = omap_temp > .1*max(omap_temp(:));
% above may need a squeeze depending on order of rows
% need to test later

finding = 0;
%finding = found;
for u=1:ncoil
    %if(isempty(find(found_matrix(:) == u)))
        finding = finding+1;
        alpha_sum(finding,:) = alpha(u,:);
        absz_estimate_sum(finding,:) = real(asin(abs(yy(:,u)) ./ omap_temp ));
        absz_estimate_sum(finding,omap_temp<.1*max(omap_temp(:))) = 0;
        anglez_sum(finding,:) = mean(anglez_sum(1:found,:),1);
    %end
end

%keyboard

end

if (found < ncoil && 0)
omap_temp = mri_b1map_omap(yy(:,found_matrix(:)),...
                [eye(size(found_matrix,1)); 2*eye(size(found_matrix,1))], ...
                real(absz_estimate_sum.* exp(1i * anglez_sum))', ...
                imag(absz_estimate_sum.*exp(1i*anglez_sum))',h);
omap_temp_mask = logical(omap_temp > .1*max(omap_temp(:)));

% Note - this only works for 4 maps with 90 degree rotation
%nx=12;
%ny=12;
b = absz_estimate_sum(1,:);
[xx yy] = meshgrid(-nx/2+.5:1:nx/2-.5, -ny/2+.5:1:ny/2-.5);
A = zeros(nx*ny);
for i=1:nx*ny
if(omap_temp_mask(i))

%% Fill in this row of the A matrix
xprime = xx(i)
yprime = yy(i)

   index1 = (yprime + ceil(ny/2)-.5)*nx + xprime + ceil(nx/2)+.5
   index4 = (-xprime + ceil(ny/2)-.5)*nx + yprime + ceil(nx/2)+.5
   index3 = (-yprime + ceil(ny/2)-.5)*nx - xprime + ceil(nx/2)+.5
   index2 = (xprime + ceil(ny/2)-.5)*nx - yprime + ceil(nx/2)+.5
   A(i,index1) = alpha(1,1);
   A(i,index2) = alpha(1,2);
   A(i,index3) = alpha(1,3);
   A(i,index4) = alpha(1,4);

%keyboard
end
end

zmap = A\b;

end
%% NOTE - this is what I used in paper
if (found < ncoil && 0)

% Estimate object based on data
omap_temp = mri_b1map_omap(yy(:,found_matrix(:)),...
                [eye(size(found_matrix,1)); 2*eye(size(found_matrix,1))], ...
                real(absz_estimate_sum.* exp(1i * anglez_sum))', ...
                imag(absz_estimate_sum.*exp(1i*anglez_sum))',h);

% Mask out the lowest pixels
% Enlarge the mask to include the whole brain (eg, the sinuses)
omap_temp_mask = logical(omap_temp > .1*max(omap_temp(:)));
mask_temp = imdilate(omap_temp_mask,strel('disk',5));
omap_temp_mask = imerode(mask_temp,strel('disk',5));
omap_temp_mask = reshape(omap_temp_mask,[nx ny]);

% Note - from below we are assuming that
% we have double angle from measurement 1

b = absz_estimate_sum(1,:);
c = anglez_sum(1,:);

smallA = zeros(ncoil);
for i=1:ncoil
    smallA(i,:) = circshift(fliplr(alpha(1,:)),[1 i]);
end

bigz = zeros(nx*ny,ncoil);
bigz_angle = zeros(nx*ny,ncoil);

nmin = min(nx,ny);
diff = abs(nx-ny)/2;
if(nx > ny)
    x = diff;
else
    x = 0;
end

for xprime=-nmin/2+.5:1:-.5
    x = x+1;
    if(ny > nx)
        y = diff;
    else
        y = 0;
    end
    for yprime = -nmin/2+.5:1:-.5
    y = y+1;

for i=1:ncoil
% Probably need some interpolation in here for later
% Right now works fine for 4 coils
    theta = (i-1)* 360 / ncoil;
    rot_matrix_clock = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
    new_coords = rot_matrix_clock^-1 * [xprime; yprime];
    index(i) = (new_coords(2) + ceil(ny/2)-.5)*nx + (new_coords(1) + ceil(nx/2) + .5);
end

   %if 1
   if(sum(omap_temp_mask(index)) == 4)

% Let's add in some "whitening" or weighting based on object magnetization

   smallb = b(index);
   smallc = c(index);
   smallweights = diag(omap_temp(index));

   smallA_weighted = smallweights^(.5) * smallA;
   smallb_weighted = smallweights^(.5) * smallb';

   %smallz = smallA \ smallb';
   smallz = smallA_weighted \ smallb_weighted;

   else
       smallz = [0 0 0 0];
       smallc = [0 0 0 0];
   end

   for i=1:ncoil
      bigz(circshift(index,[1 -i+1]),i) = smallz;
      bigz_angle(circshift(index,[1 -i+1]),i) = smallc;
   end

   if(sum(omap_temp_mask(index)) > 0 && sum(omap_temp_mask(index)) < 4 && 1)
      for i=1:ncoil
	   bigz(index,i) = b(index)/sum(alpha(1,:));
	   bigz_angle(index,i) = c(index);
      end
   end

end
end

bigz = reshape(bigz, [nx ny 4]);
omap_map_ext = repmat(reshape(omap_temp_mask,[nx ny]),[1 1 ncoil]);
bigz = bigz .* omap_map_ext;
bigz(bigz(:) == 0) = mean(bigz(bigz(:) > 0));

bigz_angle = reshape(bigz_angle, [nx ny 4]);
bigz_angle = bigz_angle .* omap_map_ext;
bigz_angle(bigz_angle(:) == 0) = mean(bigz_angle(bigz_angle(:) > 0));

figure
im pl 1 2
im(1,bigz,[0 .84]),cbar
im(2,bigz_angle,[.803 1.56]),cbar

%keyboard

zmap = bigz .* exp(1i * bigz_angle);
zmap = reshape(zmap,[nx*ny 4]);

end


%% TESTING
%% Okay - works but let's try something else
%% Let's try to estimate both abs and angle at same time
if (found < ncoil && 1)

% Estimate object based on data
omap_temp = mri_b1map_omap(yy(:,found_matrix(:)),...
                [eye(size(found_matrix,1)); 2*eye(size(found_matrix,1))], ...
                real(absz_estimate_sum.* exp(1i * anglez_sum))', ...
                imag(absz_estimate_sum.*exp(1i*anglez_sum))',h);

% Mask out the lowest pixels
% Enlarge the mask to include the whole brain (eg, the sinuses)
omap_temp_mask = logical(omap_temp > .1*max(omap_temp(:)));
mask_temp = imdilate(omap_temp_mask,strel('disk',5));
omap_temp_mask = imerode(mask_temp,strel('disk',5));
omap_temp_mask = reshape(omap_temp_mask,[nx ny]);

% Note - from below we are assuming that
% we have double angle from measurement 1

b = absz_estimate_sum;
c = anglez_sum;

% Num of DAM measures
dam_meas = size(b,1);

smallA = zeros(ncoil*dam_meas,ncoil);
for j=1:dam_meas
    for i=1:ncoil
        smallA(ncoil*(j-1)+i,:) = circshift(fliplr(alpha(j,:)),[1 i]);
    end
end

bigz = zeros(nx*ny,ncoil);
bigz_angle = zeros(nx*ny,ncoil);
bigz_total = zeros(nx*ny,ncoil);

nmin = min(nx,ny);
diff = abs(nx-ny)/2;
if(nx > ny)
    x = diff;
else
    x = 0;
end

for xprime=-nmin/2+.5:1:-.5
    x = x+1;
    if(ny > nx)
        y = diff;
    else
        y = 0;
    end
    for yprime = -nmin/2+.5:1:-.5
    y = y+1;

for i=1:ncoil
% Probably need some interpolation in here for later
% Right now works fine for 4 coils
    theta = (i-1)* 360 / ncoil;
    rot_matrix_clock = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
    new_coords = rot_matrix_clock^-1 * [xprime; yprime];
    index(i) = (new_coords(2) + ceil(ny/2)-.5)*nx + (new_coords(1) + ceil(nx/2) + .5);
end

   %if 1
   if(sum(omap_temp_mask(index)) > 0)

% Let's add in some "whitening" or weighting based on object magnetization

   smallb = b(:,index)';
   smallc = c(:,index)';
   smallbc = smallb .* exp(1i * smallc);

   smallweights = diag(repmat(omap_temp(index),[dam_meas 1]));

   smallA_weighted = smallweights^(.5) * smallA;
   smallbc_weighted = smallweights^(.5) * smallbc(:);

   %smallz = smallA \ smallb(:);
   %smallz = smallA_weighted \ smallb_weighted;
   smallz = smallA_weighted \ smallbc_weighted;

   else
       smallz = [0 0 0 0];
       smallc = [0 0 0 0];
   end

   for i=1:ncoil
      bigz_total(circshift(index,[1 -i+1]),i) = smallz;
   end

   if(sum(omap_temp_mask(index)) > 0 && sum(omap_temp_mask(index)) < 4 && 0)
      for i=1:ncoil
          %keyboard
	   bigz_total(circshift(index,[1 -i+1]),i) = mean(b(:,index))'.* exp(1i * mean(c(1,index))') /sum(alpha(1,:));
      end
   end

end
end

bigz_total = reshape(bigz_total, [nx ny 4]);
omap_map_ext = repmat(reshape(omap_temp_mask,[nx ny]),[1 1 ncoil]);
bigz_total = bigz_total .* omap_map_ext;
bigz_total(bigz_total(:) == 0) = mean(bigz_total(bigz_total(:) > 0));

zmap = reshape(bigz_total,[nx*ny 4]);

end

if(~exist('zmap'))
zmap_temp = absz_estimate_sum' .* exp(1i * anglez_sum');
zmap = alpha_sum^-1 * zmap_temp';
zmap = zmap';
end

else
   error 'b1map initialization requires double angle'
   zmap = ones(np,ncoil);
end

else
   zmap = arg.init_zmap;
   isinitial_zmap = 1;
end

if isempty(arg.init_omap)

omap = mri_b1map_omap(yy,alpha,real(zmap),imag(zmap),h);
%omap = abs(omap);

else
	omap = arg.init_omap;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% mri_b1map_iterate
%
function [zmaps omaps cost] = mri_b1map_iterate(yy, alpha, arg, nx, ny, ncoil, first_pass, z1pass, o1pass, ...
        l2beta, niter,scale_factor)

       cost = 0;
if(~first_pass)
	isave = arg.isave;
else
	isave = niter;
end

[nmeas ncoil] = size(alpha);

% load tables
[h h_deriv h_deriv2] = load_tables(arg.profile);
[c y_range] = load_denom_table(arg.profile);

% initialize all the maps
[zmap omap] = mri_b1map_init(yy, alpha, arg, h, nx, ny);

if ~isempty(arg.init_omap)
   omap = omap / scale_factor;
end

%keyboard

amap = reshape(real(zmap),[nx*ny ncoil]);
bmap = reshape(imag(zmap),[nx*ny ncoil]);

yR = real(yy);
yI = imag(yy);

% archive
zmaps = zeros(nx*ny,ncoil, length(isave));
omaps = zeros(nx*ny, length(isave));
if any(isave == 0)
    zmaps(:,:, isave == 0) = zmap;
	omaps(:, isave == 0) = omap;
end


% Regularizers
[mask1 mask2] = mri_b1map_kappa_calc(alpha, arg, nx, ny, ncoil, first_pass, z1pass, o1pass);

for i=1:ncoil
    Rzmap_R(:,i) = Robject(mask1(:,:,i), 'beta', 2^l2beta, 'order', arg.order, ...
        'type_denom', 'matlab');
    Rzmap_I(:,i) = Robject(mask2(:,:,i), 'beta', 2^l2beta, 'order', arg.order, ...
        'type_denom', 'matlab');
end

%
% iterations
%
ticker reset
iter = 1;
for iter=1:niter
    ticker(mfilename, iter, niter)

    amap_old = amap;
    bmap_old = bmap;

    [agrad bgrad] = mri_b1map_kgrad_zmap3(yR, yI, alpha, amap, bmap, omap, h, h_deriv);

    denom_value = get_denom(yy,omap,arg,alpha,ncoil,nmeas,nx*ny);

        for i=1:ncoil
            rgrad_a(:,i) = Rzmap_R(:,i).cgrad(Rzmap_R(:,i), amap(:,i));
            rgrad_b(:,i) = Rzmap_I(:,i).cgrad(Rzmap_I(:,i), bmap(:,i));

            rdenom_a(:,i) = Rzmap_R(:,i).denom(Rzmap_R(:,i),amap(:,i));
            rdenom_b(:,i) = Rzmap_I(:,i).denom(Rzmap_I(:,i),bmap(:,i));
        end

        amap = amap - (agrad + rgrad_a) ./ (denom_value + rdenom_a);
        bmap = bmap - (bgrad + rgrad_b) ./ (denom_value + rdenom_b);

		%plot_function_and_surrogate(yy,alpha,amap,bmap,omap,agrad,bgrad,...
            %rgrad_a,rgrad_b,Rzmap_R,Rzmap_I,h,denom_value,rdenom_a,rdenom_b,arg,1936,2);

        if(~rem(iter,arg.updateo_length))
            omap = mri_b1map_omap(yy, alpha, amap, bmap, h);

        end

	if any(isave == iter)
		zmaps(:,:, isave == iter) = amap+1i*bmap;
		omaps(:,isave == iter) = omap;
    end

    % Calculate the cost
    cost(iter) = mri_b1map_cost(yy,alpha,amap,bmap,omap,h,Rzmap_R,Rzmap_I);
    if(iter > 1)
        if(cost(iter)>cost(iter-1))
            plot(1:iter,cost)
        end
    end

    if ~isempty(arg.userfun), arg.userfun(arg.userarg{:}), end


%    if(~rem(iter,10)), keyboard, end

end

zmaps = reshape(zmaps, nx, ny, ncoil, []);
omaps = reshape(omaps, nx, ny, []);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pi_k = alpha_pi_k_calc(alpha,ncoil,nmeas)
for imeas = 1:nmeas
    alpha_m(imeas) = sum(abs(alpha(imeas,:)));
end
pi_k = abs(alpha) ./ repmat(alpha_m',1,ncoil);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function denom_value_new = get_denom(yy,omap,arg,alpha,ncoil,nmeas,np)

[c y_range] = load_denom_table(arg.profile);

denom_value = repmat(abs(omap).^2,1,nmeas) .* c_bound(abs(yy ./ repmat(omap,1,nmeas)),c,y_range);
denom_value(isnan(denom_value)) = c_bound(abs(yy(isnan(denom_value))),c,y_range);
pi_k = alpha_pi_k_calc(alpha,ncoil,nmeas);
denom_value_new = zeros(np,ncoil);

temp = 0;
for imeas = 1:nmeas
  temp = repmat(alpha(imeas,:).^2 ./ pi_k(imeas,:),np,1) .* repmat(denom_value(:,imeas),1,ncoil);
  temp(isnan(temp)) = 0;
  denom_value_new = denom_value_new + temp;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function value = c_bound(a,c,y_range)
value = interp1(y_range,c,a,'spline','extrap');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H = h_bloch(a, h)
r = -10*pi:pi/300:10*pi;
H = interp1(r,real(h),a,'spline','extrap') ...
	+ 1i * interp1(r,imag(h),a,'spline','extrap');

if(sum(isnan(H(:)))>0)
%	save error.mat a h
       	fail('Out of Range h_bloch %d or %d',max(a(:)),min(a(:)))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function H_prime = h_bloch_prime(a, h_prime)
r = -10*pi:pi/300:10*pi;
H_prime = interp1(r,real(h_prime),a,'spline','extrap') + i * interp1(r,imag(h_prime),a,'spline','extrap');

if (sum(isnan(H_prime(:)))>0)
	fail('Out of Range h_bloch_prime %d or %d',max(a(:)), min(a(:)))
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fprime] = F_2_prime(z, h_prime)
a = real(z);
b = imag(z);

z(z == 0) = eps;
Fprime = z.^2 ./ abs(z) .* h_bloch_prime(abs(z),h_prime);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% mri_b1map_ybar()
% forward model
% in
%	alpha [nmeasure ncoil]
%	zmap [np,ncoil]
%	omap [np]
% out
%	ybar [np,ncoil,ntip]
%
function ybar = mri_b1map_ybar(alpha, amap, bmap, omap, h);

[nmeasure ncoil] = size(alpha);
ybar = zeros(size(amap,1), nmeasure);
zmap = amap + 1i * bmap;

for it = 1:nmeasure
   ybar(:,it) = omap .* F_2(alpha(it,:) * zmap', h)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = F_2(z,h)
F = z .* h_bloch(abs(z),h);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% mri_b1map_kgrad_zmap3()
function [agrad, bgrad] = mri_b1map_kgrad_zmap3(yR, yI, alpha, a, b, omap, h, h_deriv);

[nmeasure ncoil] = size(alpha);

agrad = zeros(size(a));
bgrad = zeros(size(b));

z = a + 1i * b;

for ik=1:ncoil
for im=1:nmeasure
    g = omap .* F_2( alpha(im,:) * z', h)';
    %g = omap .* F_2( alpha(im,:) * a', alpha(im,:) * b', h)';
    gR = real(g);
    gI = imag(g);
    [da_F2R db_F2R da_F2I db_F2I] = F_2RI_prime2( alpha(im,:) * a', alpha(im,:) * b', h_deriv, h);

    temp_a1 = (yR(:,im) - gR) .* omap .* alpha(im,ik) .* da_F2R';
    temp_a2 = (yI(:,im) - gI) .* omap .* alpha(im,ik) .* da_F2I';
    agrad(:,ik) = agrad(:,ik) - temp_a1 - temp_a2;

    temp_b1 = (yR(:,im) - gR) .* omap .* alpha(im,ik) .* db_F2R';
    temp_b2 = (yI(:,im) - gI) .* omap .* alpha(im,ik) .* db_F2I';
    bgrad(:,ik) = bgrad(:,ik) - temp_b1 - temp_b2;
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% mri_b1map_omap()
%
function omap = mri_b1map_omap(yy, alpha, amap, bmap, h)
num = 0;
den = 0;

zmap = amap + 1i * bmap;

[nmeasure ncoil] = size(alpha);

for im = 1:nmeasure
    %num = num + real( conj(yy(:,im)) .* F_2( alpha(im,:) * amap', alpha(im,:) * bmap', h)');
    %den = den + abs(F_2(alpha(im,:) * amap', alpha(im,:) * bmap', h)').^2;
    num = num + real( conj(yy(:,im)) .* F_2(alpha(im,:) * zmap',h)');
    den = den + abs(F_2(alpha(im,:) * zmap',h)').^2;
end
den = den + eps;  % slight "regularization"
omap = num ./ den;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% mri_b1map_cost()
%
% Calculates the current cost of the function
%
		    function cost = mri_b1map_cost(yy, alpha, amap, bmap, omap, h, Rzmap_R, Rzmap_I)
cost = 0;
temp = 0;
temp2 = 0;
[nmeasure ncoil] = size(alpha);

zmap = amap + 1i * bmap;

for im=1:nmeasure
    temp = temp + .5 * abs( yy(:,im) - omap .* F_2( alpha(im,:) * zmap', h)').^2;
end
cost = sum(temp(:));

for ik=1:ncoil
	temp2 = temp2 + Rzmap_R(:,ik).penal(Rzmap_R(:,ik),amap(:,ik)) + ...
         Rzmap_I(:,ik).penal(Rzmap_I(:,ik),bmap(:,ik));
end
cost = cost + sum(temp2(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% plot_function_and_surrogate()
%
%
%
			    function plot_function_and_surrogate(yy,alpha,amap,bmap,omap,agrad,bgrad,rgrad_a,rgrad_b,Rzmap_R, Rzmap_I,h,denom_values,rdenom_a,rdenom_b,arg,pt,pt_coil)

% First plot the function
pixel_plot = [pt pt_coil];

current_cost = mri_b1map_cost(yy,alpha,amap,bmap,omap,h,Rzmap_R, Rzmap_I);

current_a = amap(pixel_plot(1),pixel_plot(2));
current_b = bmap(pixel_plot(1),pixel_plot(2));

if(current_a ~= 0)
    add_a = -.1*current_a:current_a*.01:.1*current_a;
    add_a_more = -current_a:current_a*.1:current_a;
else
    add_a = -.1:.01:.1;
    add_a_more = -1:.1:1;
end

for step_a = 1:size(add_a,2)
    amap(pixel_plot(1),pixel_plot(2)) = current_a + add_a(step_a);
    cost_function_a(step_a) = mri_b1map_cost(yy,alpha,amap,bmap,omap,h,Rzmap_R,Rzmap_I);
    amap(pixel_plot(1),pixel_plot(2)) = current_a + add_a_more(step_a);
    cost_function_a_more(step_a) = mri_b1map_cost(yy,alpha,amap,bmap,omap,h,Rzmap_R,Rzmap_I);
end

denom_value = denom_values(pixel_plot(1),pixel_plot(2)) + rdenom_a(pixel_plot(1),pixel_plot(2));
grad_value = agrad(pixel_plot(1),pixel_plot(2)) + rgrad_a(pixel_plot(1),pixel_plot(2));

           surrogate_a = current_cost* ones(size(add_a)) + add_a * grad_value + ...
               .5 * add_a.^2 * denom_value;

           surrogate_a_more = current_cost*ones(size(add_a_more)) + add_a_more * grad_value + ...
               .5 * add_a_more.^2 * denom_value;

grad_line_a = grad_value * add_a + current_cost;

% Reset amap
amap(pixel_plot(1),pixel_plot(2)) = current_a;

if(current_b ~= 0)
   add_b = -.1*current_b:current_b*.01:.1*current_b;
   add_b_more = -current_b:current_b*.1:current_b;
else
    add_b = -.1:.01:.1;
    add_b_more = -1:.1:1;
end

for step_b = 1:size(add_b,2)
    bmap(pixel_plot(1),pixel_plot(2)) = current_b + add_b(step_b);
    cost_function_b(step_b) = mri_b1map_cost(yy,alpha,amap,bmap,omap,h,Rzmap_R,Rzmap_I);
    bmap(pixel_plot(1),pixel_plot(2)) = current_b + add_b_more(step_b);
    cost_function_b_more(step_b) = mri_b1map_cost(yy,alpha,amap,bmap,omap,h,Rzmap_R,Rzmap_I);
end

denom_value = denom_values(pixel_plot(1),pixel_plot(2)) + rdenom_b(pixel_plot(1),pixel_plot(2));
grad_value = bgrad(pixel_plot(1),pixel_plot(2)) + rgrad_b(pixel_plot(1),pixel_plot(2));

surrogate_b = current_cost*ones(size(add_b)) + add_b * grad_value + ...
    .5 * add_b.^2 * denom_value;

surrogate_b_more = current_cost*ones(size(add_b_more)) + add_b_more * grad_value + ...
    .5 * add_b_more.^2 * denom_value;

grad_line_b = grad_value * add_b + current_cost;

%if(sum(surrogate_b<cost_function_b) > 0 || sum(surrogate_a<cost_function_a) > 0)
%if(sum(surrogate_a<cost_function_a) > 0)

%sum(surrogate_a<cost_function_a)
%sum(surrogate_b<cost_function_b)

if 1
prompt

%a_lim = [min(cost_function_a) max(cost_function_a)]
%b_lim = [min(cost_function_b) max(cost_function_b)]
%nrms(a_lim(1),a_lim(2))
%nrms(a_lim(1),a_lim(2)) > .01
%if (nrms(a_lim(1),a_lim(2)) > .01)
figure
subplot(2,2,1)
plot(add_a+current_a*ones(size(add_a)),cost_function_a,'r')
hold on
plot(add_a+current_a*ones(size(add_a)),surrogate_a,'g')
hold on
plot(add_a + current_a*ones(size(add_a)),grad_line_a,'c')
%ylim(a_lim)
title('Cost Function on the a axis (constant b)')
subplot(2,2,2)
plot(add_b+current_b*ones(size(add_b)),cost_function_b,'r')
hold on
plot(add_b+current_b*ones(size(add_b)),surrogate_b,'g')
hold on
plot(add_b+current_b*ones(size(add_b)),grad_line_b,'c')
%ylim(b_lim)
title('Cost Function on the b axis (constant a)')
subplot(2,2,3)
plot(add_a_more+current_a*ones(size(add_a_more)),cost_function_a_more,'r')
%hold on
%plot(add_a_more+current_a*ones(size(add_a_more)),surrogate_a_more)
subplot(2,2,4)
plot(add_b_more+current_b*ones(size(add_b_more)),cost_function_b_more,'r')
%hold on
%plot(add_b_more+current_b*ones(size(add_b_more)),surrogate_b_more)
%end
prompt
return

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% user function
% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default
iter = evalin('caller', 'iter');
niter = evalin('caller', 'arg.niter');
if iter == niter, disp 'done!', end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% showmaps
function showmaps(zlim_abs, zlim_phase, olim, ig)
iter = evalin('caller', 'iter');
% make line below say 10 for the standard, every ten times
if rem(iter,10), return, end
zmap = evalin('caller', 'zmap');
omap = evalin('caller', 'omap');
zmap = ig.embed(zmap);
omap = ig.embed(omap);
%blim = [0 0.1];
%plim = [-1 1] * 1.6;
%olim = [0 200];
im pl 3 1
im row 1
im(1, abs(zmap), 'B1map estimate (abs)', zlim_abs), cbar
im(2, angle(zmap), 'B1map estimate (phase)', zlim_phase), cbar
im(3, omap, 'Object estimate', olim), cbar
drawnow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% load_tables(profile)
% We have four tables - real & imaginary part of both first and second derivative
%
function [h h_deriv_rad h_deriv2_rad] = load_tables(profile_type)
if(strcmp(profile_type,'gauss'))
	load bloch_h_gauss_newnew_ub.mat;
elseif(strcmp(profile_type,'real_rect1'))
	load bloch_h_trunc_sinc_ub.mat;
elseif(strcmp(profile_type,'short_trunc_sinc'))
	load bloch_h_short_trunc_sinc.mat;
else
	r = -10*pi:pi/300:10*pi;
	h = sin(r)./r;
	middle = find(r == 0);
	h(middle) = 1;
	h_deriv_rad = (r .* cos(r) - sin(r)) ./ r.^2;
	h_deriv_rad(middle) = 0;
	%h_deriv2_rad = -sin(r) ./ r - 2 * cos(r) ./ r.^2 + 2*sin(r) ./ r.^3;
	h_deriv2_rad = (r.^2 .* -sin(r) - 2 * r .* cos(r) + 2 * sin(r)) ./ r.^3;
	h_deriv2_rad(middle) = -1/3;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%% load_denom_table(profile_type)
%
% Note - these are derived in testing_sin.m and testing2.m
function [c absy] = load_denom_table(profile_type)
if(strcmp(profile_type,'gauss'))
	load gauss_data.mat;
	c = norm_max;
elseif(strcmp(profile_type,'real_rect1'))
	load sinc_data.mat
	c = norm_max;
elseif(strcmp(profile_type,'short_trunc_sinc'))
	load short_trunc_sinc_data.mat
	c = norm_max;
else
	load sin_data.mat
	c = norm_max;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [da_F2R db_F2R da_F2I db_F2I] = F_2RI_prime2(a, b, h_prime, h)

absz = sqrt(a.^2 + b.^2);
absz(absz == 0) = eps;

h = h_bloch(absz,h);
hR = real(h);
hI = imag(h);

h_prime = h_bloch_prime(absz,h_prime);
d_hR = real(h_prime);
d_hI = imag(h_prime);

da_F2R = hR + a.^2 ./ absz .* d_hR - a.*b ./ absz .* d_hI;
db_F2R = a .* b ./ absz .* d_hR - hI - b.^2 ./ absz .* d_hI;

da_F2I = hI + a.^2 ./ absz .* d_hI + a .* b ./ absz .* d_hR;
db_F2I = a .* b./ absz .* d_hI + hR + b.^2 ./ absz .* d_hR;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mri_b1map_test_spatial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mri_b1map_test_spatial()
% tests spatial resolution of algorithm in two separate ways
%
function mri_b1map_test_spatial

%% Profile used to generate Data
arg.profile = 'real_rect1'
[h h_deriv h_deriv2] = load_tables(arg.profile);

f.dir = [path_find_dir('mri') '/../data/mri/'];
f.xtrue = [f.dir 'brainweb_t1.jpg'];
omap = double6(imread(f.xtrue)');
omap = omap(2:end-1,2:end-1); % make it 256^2
%omap = downsample2(omap, 16);
%omap = downsample2(omap,4);
omap = omap(:,2:end-1); % [64,62]

omap = ones(100);

nx = size(omap,1);
ny = size(omap,2);

ig = image_geom('nx', size(omap,1), 'ny', size(omap,2), 'fov', 25);

ncoil = 1;
alpha = [1; 2; 3];

[nmeasure ncoil] = size(alpha);

    [bmap pmap] = mri_sensemap_sim2('chat', 0, 'nx', ig.nx, 'ny', ig.ny, ...
              'dx', ig.dx, 'rcoil', 0.6* ig.fov);

bmap = bmap / max(sum(alpha,2));

zmap = bmap .* exp(1i .* pmap);

% Use just needed coils
zmap = zmap(:,:,1:ncoil);

zmap_true = zmap;
omap = ones(size(omap));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the matrices

A = kron(sparse(eye(nx*ny)),kron(sparse(alpha),eye(2)));
[nmeas ncoil] = size(alpha);

Q_j = kron(eye(nmeas),eye(2));
W = kron(sparse(eye(nx*ny)),Q_j);

i = zeros(4*nmeas*nx*ny,1);
i(1:2:end-1) = 1:2*nmeas*nx*ny;
i(2:2:end) = 1:2*nmeas*nx*ny;
j = zeros(size(i));
j(1:4:end) = 1:2:2*nmeas*nx*ny-1;
j(3:4:end) = 1:2:2*nmeas*nx*ny-1;
j(2:4:end) = 2:2:2*nmeas*nx*ny;
j(4:4:end) = 2:2:2*nmeas*nx*ny;

s = zeros(size(i));
zmap = reshape(zmap,[nx*ny ncoil]);
for imeas=1:size(alpha,1)
   [da_F2R db_F2R da_F2I db_F2I] = F_2RI_prime2(alpha(imeas,:) * real(zmap)', alpha(imeas,:) * imag(zmap)', h_deriv, h);
   s((imeas-1)*2 + 1:4*nmeas:end) = da_F2R .^2 + da_F2I .^2;
   s((imeas-1)*2 + 2:4*nmeas:end) = da_F2R .*  db_F2R + da_F2I .* db_F2I;
   s((imeas-1)*2 + 3:4*nmeas:end) = da_F2R .* db_F2R + da_F2I .* db_F2I;
   s((imeas-1)*2 + 4:4*nmeas:end) = db_F2R.^2 + db_F2I.^2;
end

W_new = sparse(i,j,s);

if 0  % old way to calculate W_new
   W_new = zeros(size(W));
   Q_j = zeros(2*nmeas,2*nmeas,nx*ny);
   zmap = reshape(zmap,[nx*ny ncoil]);
   for imeas=1:size(alpha,1)
     [da_F2R db_F2R da_F2I db_F2I] = F_2RI_prime2(alpha(imeas,:) * real(zmap)', alpha(imeas,:) * imag(zmap)', h_deriv, h);
     Q_j((imeas-1)*2+1,(imeas-1)*2+1,:) = da_F2R .^2 + da_F2I .^2;
     Q_j((imeas-1)*2+1,imeas*2,:) = da_F2R .*  db_F2R + da_F2I .* db_F2I;
     Q_j(imeas*2,(imeas-1)*2+1,:) = da_F2R .* db_F2R + da_F2I .* db_F2I;
     Q_j(imeas*2,imeas*2,:) = db_F2R.^2 + db_F2I.^2;
   end

   for ip = 1:nx*ny
     W_new((ip-1)*2*nmeas+1:ip*2*nmeas,(ip-1)*2*nmeas+1:ip*2*nmeas) = Q_j(:,:,ip);
   end
end

%% Okay - we want to try two sets of kappa and compare
% Set one - what I calculated in the algorithm

arg.kappa = 1;

[kappaR kappaI] = mri_b1map_kappa_calc(alpha, arg, nx, ny, ncoil, 0, reshape(zmap,[nx*ny ncoil]), omap(:));

% This kappa compensates for the alpha matrix and should match kappa2_new2
kappa2 = zeros(2,nx,ny);
kappa2(1,:,:) = kappaR;
kappa2(2,:,:) = kappaI;

%awa_new = A' * W_new * A;
%aa = A' * A;

%% kappa2_new uses the expected kappa, but will have varying FWHM based on alpha matrix used
%kappa2_new_temp = sqrt(awa_new ./ aa);
%kappa2_new = diag(kappa2_new_temp);

%% kappa2_new2 has a consistent FWHM regardless of alpha matrix and matches calculated one
%kappa2_new_temp2 = sqrt(awa_new);
%kappa2_new2 = diag(kappa2_new_temp2);

%kappa2_new = reshape(kappa2_new, [2 nx ny]);
%kappa2_new2 = reshape(kappa2_new2, [2 nx ny]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run through several betas

l2_beta_range = -6:.5:-1;
%l2_beta_range = [-3];

for each_beta = 1:length(l2_beta_range)

l2_beta = l2_beta_range(each_beta)
beta = 2^(l2_beta_range(each_beta))

  arg.order = 2;
checkReal = 1;

kappa1 = ones(2,nx,ny);

%% Just looking at 1 coil
Rzmap = Reg1(kappa1,'beta',beta,'order',arg.order, 'offsets',2*[1 nx nx+1 nx-1]);

our_mask = ones(2,nx,ny);

if checkReal
 arg.offset = [-1 0 0];
else
 arg.offset = [0 0 0];
end

%[psf var fwhm mtf con] = qpwls_psf(A,Rzmap,1,logical(our_mask),W,'offset',arg.offset,'fwhmtype','none');

A_ident = diag_sp(our_mask(:));

[psf var fwhm mtf con] = qpwls_psf(A_ident,Rzmap,1,logical(our_mask),1,'offset',arg.offset,'fwhmtype','none');

% Put an impulse into real and got out two peaks - real & imaginary ...

if checkReal
    fwhm_real_real_impulse(each_beta) = fwhm2(squeeze(psf(1,:,:)));
    %fwhm_imag_real_impulse(each_beta) = fwhm2(squeeze(psf(2,:,:)));
else
    %fwhm_real_imag_impulse(each_beta) = fwhm2(squeeze(psf(1,:,:)));
    fwhm_imag_imag_impulse(each_beta) = fwhm2(squeeze(psf(2,:,:)));
end

%% Okay let's calculate for our kappa corrected map .... and see how it looks

arg.kappa = 1;

%Rzmap2 = Reg1(kappa2_new,'beta', beta,'order',arg.order,'offsets',2*[1 nx nx+1 nx-1]);
Rzmap3 = Reg1(kappa2,'beta', beta,'order',arg.order,'offsets',2*[1 nx nx+1 nx-1]);

%[psf_compensated var fwhm mtf con] = qpwls_psf(A,Rzmap2,1,logical(our_mask),W_new,'offset',arg.offset,'fwhmtype','none');

[psf_uncompensated var fwhm mtf con] = qpwls_psf(A,Rzmap,1,logical(our_mask),W_new,'offset',arg.offset,'fwhmtype','none');

[psf_alpha_compensated var fwhm mtf con] = qpwls_psf(A,Rzmap3,1,logical(our_mask),W_new,'offset',arg.offset,'fwhmtype','none');

% Let's calculate our own psf ...

% Put an impulse into real and got out two peaks - real & imaginary ...

if checkReal
    %fwhm_real_real_impulse_compensated(each_beta) = fwhm2(squeeze(psf_compensated(1,:,:)));
    %fwhm_imag_real_impulse_compensated(each_beta) = fwhm2(squeeze(psf_compensated(2,:,:)));
    fwhm_real_real_impulse_uncompensated(each_beta) = fwhm2(squeeze(psf_uncompensated(1,:,:)));
    %fwhm_imag_real_impulse_uncompensated(each_beta) = fwhm2(squeeze(psf_uncompensated(2,:,:)));
    fwhm_real_real_impulse_alpha_compensated(each_beta) = fwhm2(squeeze(psf_alpha_compensated(1,:,:)));
else
    %fwhm_real_imag_impulse_compensated(each_beta) = fwhm2(squeeze(psf_compensated(1,:,:)));
    %fwhm_imag_imag_impulse_compensated(each_beta) = fwhm2(squeeze(psf_compensated(2,:,:)));
    %fwhm_real_imag_impulse_uncompensated(each_beta) = fwhm2(squeeze(psf_uncompensated(1,:,:)));
    fwhm_imag_imag_impulse_uncompensated(each_beta) = fwhm2(squeeze(psf_uncompensated(2,:,:)));
    fwhm_imag_imag_impulse_alpha_compensated(each_beta) = fwhm2(squeeze(psf_alpha_compnesated(2,:,:)));
end

if 1
%%%%%%%%% Impulse approach
printm('Now with impulses')

niter1pass = 1;
niter = 10;

zmap = reshape(zmap,[nx ny 1]);

%% Just do some testing with uniform zmap

%     zmap(60,60)
%     zmap = ones(size(zmap)) * zmap(nx/2,ny/2);

%ig = image_geom('nx', size(omap,1), 'ny', size(omap,2), 'fov', 25);
ig = image_geom('nx',size(omap,1), 'ny', size(omap,2),'dx', 1);

% Generate random data
ybar = mri_b1map_ybar(alpha, ig.maskit(real(zmap)), ig.maskit(imag(zmap)), omap(ig.mask), h);
ybar = ig.embed(ybar);
%rng(0)
sig = .00002;
yi = ybar + sig * randn(size(ybar)) + 1i * sig * randn(size(ybar));

%mri_b1map_omap(reshape(yi,[64*62 3]),alpha,reshape(zmap,[64*62 1]),h);

arg.kappa = 1;

%% without impulses

zmap_init = zmap + .01*randn(size(zmap));

[zh oh cost] = mri_b1map_rewrite3(yi, alpha, 'niter', niter, 'niter_1pass', niter1pass, ...
  'profile', arg.profile, 'l2b_zmap', l2_beta, 'l2b_zmap_1pass', -100, ...
  'kappa',arg.kappa,'scale',1);

[zhi ohi costi] = mri_b1map_rewrite3(yi, alpha, 'niter', 0, 'niter_1pass', 0, ...
   'profile', arg.profile, 'l2b_zmap', l2_beta, 'l2b_zmap_1pass', -100, ...
   'kappa',arg.kappa,'scale',1);

% save zhi.mat zhi

% Create impulses
[x y] = meshgrid(1:nx,1:ny);
%impulses = ( mod(x-7,10) == 0) & ( mod(y-7,10) == 0);
impulses = ( mod(x-10,20) == 0) & ( mod(y-10,20) == 0);
%impulses = (x == nx/2) & (y == ny/2);

%% Add in impulses to a temporary zmap (so we don't keep adding impulses over and over)
zmap_r = real(zmap);
zmap_i = imag(zmap);
zmap_temp(:,:,1) = (zmap_r(:,:,1) + .7 * max(zmap_r(:)) * impulses') + 1i * zmap_i(:,:,1);
%zmap_temp(:,:,1) = zmap_r(:,:,1) + .2 * impulses' + 1i*zmap_i(:,:,1);

ybar2 = mri_b1map_ybar(alpha, ig.maskit(real(zmap_temp)), ig.maskit(imag(zmap_temp)), omap(ig.mask), h);
ybar2 = ig.embed(ybar2);
yi2 = ybar2 + sig * randn(size(ybar2)) + 1i * sig * randn(size(ybar2));

zmap_temp_init = zmap_temp + .01*randn(size(zmap_temp));

[zh1 oh1 cost1] = mri_b1map_rewrite3(yi2, alpha, 'niter', niter, 'niter_1pass', niter1pass, ...
  'profile', arg.profile, 'l2b_zmap', l2_beta, 'l2b_zmap_1pass', -100, ...
   'kappa', arg.kappa,'scale',0);

o_diff = oh1 - oh;
zr_diff = real(zh1) - real(zh);
zi_diff = imag(zh1) - imag(zh);

num_impulses = sum(impulses(:))
num_impulses_x = max(sum(impulses,1))
num_impulses_y = max(sum(impulses,2))
[locX locY] = find(impulses == 1);

%keyboard

amt = 10;
num0 = 0;
for i=1:num_impulses
    if(zr_diff(locX(i),locY(i)) < 0)
        zr_diff = -zr_diff;
    end
    try
       fwhm_zreal(i) = fwhm2(zr_diff(locX(i)-amt:locX(i)+amt,locY(i)-amt:locY(i)+amt))
    catch
       num0 = num0 + 1
       fwhm_zreal(i) = 0;
    end
end

%keyboard

%save fwhm_zreal.mat fwhm_zreal each_beta

    fwhm_zreal_mean(each_beta) = mean(fwhm_zreal(fwhm_zreal(:)>0))
    fwhm_zreal_total(each_beta,:) = fwhm_zreal;

%save fwhm_zreal.mat fwhm_zreal_mean fwhm_zreal_total
end % end the impulse "if 0"

end % end the beta loop
%save fwhm_data.mat fwhm_zreal_total

figure
plot(l2_beta_range,fwhm_real_real_impulse,'o-')
hold on
%plot(beta_range,fwhm_real_real_impulse_compensated,':')
%hold on
plot(l2_beta_range,fwhm_real_real_impulse_uncompensated,'o--')
hold on
plot(l2_beta_range,fwhm_real_real_impulse_alpha_compensated,'o-.')
hold on
plot(l2_beta_range,fwhm_zreal_mean,'o:')
%legend('Ideal','No Kappa','Kappa adjusted')
legend('Ideal','No Kappas','Kappa adjusted','Kappa with impulses')
%legend('Ideal','Kappa Compensated','No Kappas','Kappa / alpha adjusted','Kappa with Impulses')
%legend('Ideal','Kappa Compensated','No Kappas','Kappa / alpha adjusted')
xlabel('log_2 \beta')
ylabel('FWHM for a real impulse')
title('\beta vs FWHM')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mri_b1map_test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mri_b1map_test()
% built-in test/example
%
function mri_b1map_test

niter1pass = 5;
niter = 50;
nreal = 1;

arg.filt = 1;

%load data_20_reals_500iters_30snr_k1measures

%for l2beta = [-8 -6 -4 -2 0]
for my_method = [2]
for my_pulse = [4]

savef = 1;
plotf = 1;
l2beta = -1
init = 1;
new_maps = 1;
updateo_length = 1;
compare_coil = 1;
old_inverse_method = 0;
argrotate = 1;

switch my_method
    case 1
        one_at_atime = 1
        leave_one_out = 0;
    case 2
        leave_one_out = 1
        one_at_atime = 0;
    otherwise
        printm('error')
        return
end

meas_2n = 1;
meas_n1 = 0;

scale = 1;

%% Profile used to generate Data

switch my_pulse
    case 1
        arg.profile = 'real_rect1'
    case 2
        arg.profile = 'sin'
    case 3
        arg.profile = 'gauss'
    case 4
        arg.profile = 'short_trunc_sinc'
    otherwise
        printm('Error - must pick proper method')
        return
end

f.dir = [path_find_dir('mri') '/../data/mri/'];
f.xtrue = [f.dir 'brainweb_t1.jpg'];
omap = double6(imread(f.xtrue)');
omap = omap(2:end-1,2:end-1); % make it 256^2
omap = downsample2(omap, 4);
omap = omap(:,2:end-1); % [64,62]

%% Mask and Bmap
olim = [0 200];
if(plotf)
    im plc 3 1
    im(3, omap, 'object', olim), cbar
end
mask = (omap>.1*max(omap(:)));

       masktemp = reshape(mask,size(omap));
       masktemp = imdilate(masktemp,strel('disk',5));
       mask_full = imerode(masktemp,strel('disk',6));


mask3 = omap < .3*max(omap(:));
mask3 = mask3 .* mask_full;


nx = size(omap,1);
ny = size(omap,2);

ig = image_geom('nx', size(omap,1), 'ny', size(omap,2), 'fov', 25);

if one_at_atime && meas_2n
    alpha = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 2 0 0 0; 0 2 0 0; 0 0 2 0; 0 0 0 2];
    alpha = alpha * 3
elseif one_at_atime && meas_n1
    alpha = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 2 0 0 0];
    alpha = alpha * 3;
elseif leave_one_out && meas_2n
    %alpha = [ones(4) - eye(4); 2*ones(4) - 2*eye(4)]
    alpha = [ones(4) - 2*eye(4); 2*ones(4) - 4*eye(4)]
    % Note - below would work, but code is not currently written to accomodate imaginary alpha
    %    alpha = alpha .* repmat([1 1i -1 -1i],[8 1]);
else
    alpha = [ones(4) - eye(4); 0 2 2 2];
    alpha = alpha .* repmat([1 1i -1 -1i],[5 1]);
end

[nmeasure ncoil] = size(alpha);

bmap = ir_mri_sensemap_sim('chat', 0, 'nx', ig.nx, 'ny', ig.ny, 'dx', ig.dx, ...
                        'rcoil', 0.6 * ig.fov);

big_full_mask = repmat(mask_full,[1 1 ncoil]);

blim = [0 max(abs(bmap(:).*big_full_mask(:)))];

if(plotf)
    im row 1
    im(1, abs(bmap), 'B1 maps', blim), cbar
end

      plim = [min(angle(bmap(:))) max(angle(bmap(:)))];

    bmap_temp = angle(bmap);
    bmap_temp(:,:,1) = unwrap(bmap_temp(:,:,1),[],2);

    for i=1:ncoil
        bmap_temp(:,:,i) = bmap_temp(:,:,i) + pi/2*(ncoil-i-1);
    end

if (plotf)
	im(2,bmap_temp,'Phase maps'), cbar
%	im(2, unwrap(angle(bmap),[],2), 'Phase maps'), cbar
	if (savef)
	%	print -deps -r600 fig_mri_b1map_sim1_true
	%	saveas(gcf,'fig_mri_b1map_sim1_true.fig')
	end
end

%keyboard

plim = [min(bmap_temp(:)) max(bmap_temp(:))];

zmap = bmap;
zlim_abs = blim;
zlim_phase = plim;

zmap = zmap(:,:,1:ncoil);

zmap_true = zmap;

%% Trick - add in the complex part of alpha here ...
for i=1:ncoil
   zmap(:,:,i) = zmap(:,:,i) * exp(1i * pi / 2 * (i-1));
end

[h h_deriv h_deriv2] = load_tables(arg.profile);

ybar = mri_b1map_ybar(alpha, ig.maskit(real(zmap)), ig.maskit(imag(zmap)), ...
   omap(ig.mask), h);
ybar = ig.embed(ybar);

% (complex) noise
rng(315);

sig = 3;

start_real = 1;
for ireal = start_real:nreal

yi = ybar + sig * randn(size(ybar)) + 1i * sig * randn(size(ybar));
snr = 10 * log10( sum(abs(ybar(:)).^2) / sum(abs(yi(:) - ybar(:)).^2) )
pr snr

if 1
%if(plotf)
	im clf
	im('pl', nmeasure/4+2, 4)
	clim = [0 max(abs(ybar(:)))];
	for it = 1:nmeasure
		im('notick',it,abs(yi(:,:,it)), clim)
		title(['Scan for y_',int2str(it)])
		if (it == 5)
			xlabelf('SNR = %4.1f dB', snr)
		end
	end
	if 0
	%	print -deps -r600 fig_mri_b1map_sim1_yi
	%	saveas(gcf,'fig_mri_b1map_sim1_yi.fig')
	end
end

keyboard

for j=1:ncoil
	masks(:,:,j) = mask;
	masks3(:,:,j) = mask3;
	masks_full(:,:,j) = mask_full;
end
numPixels = sum(masks(:));
numPixel = sum(mask(:));

tic
[zh1 oh1 cost1] = mri_b1map_rewrite3(yi, alpha, 'niter', niter, 'niter_1pass', niter1pass, ...
	'profile', arg.profile, 'l2b_zmap', l2beta, 'l2b_zmap_1pass', -10,'rotate',argrotate,...
	'isave',[niter/5 niter/2.5 niter]);
% 'userfun', @showmaps, 'userarg', {zlim_abs, zlim_phase, olim, ig});
toc


[zhi1 ohi1 costi1] = mri_b1map_rewrite3(yi, alpha, 'niter', 0, 'niter_1pass', 0, ...
					'profile', arg.profile, 'l2b_zmap', l2beta, 'l2b_zmap_1pass', -10, 'rotate',argrotate, ...
	'userfun', @showmaps, 'userarg', {zlim_abs, zlim_phase, olim, ig});

%keyboard

zh1_niter4 = zh1(:,:,:,1);
zh1_niter2 = zh1(:,:,:,2);
zh1_niter = zh1(:,:,:,3);

%% Trick - remove the complex part of alpha here ...
for i=1:ncoil
    zh1_niter(:,:,i) = zh1_niter(:,:,i) * exp(-1i * pi / 2 * (i-1));
    zh1_niter2(:,:,i) = zh1_niter2(:,:,i) * exp(-1i * pi / 2 * (i-1));
    zh1_niter4(:,:,i) = zh1_niter4(:,:,i) * exp(-1i * pi / 2 * (i-1));
    zhi1(:,:,i) = zhi1(:,:,i) * exp(-1i * pi / 2 * (i-1));
%    zmap(:,:,i) = zmap(:,:,i) * exp(-1i * pi / 2 * (i-1));
end


%% Created filtered initial estimate
switch my_method
    case 1
        width = 7;
    case 2
        width = 14;
end

h = fspecial('gaussian',[5 5], .0625*width);
for i=1:ncoil
   zhi1_filt(:,:,i) = imfilter(zhi1(:,:,i),h,'replicate');
end

if 0

%% I tested below using sin data
%% For OAAT, i = 7 was best
%% For INV, i = 14 was best

test = zeros(64,64);
test(32,32) = 1;

   for i=1:50
     h = fspecial('gaussian',[5 5], .0625*i);
     for j=1:4
        filtered_map(:,:,j,i) = imfilter(zhi1(:,:,j),h,'replicate');
     end
     psf(:,:,i) = imfilter(test,h,'replicate');
%     fwhm2(psf(:,:,i))
     nrms_err(i) = nrms(abs(filtered_map(:,:,:,i)).*masks,abs(zmap).*masks);
     roi_err(i) = nrms(abs(filtered_map(:,:,:,i)).*masks3,abs(zmap).*masks3);
   end

[y i] = min(nrms_err);
sigma = .0625*i;
i
nrms_err(i)
roi_err(i)
fwhm2(psf(:,:,i))
prompt

end


close all

%% Hand calculated DAM estimate
if meas_2n
     yy_alpha = yi(:,:,1:4);
     yy_alpha2 = yi(:,:,5:8);
     estimate_temp = acos(.5 * abs(yy_alpha2 ./ yy_alpha));

else
  yy_alpha = yi(:,:,1);
  yy_alpha2 = yi(:,:,5);
  estimate_temp(:,:,1) = acos(.5 * abs(yy_alpha2 ./ yy_alpha));
for i=2:4
  estimate_temp(:,:,i) = estimate_temp(:,:,1);
end
end

if meas_n1 || leave_one_out
	printm('here')
    for u=1:ncoil
        alpha_inv(u,:) = alpha(u,:);
    end

    estimate = alpha_inv^-1 * reshape(estimate_temp,[size(estimate_temp,1)*size(estimate_temp,2) size(estimate_temp,3)])';
    estimate = estimate';
    estimate = reshape(estimate,[size(estimate_temp,1) size(estimate_temp,2) size(estimate_temp,3)]);
else
    estimate = estimate_temp;
end

error_abs_bar = [-.07 .07];
error_phase_bar = [-pi/8 pi/8];


if compare_coil && (nreal==1) && 1

    zh1 = zh1_niter;

    %% Let's adjust the phase a bit for display purposes
    %% Put all on same "scale" for display
    %% Also unwrap if needed

    tempi = angle(zhi1);
    temp1 = angle(zh1);
    temp = angle(zmap_true);

    % First unwrap (this is messy because of noise - let's do my pseudo
    % unwrap)
    for i=1
        %tempi(:,:,i) = unwrap(angle(zhi1(:,:,i)),[],2);
        %temp1(:,:,i) = unwrap(angle(zh1(:,:,i)),[],2);
        temp(:,:,i) = unwrap(angle(zmap_true(:,:,i)),[],2);
        temp_i = tempi(:,:,i);
        temp_1 = temp1(:,:,i);
        temp_i(temp_i > - pi/2) = temp_i(temp_i > -pi/2) - 2*pi;
        temp_1(temp_1 > - pi/2) = temp_1(temp_1 > -pi/2) - 2*pi;
        tempi(:,:,i) = temp_i;
        temp1(:,:,i) = temp_1;
    end

    for i=1:ncoil
        tempi(:,:,i) = tempi(:,:,i) + pi/2*(ncoil-i-1);
        temp1(:,:,i) = temp1(:,:,i) + pi/2*(ncoil-i-1);
        temp(:,:,i) = temp(:,:,i) + pi/2*(ncoil-i-1);
    end



% save data_for_picture.mat masks zlim_abs error_abs_bar error_phase_bar zlim_phase estimate zh1 zhi1 zmap

     figure
     im pl 5 5
     im row 2
   im('notick', 1, abs(zhi1), 'Initial est.', zlim_abs)
   im('notick',2, abs(zh1), 'Reg. est. (abs)', zlim_abs)

  im('notick',6,(abs(zhi1)-abs(zmap_true)).*masks, 'Masked init. error', error_abs_bar)
  im('notick',7,(abs(zh1)-abs(zmap_true)).*masks, 'Masked reg. error', error_abs_bar)

   im('notick',11,angle(zhi1),'Initial (phase) est.',zlim_phase)
   im('notick',12,angle(zh1),'Reg. est. (phase)', zlim_phase)

  im('notick',16,(angle(zhi1) - angle(zmap_true)).*masks,'Masked init. error', error_phase_bar)
     im('notick',17,(angle(zh1) - angle(zmap_true)).*masks,'Masked reg. error', error_phase_bar)

     figure
     im row 2
     nrow = 4;
     ncol = 2;
     hmargin = .05;
     wmargin = 0;
     height = (1-hmargin*(nrow+1))/nrow;
     %width = (1-wmargin*(ncol+1))/ncol;
     width = height;

     position = [1 1];
     subplot('position',[wmargin*position(2)+width*(position(2)-1) hmargin*(nrow-position(1)+1)+height*(nrow-position(1)) width height])
     im('notick',abs(zhi1),'Initial abs est.',zlim_abs)

     position = [1 2];
     subplot('position',[wmargin*position(2)+width*(position(2)-1) hmargin*(nrow-position(1)+1)+height*(nrow-position(1)) width height])
     im('notick',abs(zh1),'Reg. abs est.',zlim_abs)

     position = [2 1];
     subplot('position',[wmargin*position(2)+width*(position(2)-1) hmargin*(nrow-position(1)+1)+height*(nrow-position(1)) width height])
     im('notick',(abs(zhi1)-abs(zmap_true)).*masks,'Masked init. abs error',error_abs_bar)

     position = [2 2];
     subplot('position',[wmargin*position(2)+width*(position(2)-1) hmargin*(nrow-position(1)+1)+height*(nrow-position(1)) width height])
     im('notick',(abs(zh1)-abs(zmap_true)).*masks,'Masked reg. abs error',error_abs_bar)

     position = [3 1];
     subplot('position',[wmargin*position(2)+width*(position(2)-1) hmargin*(nrow-position(1)+1)+height*(nrow-position(1)) width height])
     im('notick',tempi,'Init. phase est.',zlim_phase)

     position = [3 2];
     subplot('position',[wmargin*position(2)+width*(position(2)-1) hmargin*(nrow-position(1)+1)+height*(nrow-position(1)) width height])
     im('notick',temp1,'Reg. phase est.',zlim_phase)

     position = [4 1];
     subplot('position',[wmargin*position(2)+width*(position(2)-1) hmargin*(nrow-position(1)+1)+height*(nrow-position(1)) width height])
     im('notick',(tempi - temp).*masks,'Masked init. phase error.',error_phase_bar)

     position = [4 2];
     subplot('position',[wmargin*position(2)+width*(position(2)-1) hmargin*(nrow-position(1)+1)+height*(nrow-position(1)) width height])
     im('notick',(temp1 - temp).*masks,'Masked reg. phase error',error_phase_bar)


if 0
	if one_at_atime && meas_2n
	%	ir_savefig('fig_oaat_2nx3_dam')
	%	saveas(gcf,'fig_oaat_2nx3_dam')
	%	ir_savefig('fig_oaat_2n_beta-1_150iters')
	%	saveas(gcf,'fig_oaat_2n_beta-1_150iters.fig')
	%	save oaat_2n_beta-1.mat estimate zh1 zhi1 niter alpha ...
	%	zlim_abs zlim_phase error_abs_bar error_phase_bar zmap masks
	elseif one_at_atime && meas_n1
	%	ir_savefig('fig_oaat_n1_beta-4x3_1000iters')
	%	saveas(gcf,'fig_oaat_n1_beta-4x3_1000iters.fig')
	%	save oaat_n1_beta-4x3_1000iters.mat estimate zh1 zhi1 ...
	%	niter alpha zlim_abs zlim_phase error_abs_bar ..
	%	error_phase_bar zmap masks
	elseif leave_one_out && meas_2n
	%	ir_savefig('fig_loo_2n_beta-1_150iters')
	%	saveas(gcf,'fig_loo_2n_beta-1_150iters.fig')
	%	save loo_2n.mat estimate zh1 zhi1 niter alpha zlim_abs ...
	%	zlim_phase error_abs_bar error_phase_bar zmap masks
	else
	%	ir_savefig('fig_loo_n1_beta-1_500iters')
	%	saveas(gcf,'fig_loo_n1_beta-1_500iters.fig')
	%	save loo_n1_beta-1_500iters.mat estimate zh1 zhi1 niter alpha zlim_abs zlim_phase error_abs_bar error_phase_bar zmap masks
	end
end

% save image_data.mat estimate zh1 zhi1 niter alpha zlim_abs zlim_phase error_abs_bar error_phase_bar zmap masks omap oh1 ohi1 masks3
end


printm('Error initial estimate abs')
	error_abs_init(my_method,my_pulse,ireal) = nrms(abs(zhi1).*masks,abs(zmap_true).*masks)
printm('Error filtered init abs')
	error_abs_init_filt(my_method,my_pulse,ireal) = nrms(abs(zhi1_filt).*masks,abs(zmap_true).*masks)
	printm('Error coil estimate abs')
	error_abs(my_method,my_pulse,ireal,1) = nrms(abs(zh1_niter4).*masks,abs(zmap_true).*masks);
	error_abs(my_method,my_pulse,ireal,2) = nrms(abs(zh1_niter2).*masks,abs(zmap_true).*masks);
	error_abs(my_method,my_pulse,ireal,3) = nrms(abs(zh1_niter).*masks,abs(zmap_true).*masks)
	printm('Error initial estimate phase')
	error_angle_init(my_method,my_pulse,ireal) = nrms(exp(1i * angle(zhi1)).*masks,exp(1i * angle(zmap_true)).*masks)
printm('Error filtered init phase')
  error_angle_init_filt(my_method,my_pulse,ireal) = nrms(exp(1i * angle(zhi1_filt)).*masks,exp(1i * angle(zmap_true)).*masks)
  error_angle(my_method,my_pulse,ireal,1) = nrms(exp(1i * angle(zh1_niter4)).*masks,exp(1i * angle(zmap_true)).*masks);
  error_angle(my_method,my_pulse,ireal,2) = nrms(exp(1i * angle(zh1_niter2)).*masks,exp(1i * angle(zmap_true)).*masks);
  printm('Error coil estimate phase')
  error_angle(my_method,my_pulse,ireal,3) = nrms(exp(1i * angle(zh1_niter)).*masks,exp(1i * angle(zmap_true)).*masks)

% Low mag error
printm('Error initial estimate abs - low mag')
  error_abs_init_low(my_method,my_pulse,ireal) = nrms(abs(zhi1).*masks3,abs(zmap_true).*masks3)
printm('Error filtered init abs - low mag')
  error_abs_init_filt_low(my_method,my_pulse,ireal) = nrms(abs(zhi1_filt).*masks3,abs(zmap_true).*masks3)
  error_abs_low(my_method,my_pulse,ireal,1) = nrms(abs(zh1_niter4).*masks3,abs(zmap_true).*masks3);
  error_abs_low(my_method,my_pulse,ireal,2) = nrms(abs(zh1_niter2).*masks3,abs(zmap_true).*masks3);
  printm('Error coil estimate abs - low mag')
  error_abs_low(my_method,my_pulse,ireal,3) = nrms(abs(zh1_niter).*masks3,abs(zmap_true).*masks3)

printm('Error initial estimate phase - low mag')
  error_angle_init_low(my_method,my_pulse,ireal) = nrms(exp(1i * angle(zhi1)).*masks3,exp(1i * angle(zmap_true)).*masks3)
printm('Error filtered init phase - low mag')
   error_angle_init_filt_low(my_method,my_pulse,ireal) = nrms(exp(1i * angle(zhi1_filt)).*masks3,exp(1i * angle(zmap_true)).*masks3)
  error_angle_low(my_method,my_pulse,ireal,1) = nrms(exp(1i * angle(zh1_niter4)).*masks3,exp(1i * angle(zmap_true)).*masks3);
  error_angle_low(my_method,my_pulse,ireal,2) = nrms(exp(1i * angle(zh1_niter2)).*masks3,exp(1i * angle(zmap_true)).*masks3);
  printm('Error coil estimate phase - low mag')
  error_angle_low(my_method,my_pulse,ireal,3) = nrms(exp(1i * angle(zh1_niter)).*masks3,exp(1i * angle(zmap_true)).*masks3)


sprintf('ireal: %d',ireal)
sprintf('my_pulse: %s',arg.profile)
sprintf('my_method: %d', my_method)

if(nreal > 1)
%	save data_20_reals_beta-1_051509.mat error_abs_init error_abs error_abs_init_low error_abs_low error_angle error_angle_init error_angle_init_low error_angle_low error_abs_init_filt error_angle_init_filt error_abs_init_filt_low error_angle_init_filt_low
end

end % end of ireal loop
end % end of my_pulse loop
end % end of my_method loop
%end % end of l2_beta loop

%save data_many_beta_30snr_k1measures_loo.mat error_abs_init error_abs error_abs_init_low error_abs_low error_angle error_angle_init error_angle_init_low error_angle_low zhi1_all zh1_all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mri_b1map_test_spatial2()
% tests spatial resolution of algorithm in two separate ways
%
function mri_b1map_test_spatial2

%% Profile used to generate Data
arg.profile = 'real_rect1'
[h h_deriv h_deriv2] = load_tables(arg.profile);

f.dir = [path_find_dir('mri') '/../data/mri/'];
f.xtrue = [f.dir 'brainweb_t1.jpg'];
omap = double6(imread(f.xtrue)');
omap = omap(2:end-1,2:end-1); % make it 256^2
%omap = downsample2(omap, 16);
%omap = downsample2(omap,4);
omap = omap(:,2:end-1); % [64,62]

omap = ones(100);

nx = size(omap,1);
ny = size(omap,2);

ig = image_geom('nx', size(omap,1), 'ny', size(omap,2), 'fov', 25);

ncoil = 4;
%ncoil = 1;
%alpha = [1; 2; 3];
alpha = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 2 0 0 0];
%alpha = [1 1 1 0; 1 1 0 1; 1 0 1 1; 0 1 1 1; 2 2 2 0];
%alpha = [1 1 1 0; 1 1 0 1; 1 0 1 1; 0 1 1 1; 2 2 2 0; 2 2 0 2; 2 0 2 2; 0 2 2 2];
%alpha = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 2 0 0 0; 0 2 0 0; 0 0 2 0; 0 0 0 2];

[nmeasure ncoil] = size(alpha);

% see dat/arch/student/funai,amanda/b1map_slice/bloch
[bmap pmap] = mri_sensemap_sim2('chat', 0, 'nx', ig.nx, 'ny', ig.ny, ...
	'dx', ig.dx, 'rcoil', 0.6* ig.fov);

bmap = bmap / max(sum(alpha,2));

zmap = bmap .* exp(1i .* pmap);

% Use just needed coils
zmap = zmap(:,:,1:ncoil);

zmap_true = zmap;
omap = ones(size(omap));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run through several betas

%l2_beta_range = -6:.5:-1;
l2_beta_range = [-1];

for each_beta = 1:length(l2_beta_range)

l2_beta = l2_beta_range(each_beta)
beta = 2^(l2_beta_range(each_beta))

arg.order = 2;

printm('Impulses')

niter1pass = 5;
niter = 100;

%zmap = reshape(zmap,[nx ny 4]);
size(zmap)
size(alpha)

% Generate random data
ybar = mri_b1map_ybar(alpha, ig.maskit(real(zmap)), ig.maskit(imag(zmap)), omap(ig.mask), h);
ybar = ig.embed(ybar);
rng(0)
sig = .00002;
yi = ybar + sig * randn(size(ybar)) + 1i * sig * randn(size(ybar));

%mri_b1map_omap(reshape(yi,[64*62 3]),alpha,reshape(zmap,[64*62 1]),h);

arg.kappa = 1;

%% without impulses

zmap_init = zmap + .01*randn(size(zmap));

[zh oh cost] = mri_b1map_rewrite3(yi, alpha, 'niter', niter, 'niter_1pass', niter1pass, ...
  'profile', arg.profile, 'l2b_zmap', l2_beta, 'l2b_zmap_1pass', -100, ...
   'kappa',arg.kappa,'scale',1);
%  'kappa',arg.kappa,'scale',1,'init_omap',omap(:));

[zhi ohi costi] = mri_b1map_rewrite3(yi, alpha, 'niter', 0, 'niter_1pass', 0, ...
   'profile', arg.profile, 'l2b_zmap', l2_beta, 'l2b_zmap_1pass', -100, ...
   'kappa',arg.kappa,'scale',1);
%   'kappa',arg.kappa,'scale',1,'init_omap',omap(:));

% Create impulses
[x y] = meshgrid(1:nx,1:ny);
%impulses = ( mod(x-7,10) == 0) & ( mod(y-7,10) == 0);
impulses = ( mod(x-10,20) == 0) & ( mod(y-10,20) == 0);
%impulses = (x == nx/2) & (y == ny/2);

%% Add in impulses to a temporary zmap (so we don't keep adding impulses over and over)
zmap_r = real(zmap);
zmap_i = imag(zmap);
zmap_temp(:,:,1) = (zmap_r(:,:,1) + .2 * max(zmap_r(:)) * impulses') + 1i * zmap_i(:,:,1);
%zmap_temp(:,:,1) = zmap_r(:,:,1) + .2 * impulses' + 1i*zmap_i(:,:,1);
for j=2:ncoil
    zmap_temp(:,:,j) = zmap_r(:,:,j) + 1i*zmap_i(:,:,j);
end

size(zmap_temp)
size(alpha)

ybar2 = mri_b1map_ybar(alpha, ig.maskit(real(zmap_temp)), ig.maskit(imag(zmap_temp)), omap(ig.mask), h);
ybar2 = ig.embed(ybar2);
yi2 = ybar2 + sig * randn(size(ybar2)) + 1i * sig * randn(size(ybar2));

zmap_temp_init = zmap_temp + .01*randn(size(zmap_temp));

[zh1 oh1 cost1] = mri_b1map_rewrite3(yi2, alpha, 'niter', niter, 'niter_1pass', niter1pass, ...
  'profile', arg.profile, 'l2b_zmap', l2_beta, 'l2b_zmap_1pass', -100, ...
'kappa',arg.kappa,'scale',1);
%   'kappa', arg.kappa,'scale',1,'init_omap',omap(:));

nrms(abs(zhi),abs(zmap))
nrms(abs(zh1),abs(zmap))
nrms(abs(zh),abs(zmap))

o_diff = oh1 - oh;
zr_diff = real(zh1) - real(zh);
zi_diff = imag(zh1) - imag(zh);

num_impulses = sum(impulses(:))
num_impulses_x = max(sum(impulses,1))
num_impulses_y = max(sum(impulses,2))
[locX locY] = find(impulses == 1);

%keyboard

amt = 10;
num0 = 0;
for i=1:num_impulses
    if(zr_diff(locX(i),locY(i)) < 0)
        zr_diff = -zr_diff;
    end
    try
       fwhm_zreal(i) = fwhm2(zr_diff(locX(i)-amt:locX(i)+amt,locY(i)-amt:locY(i)+amt))
    catch
       num0 = num0 + 1
       fwhm_zreal(i) = 0;
    end
end

%keyboard
%save fwhm_zreal.mat fwhm_zreal each_beta

    fwhm_zreal_mean(each_beta) = mean(fwhm_zreal(fwhm_zreal(:)>0))
    fwhm_zreal_total(each_beta,:) = fwhm_zreal;

%save fwhm_zreal.mat fwhm_zreal_mean fwhm_zreal_total

end % end the beta loop
%save fwhm_data.mat fwhm_zreal_total

plot(l2_beta_range,fwhm_zreal_mean)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% mri_b1map_test_less
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% mri_b1map_test_less()
% built-in test/example
%
function mri_b1map_test_less

niter1pass = 5;
niter = 250;
nreal = 1;

l2beta = -3
init = 1;
new_maps = 1;
updateo_length = 1;
compare_coil = 1;
old_inverse_method = 0;
argrotate = 1;

leave_one_out = 1
one_at_atime = 0;

meas_2n = 1;
meas_n1 = 0;

scale = 1;

%% Profile used to generate Data


arg.profile = 'real_rect1'

f.dir = [path_find_dir('mri') '/../data/mri/'];
f.xtrue = [f.dir 'brainweb_t1.jpg'];
omap = double6(imread(f.xtrue)');
omap = omap(2:end-1,2:end-1); % make it 256^2
omap = downsample2(omap, 4);
omap = omap(:,2:end-1); % [64,62]

%% Mask and Bmap
olim = [0 200];
mask = (omap>.1*max(omap(:)));

       masktemp = reshape(mask,size(omap));
       masktemp = imdilate(masktemp,strel('disk',5));
       mask_full = imerode(masktemp,strel('disk',6));


mask3 = omap < .3*max(omap(:));
mask3 = mask3 .* mask_full;


nx = size(omap,1);
ny = size(omap,2);

ig = image_geom('nx', size(omap,1), 'ny', size(omap,2), 'fov', 25);

if one_at_atime && meas_2n
    alpha = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 2 0 0 0; 0 2 0 0; 0 0 2 0; 0 0 0 2];
   alpha = alpha * 3;
elseif one_at_atime && meas_n1
    alpha = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 2 0 0 0];
    alpha = alpha * 3;
elseif leave_one_out && meas_2n
    alpha = [ones(4) - eye(4); 2*ones(4) - 2*eye(4)];
    %% Let's just use 4 scans (then we'll see about 2 ... have my doubts)
    alpha = alpha([1:6],:);
else
    alpha = [ones(4) - eye(4); 0 2 2 2];
end

[nmeasure ncoil] = size(alpha);

if(new_maps)
    [bmap pmap] = mri_sensemap_sim2('chat', 0, 'nx', ig.nx, 'ny', ig.ny, ...
              'dx', ig.dx, 'rcoil', 0.6* ig.fov);
else
   bmap = mri_sensemap_sim('chat', 0, 'nx', ig.nx, 'ny', ig.ny, 'dx', ig.dx, ...
		'rcoil', 0.6 * ig.fov); % make amplitude drop off a lot
end

blim = [0 max(bmap(:))];
blim = [0 .84];

plim = [1.0 pi/2];
plim = [min(pmap(:)) max(pmap(:))];

zmap = bmap .* exp(1i .* pmap);
zlim_abs = [0 .84];
zlim_phase = [1.0 pi/2];

zmap = zmap(:,:,1:ncoil);

[h h_deriv h_deriv2] = load_tables(arg.profile);


ybar = mri_b1map_ybar(alpha, ig.maskit(real(zmap)), ig.maskit(imag(zmap)), ...
   omap(ig.mask), h);
ybar = ig.embed(ybar);

% (complex) noise
rng(315)

%sig = 4.04;  %20dB
%sig = 2.275  %25dB
sig = 1.2803;  %30 dB
ireal = 1;

yi = ybar + sig * randn(size(ybar)) + 1i * sig * randn(size(ybar));
snr = 10 * log10( sum(abs(ybar(:)).^2) / sum(abs(yi(:) - ybar(:)).^2) )
pr snr

for j=1:ncoil
    masks(:,:,j) = mask;
    masks3(:,:,j) = mask3;
    masks_full(:,:,j) = mask_full;
end
numPixels = sum(masks(:));
numPixel = sum(mask(:));

%yy_alpha = yi(:,:,1:2);
%yy_alpha2 = yi(:,:,3:4);
%estimate_temp = acos(.5 * abs(yy_alpha2 ./ yy_alpha));

%keyboard

tic
[zh1 oh1 cost1] = mri_b1map_rewrite3(yi, alpha, 'niter', niter, 'niter_1pass', niter1pass, ...
				     'profile', arg.profile, 'l2b_zmap', l2beta, 'l2b_zmap_1pass', -10,'rotate',argrotate);
toc


[zhi1 ohi1 costi1] = mri_b1map_rewrite3(yi, alpha, 'niter', 0, 'niter_1pass', 0, ...
					'profile', arg.profile, 'l2b_zmap', l2beta, 'l2b_zmap_1pass', -10, 'rotate',argrotate, ...
   'userfun', @showmaps, 'userarg', {zlim_abs, zlim_phase, olim, ig});


close all

error_abs_bar = [-.07 .07];
error_phase_bar = [-pi/8 pi/8];

figure
im pl 5 5
im row 2
im('notick', 1, abs(zhi1), 'Initial est.', zlim_abs)
im('notick',2, abs(zh1), 'Reg. est. (abs)', zlim_abs)

im('notick',6,(abs(zhi1)-abs(zmap)).*masks, 'Masked init. error', error_abs_bar)
im('notick',7,(abs(zh1)-abs(zmap)).*masks, 'Masked reg. error', error_abs_bar)

im('notick',11,angle(zhi1),'Initial (phase) est.',zlim_phase)
im('notick',12,angle(zh1),'Reg. est. (phase)', zlim_phase)

im('notick',16,(angle(zhi1) - angle(zmap)).*masks,'Masked init. error', error_phase_bar)
im('notick',17,(angle(zh1) - angle(zmap)).*masks,'Masked reg. error', error_phase_bar)



printm('Error initial estimate abs')
nrms(abs(zhi1).*masks,abs(zmap).*masks)
printm('Error abs')
nrms(abs(zh1).*masks,abs(zmap).*masks)
printm('Error initial estimate phase')
nrms(exp(1i * angle(zhi1)).*masks,exp(1i * angle(zmap)).*masks)
printm('Error phase')
nrms(exp(1i * angle(zh1)).*masks,exp(1i * angle(zmap)).*masks)


printm('Error initial estimate abs - LOW')
nrms(abs(zhi1).*masks3,abs(zmap).*masks3)
printm('Error abs - LOW')
nrms(abs(zh1).*masks3,abs(zmap).*masks3)
printm('Error initial estimate phase - LOW')
nrms(exp(1i * angle(zhi1)).*masks3,exp(1i * angle(zmap)).*masks3)
printm('Error phase - LOW')
nrms(exp(1i * angle(zh1)).*masks3,exp(1i * angle(zmap)).*masks3)

%save temp.mat zhi1 zh1 zmap masks masks3 omap oh1 ohi1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mri_b1map_function_test

S = strvcat('k-', 'k:', 'k--');

im clf
for my_pulse = 1:3
	switch my_pulse
	case 1
		arg.profile = 'real_rect1';
	case 2
		arg.profile = 'sin';
	case 3
		arg.profile = 'gauss';
	otherwise
		fail('Error - must pick proper method')
	end

	[h h_deriv h_deriv2] = load_tables(arg.profile);

	alpha = 0:pi/300:pi/2;

	% assume noiseless model and fj=1

	y1 = F_2(alpha,h);
	y2 = F_2(alpha*2,h);

	alpha_hat(:,my_pulse) = acos(.5 * abs(y2 ./ y1));

	if im
		plot(alpha, alpha_hat(:,my_pulse), S(my_pulse,:))
		hold on
	end
end
hold off

xlabel('True Alpha')
ylabel('Noiseless DAM alpha estimate')
legend('Truncated sinc','Ideal sinc','Gaussian','Location','NorthWest')

%save function_dam_data.mat alpha alpha_hat
