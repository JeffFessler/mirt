 function xsense = desense(xalias, smap, varargin)
%function xsense = desense(xalias, smap, varargin)
% compute 2D images using traditional cartesian SENSE recon (Pruessman 99
% paper)
% Example usage:
%     xsense = desense(xalias, smap, 'phi', phi, 'sos', sosmap);
% 
% This is the Oct 2010 update to mri_cartesian_sense_recon_new2.m
% In this version:
%   auto finding of px with high sensitivity for regularization
%   less looping (hopefully) for speedier execution
%   xalias can be ny/2 x nz x nc  -OR- ny x nz x nc,
%   - alg can deal with either correctly
%
% in:
%   xalias 	[(Nr) nc]   	aliased coil images (2D or 3D)
%   smap 	[(N) nc]        complex sensemaps
%
% options
%   'phi' 	[nc nc]         receiver noise matrix
%   'sos'   [(N)]           sum of squares of smaps 
%                           including it as an input saves time
%                           computed within this fn if not provided
% out:
%   xsense   	[(N)]      	reconstructed image
%
% K. Khalsa, Oct 2010


% establish defaults and pair optional inputs
arg.phi = 1;    
arg.sos = [];
arg = vararg_pair(arg, varargin);

if (size(arg.phi,1) ~= size(arg.phi,2)), 
    error 'receiver noise matrix must be square! nc x nc'
else
    nc1 = size(arg.phi,1);
end

% extract sizes of aliased and full images
[nya, nza, nc2] = size(xalias);
 
sdim = size(smap);
nc3 = sdim(end);
sdim = sdim(1:end-1);  	% [nx (ny) (nz)]
ny = sdim(1);       % 1st dim for smap2D (slice in x) ONLY
nz = sdim(2);       % 2nd dim for smap2D (slice in x) ONLY

% figure out offset pairs based on size of xalias
if (nya == ny)
    y_offset = [0 ny/2];
else
    y_offset = ny/4 + [0 ny/2];
end

% some error checking
if (~isequal(nc2, nc3) || ((~isequal(nc1, 1)) && (~isequal(nc1, nc2))))
    error(sprintf('number of coils must be the same in smap, y and phi!'))
else
    nc = nc2;
end

if isempty(arg.sos)
    arg.sos = sqrt(sum(abs(smap).^2, 3));
end

invphi = arg.phi \ eye(size(arg.phi));



% define lambda, regularization parameter for SENSE recon
% lam = 0.01 * max(diag(SPS)), for a pixel in high signal area
% equiv: lam = 0.01 * max(sos);
% from bydder:07:oos

lam = 0.01 * max(arg.sos(:));
% *************TMP for testing!!!!**************
lam = 0;

xsense = zeros(ny, nz);
pxLoopStart = tic;

for iz = 1:nz

    for iy = 1:ny/2
%         iyatmp = iy + ny/4 + [0 ny/2];    % maybe this?
        iyatmp = iy + y_offset;
        iya = mod(iyatmp - 1, ny) + 1;
        
        a = squeeze(xalias(iy, iz, :));
        S = squeeze(smap(iya, iz, :)).';
        SP = S' * invphi;
        SPS = S' * invphi * S + lam * eye(2);
%         U = SPS \ SP;
%         v = U * a;
        v = SPS \ (SP * a);
        xsense(iya, iz) = v;        
        
	if (iy == 1 && iz == nz/2+1)
	         disp('desense: check iya pairs and smap, xalias');
	         keyboard
        end

    end
end
pxLoopTime = toc(pxLoopStart)


end     % end function





function old_desense_parts()
% \begin{oldWay}, was very very slow....
if 0
for ip = 1:Nr   
    ticker(mfilename, ip, Nr);
%     zoffset = floor( (ip - 1)/nyr);
%     yoffset = mod( (ip-1), nyr) + 1;

    % index pixel locations on aliased images
    [yoffset, zoffset] = ind2sub([nyr nz], ip);  % same as above, but cleaner
%     zoffset = zoffset - 1;  % old: not needed with sub2ind
    
    % translate into locations of full-sized image
    %  assume same center of FOV as aliased, hence the ny/4 factor
%     ips = zoffset*ny + yoffset + floor(ny/4);     % old, less clean
    ips = sub2ind([ny nz], yoffset + floor(ny/4), zoffset);       % new: cleaner
    ips1 = ips;
    if ( (yoffset + floor(ny/4) + ny/2) <= ny) 
        ips2 = ips + ny/2;
    else
        ips2 = ips - ny/2;
    end
    ipss(ip, :) = [ips1 ips2];  % smap coords for overlapping pixels
%     sprintf('ip = %d, zoffset = %d, yoffset = %d, ips1 = %d, ips2 = %d',ip,zoffset, yoffset, ips1, ips2)

    a = xc(ip,:).';
    S = sc([ips1 ips2],:).';

    % 9/30/10: replaced longer calculation of v by v = S\a;
%     SPS = S' * phi * S + lambda * eye(2);
%     SPSinv = 1 / (det(SPS)) * [SPS(2,2), -SPS(1,2); -SPS(2,1), SPS(1,1)];
%     U = SPSinv * S' * phiinv;
%     v = U * a;
   % keyboard
  v = S \ a;	
    xsense([ips1 ips2]) = v;
end

% keyboard
% reshape xsense
xsense = reshape(xsense, sdim);  
% \end{oldWay}
end

end     % end old_parts function



    
    
