function [w0, t] = adw_fan(ob, wi, mask, Phi, varargin)
% function [w0, phi, t] = adw_fan(ob, wi)
% Compute the angular-dependent weighting for 3rd gen fan-beam
% w0(\Phi) = w(x0, y0, \Phi) 
%         for fully corrected penalty:  
%          = w(s',\beta')*J(s')|\phi'=\Phi + ... 
%            w(s',\beta')*J(s')|\phi'=\Phi-pi 
% see fessler chapter 3.5.3
%
% in: 
%    ob         a struct or object 
%               including system parameters, such as pixel-size etc.
%    wi         cov(yi)
%    varargin   default: only compute phi' = Phi
%               'all': compute phi' = Phi + phi' = Phi- pi
% out: 
%    w0     angular-dependent weighting 
%    t      computation time
%
% Copyright 2005-4-07, Yingying Zhang, The University of Michigan

all = logical(0);
if ~isempty(varargin)
    arg = varargin{1};
    if isempty(arg)
			%	do nothing for empty arguments

	elseif ~ischar(arg) 
        error 'unknown non-string argument?'
    elseif streq(arg, 'all')
        all = logical(1);
    else
        error 'unknown string argument?'
    end
end      

% default argument
ob.Df = 0; % 3rd generation CT scanner

% image domain pixel locations
% only compute within mask

yflip = (-1); % since the top is positive in y direction % 
% x = ([1:ob.nx] - (ob.nx+1)/2)* ob.dx;
% y = ([1:ob.ny] - (ob.ny+1)/2)* ob.dx * yflip;

x = ([0:ob.nx-1] - ob.nx/2)* ob.dx;
y = ([0:ob.ny-1] - ob.ny/2)* ob.dx * yflip;

[x y] = ndgrid(x, y);

rr = sqrt(x.^2 + y.^2);
theta = atan2(y,x); % use it for our purpose since atan will have pi flip

r0 = rr(mask(:));
phi0 = theta(mask(:));

clear rr theta

% dummy variables in freq. domain
phi = Phi;

% if ~isempty(varargin) & isnumeric(varargin{1})
%     phi = varargin{1};
% end

% sinogram domain locations
ob.beta = ([0:ob.na-1]'/ob.na * ob.orbit ... 
    + ob.orbit_start) * pi / 180;
ob.s = ([-(ob.nb-1)/2:(ob.nb-1)/2]' ...
			- ob.offset_s) * ob.ds;
ob.gamma = ob.s / ob.dis_src_det;
[ss, beta] = meshgrid(ob.s,ob.beta); % meshgrid and ndgrid reverse


% loop over \phi_i
Dcf = ob.dis_src_det + ob.Df;
ratio_fcf = ob.Df/Dcf;
% w0 = zeros(ob.nx, ob.ny, ob.na);
w0 = zeros([size(r0,1) length(phi)]);

tic    
for ia = 1:length(phi)
    % -------------------------
    % for \phi' = \phi
    % -------------------------
    
    % compute new variables in (3.5.2)
    r_p = r0 .* cos(phi(ia) - phi0); % use ob.beta, not outer_sum(ob.gamma, ob.beta) 
%                                         since evaluate at \phi_p = \phi(view angle)
    s_p = Dcf.*(asin(r_p./ob.dis_src_iso) - asin(ratio_fcf.*r_p./ob.dis_src_iso));
    beta_p = phi(ia) - asin(r_p./ob.dis_src_iso);
    
    % compute Jacobian determinant (2.8.9) 
    Jacob = ob.dis_src_det ./ sqrt(ob.dis_src_iso.^2 - r_p.^2);
    
%     % find corresponding wi: nearest neighbor method
%     s_cen = (ob.nb - 1)/2 + ob.offset_s; % for nb = even % <-----(yy) plus:(to me account for the right index) or (hugo) minus? 
%     s_loc = round(s_p./ob.ds + s_cen); 
%     s_loc = min(s_loc, ob.nb);
%     s_loc = max(s_loc, 1);
%     beta_loc = round((mod(beta_p,2*pi)-ob.orbit_start)./(2*pi/ob.na)) + 1; 
%     beta_loc = min(beta_loc, ob.na);
%     beta_loc = max(beta_loc, 1);
% 
%     % find the corresponding index in column-stack wi
%     loc = 1 + (s_loc(:)-1) + (beta_loc(:)-1)*ob.nb; 

    % compute w0 at \phi' = \phi
%     w0(:,ia) =wi(loc(:)) .* Jacob;
    si = min(s_p, ob.s(end)); si = max(si, ob.s(1));
    betai = min(beta_p, ob.beta(end)); betai = max(betai, ob.beta(1));
    w0(:,ia) = interp2(ss, beta, wi.', si(:), betai(:)) .* Jacob;
    
    clear s_loc beta_loc loc 
    % -------------------------
    % for \phi' = \phi - pi
    % -------------------------
    if all
        phi_pi = phi(ia) - pi;
        r_p = r0 .* cos(phi_pi - phi0);
        s_p = Dcf.*(asin(r_p./ob.dis_src_iso) - asin(ratio_fcf.*r_p./ob.dis_src_iso)); 
        beta_p = phi_pi - asin(r_p./ob.dis_src_iso);
    
% %         Jacob_p again: not needed since same
%         Jacob_p = abs(ob.dis_src_iso*cos(s_p./ob.dis_src_det) - ...
%             ob.source_offset*sin(s_p./ob.dis_src_det))./ob.dis_src_det;
%     
%         % find corresponding wi: nearest neighbor method
%         s_loc = round(s_p./ob.ds + s_cen); 
%         s_loc = min(s_loc, ob.nb);
%         s_loc = max(s_loc, 1);
%     
%         beta_loc = round((mod(beta_p,2*pi)-ob.orbit_start)./(2*pi/ob.na)) + 1; 
%         beta_loc = min(beta_loc, ob.na);
%         beta_loc = max(beta_loc, 1);
%     
%         % find the corresponding index in column-stack wi
%         loc = 1 + (s_loc(:)-1) + (beta_loc(:)-1)*ob.nb; 
% 
%         % compute w0 at \phi' = \phi - pi and sum with \phi' = \phi
%         w0(:,ia) = w0(:,ia) + wi(loc(:)).*Jacob;

        si = min(s_p, ob.s(end)); si = max(si, ob.s(1));
        betai = min(beta_p, ob.beta(end)); betai = max(betai, ob.beta(1));
        w0(:,ia) = w0(:,ia) + interp2(ss, beta, wi.', si(:), betai(:)) .* Jacob;
        clear s_loc beta_loc loc 
    end
    
%     % angle_dependent blur
%     b0(:,:,ia) = (ob.ray_spacing./Dcf) .* sqrt(ob.dis_src_iso^2 + r0^2 + 2 .* ob.dis_src_iso .* Dcf .* cos(phi(ia)-phi0)); 
end
 
t = toc;
sprintf('angular weighting computation time %.4f', t)
    
    