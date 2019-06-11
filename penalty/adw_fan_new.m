 function [w0, phi] = adw_fan_new(sg, ig, wi)
%function [w0, phi] = adw_fan_new(sg, ig, wi)
%
% Compute the angular-dependent weighting for fan-beam geometry.
% w0(\Phi) = w(x0, y0, \Phi)
% For fully corrected penalty:
%	w(s',\beta') * J(s') | \phi'=\Phi + ...
%	w(s',\beta') * J(s') | \phi'=\Phi-pi
% See fessler chapter 3.5.3 (?).
% Useful for regularization design and variance prediction.
%
% in:
%	sg	sino_geom()
%	ig	image_geom()
%	wi	[nb,na] var(yi)
% out:
%	w0	[np,na] angular-dependent weighting, np = sum(ig.mask(:))
%
% Copyright 2005-4-07, Yingying Zhang & Jeff Fessler, The University of Michigan

if nargin == 1 && streq(sg, 'test'), adw_fan_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

[x y] = ndgrid(ig.x, ig.y); % image domain pixel locations

r0 = sqrt(x.^2 + y.^2);
phi0 = atan2(y,x); % use it for our purpose since atan will have pi flip
% [phi0 r0] = cart2pol(x,y);
% image domain dummy variables for polar coordinates

% dummy variables in freq. domain
phi = sg.ar;

% sinogram domain locations
beta = sg.ar; % angles in radians
gamma = sg.s / sg.dsd;

% loop over \phi_i
Dcf = sg.dsd + sg.dfs;
ratio_fcf = sg.dfs / Dcf;
% w0 = zeros(ig.nx, ig.ny, sg.na);
w0 = zeros([size(r0(ig.mask),1) sg.na]);

ticker reset
for ia = 1:sg.na
	ticker(mfilename, ia, sg.na)

	% -------------------------
	% for \phi' = \phi
	% -------------------------
	
	% compute new variables in (3.5.2)
	r_p = r0 .* cos(phi(ia) - phi0); % use beta, not outer_sum(gamma, beta)
%	since evaluate at \phi_p = \phi(view angle)
%	r_p = x .* cos(phi(ia)) + y .* sin(phi(ia));
	s_p = Dcf .* (asin(r_p / sg.dso) - asin(ratio_fcf .* r_p / sg.dso));
	beta_p = phi(ia) - asin(r_p / sg.dso);
	
	% compute Jacobian determinant (2.7.9)
	Jacob_p = abs(sg.dso * cos(s_p / sg.dsd) - ...
	sg.offset * sin(s_p / sg.dsd)) / sg.dsd;
	
	% find corresponding wi: nearest neighbor method
	s_cen = sg.w;
	% caution: this looks wrong!!! (sg.nb + 1)/2 + ob.offset_s; % for nb = even % <-----(yy) plus:(to me account for the right index) or (hugo) minus?
	s_loc = round(s_p / sg.d + s_cen);
	s_loc = min(s_loc, sg.nb);
	s_loc = max(s_loc, 1);
	beta_loc = round((mod(beta_p,2*pi)-deg2rad(sg.orbit_start)) / (2*pi/sg.na)) + 1;
	beta_loc = min(beta_loc, sg.na);
	beta_loc = max(beta_loc, 1);

	% find the corresponding index in column-stack wi
	loc = 1 + (s_loc(:)-1) + (beta_loc(:)-1)*sg.nb;

	% compute w0 at \phi' = \phi
%	w0(:,:,ia) = w0(:,:,ia) + reshape(wi(loc(:)), [ig.nx ig.ny]) .* (1 ./ Jacob_p);
	wi_temp = wi(loc(:));
	wi_p = wi_temp(ig.mask);
	Jacob = 1 ./ Jacob_p(ig.mask);
	w0(:,ia) = wi_p .* Jacob;
	
	clear s_loc beta_loc loc

	% -------------------------
	% for \phi' = \phi - pi
	% -------------------------
	phi_pi = phi(ia) - pi;
	r_p = r0 .* cos(phi_pi - phi0);
%	r_p = x .* cos(phi(ia)) + y .* sin(phi(ia));
	s_p = Dcf .* (asin(r_p / sg.dso) - asin(ratio_fcf .* r_p / sg.dso));
	beta_p = phi_pi - asin(r_p / sg.dso);
	
	% Jacob_p again: not needed since same
%	Jacob_p = abs(sg.dso*cos(s_p / sg.dsd) - ...
%	sg.offset * sin(s_p / sg.dsd)) / sg.dsd;
	
	% find corresponding wi: nearest neighbor method
	s_loc = round(s_p / sg.d + s_cen);
	s_loc = min(s_loc, sg.nb);
	s_loc = max(s_loc, 1);
	
	beta_loc = round((mod(beta_p,2*pi)-deg2rad(sg.orbit_start)) / (2*pi/sg.na)) + 1;
	beta_loc = min(beta_loc, sg.na);
	beta_loc = max(beta_loc, 1);
	
	% find the corresponding index in column-stack wi
	loc = 1 + (s_loc(:)-1) + (beta_loc(:)-1)*sg.nb;

	% compute w0 at \phi' = \phi - pi and sum with \phi' = \phi
%	w0(:,:,ia) = w0(:,:,ia) + reshape(wi(loc(:)), [ig.nx ig.ny]) .* (1 ./ Jacob_p);
	wi_temp = wi(loc(:));
	wi_p = wi_temp(ig.mask);
	w0(:,ia) = w0(:,ia) + wi_p .* Jacob;

	clear s_loc beta_loc loc
	
	% angle_dependent blur
%	b0(:,:,ia) = (sg.ds / Dcf) .* sqrt(sg.dso^2 + r0^2 + 2 .* sg.dso .* Dcf .* cos(phi(ia)-phi0));
end

function adw_fan_test
down = 8;
sg = sino_geom('ge1', 'down', down);
ig = image_geom('nx', 512, 'fov', 500, 'down', down);
wi = sg.ones;
w0 = adw_fan_new(sg, ig, wi);
w0 = ig.embed(w0);
im clf, im(w0(:,:,1:10:end)), cbar
