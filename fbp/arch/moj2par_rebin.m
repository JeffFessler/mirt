 function psino = moj2par_rebin(msino, varargin)
%function psino = moj2par_rebin(msino, varargin)
%
% rebin mojette sinogram into parallel-beam sinogram
%
% in
%	psino	[nr,na]		parallel-beam sinogram
% option
%	(many, see arg.* below)
% out
%	msino	[nm,na]		mojette sinogram
%
% Copyright 2005-12-7, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(psino, 'test'), par2moj_rebin test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% parallel-beam defaults
arg.dr = 1;
arg.offset_r = 0;
% mojette defaults
arg.nt = [];
arg.dx = 1;
arg.offset_t = 0;
% general defaults
arg.orbit = 180;
arg.orbit_start = 0;
arg.interp = {'order', 3, 'ending', 'zero'};

arg = vararg_pair(arg, varargin);

if isempty(arg.nt), arg.nt = size(psino,1); end
arg

na = size(psino,2);
ang = deg2rad(arg.orbit_start + arg.orbit * [0:na-1] / na);

msino = moj2par_rebin_do(psino, arg.dr, arg.offset_r, ...
	arg.nt, arg.dx, arg.offset_t, ...
	arg.interp, ang);

%
% moj2par_rebin_do()
%
function sino = moj2par_rebin_do(sino, dr, offset_r, ...
	nt, dx, offset_t, interp, ang)

% radial interpolation

[nr na] = size(sino);
wr = (nr-1)/2 + offset_r;
wt = (nt-1)/2 + offset_t;
t = ([0:nt-1]' - wt) * dx;
r = t * max(abs(cos(ang)), abs(sin(ang))); % [nt,na]
r_int = r / dr + wr;

%sino = bspline_1d_interp(sino, r_int, interp{:});
sino = interp2([0:nr-1], 1:na, sino', r_int, repmat(1:na, [nt 1]));
