 function msino = par2moj_rebin(psino, varargin)
%function msino = par2moj_rebin(psino, varargin)
%
% rebin parallel-beam sinogram into mojette sinogram
%
% in
%	psino	[nr,na]		parallel-beam sinogram
% option
%	(many, see arg.* below)
% out
%	msino	[nm,na]		mojette sinogram
%
% Copyright 2005-12-7, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(psino, 'test'), par2moj_rebin_test, return, end
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

msino = par2moj_rebin_do(psino, arg.dr, arg.offset_r, ...
	arg.nt, arg.dx, arg.offset_t, ...
	arg.interp, ang);

%
% par2moj_rebin_do()
%
function sino = par2moj_rebin_do(sino, dr, offset_r, ...
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


%
% test par2moj_rebin()
%
function par2moj_rebin_test

down = 4;
nr = 1000/down;
na_p = 800/down;
ell = [30 0 160 150 0 1];
dr = 0.5*down;
offset_r = 0.25;      % quarter detector
orbit = 180;

nt = 1120/down;
na_m = na_p;
dx = down;
offset_t = 0.25;

dx = dr;

% analytical sinograms
oversample = 4;
psino = ellipse_sino(ell, 'ds', dr, 'nb', nr, 'na', na_p, 'orbit', orbit, ...
	'offset_s', offset_r, 'oversample', oversample, 'mojette', 0);

msino = ellipse_sino(ell, 'ds', [], 'nb', nt, 'na', na_m, 'orbit', orbit, ...
	'offset_s', offset_t, 'oversample', oversample, 'mojette', dx);

% rebinned sinograms
%profile on
cpu etic
msino2 = par2moj_rebin(psino, 'dr', dr, 'offset_r', offset_r, ...
	'orbit', orbit, ...
	'nt', nt, 'dx', dx, 'offset_t', offset_t);
cpu etoc 'par2moj rebin time'

%	'interp', {'order', 1, 'ending', 'periodic', 'mex', 0}, ...
%profile report, return

%cpu etic
psino2 = zeros(size(psino));
%psino2 = moj2par_rebin(fsino, 'ds', ds, 'offset_s', offset_s, ...
%	'nr', nr, 'nphi', nphi, 'dr', dr, 'offset_r', offset_r);
%cpu etoc 'moj2par rebin time'

clf, pl=230;
im(pl+1, psino, 'psino'), cbar
im(pl+4, msino, 'msino'), cbar
im(pl+2, psino2, 'psino2'), cbar
im(pl+5, msino2, 'msino2'), cbar
im(pl+3, psino-psino2, 'error'), cbar
im(pl+6, msino-msino2, 'error'), cbar

max_percent_diff(msino, msino2)
max_percent_diff(psino, psino2)
