 function es = ellipse_motion(ells, varargin)
%function es = ellipse_motion(ells, varargin)
%
% create strum object for simple ellipse motion
% to pass to ellipse_sino()
%
% in:
%	ells	[ne,6]		[centx centy radx rady angle_degrees amplitude]
%
% options:
%	'type'			'linear' (default)
%				'none' no motion, for testing
%	'ellend'		[ne,6] end ellipse for linear motion
%
% out:
%	es			strum with methods
%		es.ell(ie,na)	returns [na,6]
%
% Copyright 2007-10-12, Jeff Fessler, The University of Michigan

if ~nargin, help(mfilename), error(mfilename), end
if streq(ells, 'test'), fbp2_motion_example, return, end

% defaults
arg.type = 'linear';
arg.ellend = [];
arg = vararg_pair(arg, varargin);

switch arg.type
case 'none'
	es.ells = ells;
	es.ne = size(ells,1);
	fun = @(es, ie, na) repmat(es.ells(ie,:), [na 1]);
	es = strum(es, {'ell', fun, '(ie,na)'});

case 'linear'
	es.ell0 = ells;
	es.ell1 = arg.ellend;
	es.ne = size(ells,1);
	if (es.ne ~= size(es.ell1,1)), fail 'ellend size mismatch', end
	es = strum(es, {'ell', @ellipse_motion_linear, '(ie,na)'});

otherwise
	fail('motion "%s" not done', arg.type)
end

function ell = ellipse_motion_linear(es, ie, na)
frac = linspace(0, 1, na)';
ell0 = es.ell0(ie,:);
ell1 = es.ell1(ie,:);
ell = (1-frac) * ell0 + frac * ell1; % [na,6] 
