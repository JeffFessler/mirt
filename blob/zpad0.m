 function xpad = zpad0(x, nxp, nyp, nzp)
%function xpad = zpad0(x, nxp, nyp, nzp)
% zero pad an input signal x symmetrically around "0" (image center)
% by A. Yendiki

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(x, 'test'), zpad0_test, return, end

if nargin == 2
	xpad = ir_pad_into_center(x, nxp);
elseif nargin == 3
	xpad = ir_pad_into_center(x, [nxp nyp]);
elseif nargin == 4
	xpad = ir_pad_into_center(x, [nxp nyp nzp]);
else
	error 'Too many arguments'
end

function zpad0_test
zpad0([1 2 1], 7)
zpad0([1 2 1], 7, 3)
