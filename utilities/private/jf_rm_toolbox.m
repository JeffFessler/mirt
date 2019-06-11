function jf_rm_toolbox
% remove matlab toolboxes from path
% for testing my toolbox with vanilla matlab

rm aero
rm aeroblks
rm bioinfo
rm comm
rm commblks
rm compiler
rm control
rm curvefit
rm database
rm des
rm distcomp
rm dspblks
rm eml
rm emlcoder
rm filterdesign
rm finance
rm fixedpoint
rm fuzzy
rm geoweb
rm ident
rm idelink
rm images
rm instrument
%rm local
rm map
%rm matlab
rm nnet
rm optim
rm pde
rm physmod
rm rf
rm rfblks
rm robust
rm rtw
%rm shared
rm signal
rm simulink
rm sl3d
rm slcontrol
rm slvnv
rm splines
rm stateflow
rm stats
rm symbolic
rm vipblks
rm vision
rm wavelet

%path

function rm(box)
tmp = which('svd');
ii = strfind(tmp, 'matlab/matfun/@single/svd');
if numel(ii) ~= 1
	fail('bad "%s"', tmp)
end
tmp(ii:end) = '';
pre = strrep(tmp, 'built-in (', '');
%pre = '/Volumes/a2/pub/matlab/matlab-2010a.app/toolbox/';
tmp = genpath([pre box])
rmpath(tmp);
%keyboard
