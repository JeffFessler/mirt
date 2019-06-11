%	setup.m
%	run this file to set up matlab path etc.

p1=path;
p2=delete_path(p1,{'toolbox/map/'});
path(p2);
clear p1 p2

[stat_temp,hname_temp]=unix('hostname');
if strncmp(hname_temp,'ir42',4) 
  path(path, '../../../mex/v7/x86_64,gcc')                   % contains other mex files
  path(path, '../../../../../../bin/x86_64,gcc/mex')         % contains ddmex*.mexglx files
else
  path(path, '../../../mex/v7')                            % contains other mex files
  path(path, '../../../../../../bin/i686,gcc/mex')         % contains ddmex*.mexglx files
end
clear stat_temp hname_temp 

path(path, '/n/kepler/y/fessler/l/src/matlab/alg/ct')
path(path, '/n/kepler/y/fessler/l/src/matlab/alg/utilities')
path(path, '/n/kepler/y/fessler/l/src/matlab/alg/graph')
path(path, '/n/kepler/y/fessler/l/src/matlab/alg/systems')
%path(path, '/n/kepler/y/fessler/l/src/matlab/alg/wls')

%path(path, '/n/kepler/y/fessler/l/src/matlab/ge')
path(path, '../../../../ge')

path(path, '/n/kepler/y/fessler/l/src/matlab/mfiles')
path(path, '/n/kepler/y/fessler/l/src/matlab/mfiles/arch')
path(path, '/n/kepler/y/fessler/l/src/matlab/alg/penalty')
path(path, '/n/kepler/y/fessler/l/src/matlab/alg/general')
path(path, '/n/kepler/y/fessler/l/src/matlab/alg/fbp')
%path(path, '/n/kepler/y/fessler/l/src/matlab/alg/nufft')

path(path, '../../../../../../src/matlab/alg/wls')
path(path, '~someshs/matlab_utils')
double6 single % was previously double6 double


