function [mxy, mz] = steadyState_blochCim(Minit,b1,gx,gy,gz,sens,x,y,z,dt,fieldmap,roi, T1, T2,mode,dz)
% do the steady state simulation with crusher using the c-bloch simulator
% (preferred than the matrix version)
% Minit: inital magnetization: 0: [0,0,1]. Otherwise dimension should be npos*3
% b1: b1 pulse: complex number in (Gauss). dimension can be nb1 * 1 or ntime * ncoil
% gx, gy, gz: gradient in (Gauss)
% sens: normalized sensitivities matrix: must be ncoils * npos
% x,y,z: simulation range in (cm)
% dt: time duration of each point in (s)
% fieldmap: in (Hz)
% roi: mask
% T1, T2 in (s)
% mode: 0(default): just simu one TR. 1: simu steady state
% dz: resolution in z direction
% Hao Sun; Jul 22. 2012

if matlabpool('size') == 0
    matlabpool open
end



if ~exist('T1','var')
    T1 = 520/1000; %s
    T2 = 50/1000;
end

if ~exist('mode','var')
    mode = 0;
end

%% initial magnetization
if Minit == 0
    Minit = [zeros(size(fieldmap(:))),zeros(size(fieldmap(:))),ones(size(fieldmap(:)))];
end
mx0 = Minit(roi,1);
my0 = Minit(roi,2);
mz0 = Minit(roi,3);

%% mask the sens: the dimension of sens is ncoils * npos
sens_masked = [];
for ic = 1:size(sens,1)
    sens_i = sens(ic,:).';
    sens_masked(ic,:) = sens_i(roi);
end

%% mask the field map
freq = fieldmap(roi);

%% simu gradient crusher effect
% dz = 0.375; %cm 24/(192*0.3333)
nslice = 100; % the difference between using 100 and 180 is small 1e-3.
dzline = linspace(-dz/2, dz/2, nslice);
mx_all = zeros(size(freq));
my_all = zeros(size(freq));
mz_all = zeros(size(freq));

isread = (b1 == 0)&(gx == 0)&(gy == 0)&(gz == 0); % this point is read
isread_shift = [isread(2:end); 0]; % the next point is read
isSimu = (~isread) | (~isread_shift); % simulate all the non-read points and the points before a non-read point
tt = dt*(1:length(gx));
tt_m = tt(isSimu);
b1_m = b1(isSimu);
gx_m = gx(isSimu);
gy_m = gy(isSimu);
gz_m = gz(isSimu);
parfor iz = 1:length(dzline)
    ztmp = z + dzline(iz);
    [xx,yy,zz] = ndgrid(x,y,ztmp);
    [mx,my,mz] = blochCim_noPrint(b1,[gx,gy,gz],dt,T1,T2,freq,[xx(roi),yy(roi),zz(roi)],mode,sens_masked,mx0,my0,mz0);
    %     [mx,my,mz] = blochCim_noPrint(b1_m,[gx_m,gy_m,gz_m],tt_m,T1,T2,freq,[xx(roi),yy(roi),zz(roi)],mode,sens_masked,mx0,my0,mz0);
    mx_all = mx_all + mx;
    my_all = my_all + my;
    mz_all = mz_all + mz;
end

%% restore to original image dimension
mxy = (mx_all+1i*my_all)./nslice;
mxy = embed(mxy,roi);
mz = embed(mz_all,roi);
mz = mz(:)';
mxy = reshape(mxy, length(x),length(y),length(z));
% mz = reshape(mz, length(x),length(y),length(z));
mxy(~roi) = 0;
mz(~roi) = 0;