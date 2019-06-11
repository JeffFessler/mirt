function [mxy, mz] = parallel_blochCim(Minit,pulses,gx,gy,gz,sens,x,y,z,dt,fieldmap,roi, T1, T2,mode)
% do the parallel transmit Bloch simulation using the c-bloch simulator
% Minit: inital magnetization: 0: [0,0,1]. Otherwise dimension should be npos*3
% pulses: b1 pulse: complex number in (Tesla). dimension can be nb1 * 1 or ntime * ncoil
% gx, gy, gz: gradient in (Gauss)
% sens: normalized sensitivities matrix: must be ncoils * npos
% x,y,z: simulation range in (cm)
% dt: time duration of each point in (s)
% fieldmap: in (Hz)
% roi: mask
% T1, T2 in (s)
% mode: 0(default): just simu one TR. 1: simu steady state
% Hao Sun; Jul 22. 2012

roi = logical(roi); 
% pulses = pulses.*1e4; %convert to gauss 
freq = fieldmap(roi);

if ~isvar('T1') 
T1 = 520/1000; %s 
T2 = 50/1000; 
end

if ~isvar('mode')
    mode = 0; 
end 

[xx,yy,zz] = ndgrid(x,y,z); 

if Minit == 0
    Minit = [zeros(size(fieldmap(:))),zeros(size(fieldmap(:))),ones(size(fieldmap(:)))]; 
end
mx0 = Minit(roi,1); 
my0 = Minit(roi,2); 
mz0 = Minit(roi,3); 


sens_masked = []; 
% mask the sens: the dimension of sens is ncoils * npos 
for ic = 1:size(sens,1)
    sens_i = sens(ic,:).';
    sens_masked(ic,:) = sens_i(roi);
end
[mx,my,mz] = blochCim(pulses,[gx,gy,gz],dt,T1,T2,freq,[xx(roi),yy(roi),zz(roi)],mode,sens_masked,mx0,my0,mz0);
 
mxy = mx+1i*my;  
mxy = embed(mxy,roi);
mz = embed(mz,roi);
mz = mz(:)'; 

mxy = reshape(mxy, length(x),length(y),length(z));
% mz = reshape(mz, length(x),length(y),length(z));
M(~roi) = 0;
mz(~roi) = 0;