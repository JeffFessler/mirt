function [M,mz]=parallel_blochsim_field_roi(Minit,pulses,gx,gy,gz,sens,x,y,z,dt,fieldmap,fieldmap_start,roi, T1, T2)
%set up bz with gradients
% T1, T2 (sec)
% pulses = pulses/1e4; %convert from Gauss to Tesla
if ~exist('T1', 'var')
    T1 = 520/1000;
    T2 = 50/1000;
end
T1 = T1*1000; 
T2 = T2*1000; %convert to (msec)

roi = logical(roi); 
[X,Y,Z] = ndgrid(x,y,z);
X = X(roi(:))';
Y = Y(roi(:))';
Z = Z(roi)';
fieldmap = fieldmap(roi);
sens = sens(:,roi);


gam=26751;                   %rad/sec/g
gambar = gam/2/pi;           % Hz/g
th_grad = fieldmap/gambar;
th_grad = th_grad(:)';
%bz = (gx*X+gy*Y+gz*Z+ones(size(gx))*th_grad)*10^(-4);
bz = gx*X;
bz = bz+gy*Y;
bz = bz+gz*Z;
eff = [zeros(fieldmap_start,1);ones(length(gz)-fieldmap_start,1)];
bz = bz+eff*th_grad;
bz = bz*10^(-4);


%set up bx and by with pulses and sens
bx = 0;
by = 0;
for i = 1:size(sens,1)
    pulse_i = squeeze(pulses(:,i));
    sens_i = squeeze(sens(i,:,:));
    sens_i = sens_i(:);
    b_xy = pulse_i*(sens_i.');
    bx = bx+real(b_xy);
    by = by+imag(b_xy);
end

T1 = T1*ones(length(X),1);
T2 = T2*ones(length(X),1);
if (Minit == 0)
    Mi = [0;0;1]*ones(1,length(X));
else
        Mi = Minit(roi(:),:).';
end
[mx,my,mz] = blochsim2(Mi, bx, by, bz, T1, T2, dt*10^(3));

m_trav = mx(end,:)+1i*my(end,:);
mz = mz(end,:);


m_trav = embed(m_trav.',roi);
mz = embed(mz.',roi);
mz = mz(:).';

M = reshape(m_trav, length(x),length(y),length(z));
% mz = reshape(mz, length(x),length(y),length(z));
M(~roi) = 0;
mz(~roi) = 0;
