% This code demos:
% (0) How to use the blochCim.c code.
% (1) Compare c simulation with MATLAB simulation. 
% (2) Compare the steady state mode of blochCim.c with for loop steady
% state simulation
% setup; 
clear 

%	--- Setup parameters ---
TR = 0.3/1000;		% Sec.
Trf = 0.1/1000;		% 100 us "hard" RF pulse.
alpha = 20;		% Degrees.
gamma = 4258;		% Hz/G.
T1 = 520/1000;			% Sec.
T2 = 50/1000;		% Sec.
Tpad = (TR-Trf)/2;	% Sec.
doplot = 1;
x = [-20:20];   %cm
y = [-10:10]; 
z = [2];
[xx,yy,zz] = ndgrid(x,y,z); 
freq = linspace(-100,100,length(xx(:))); 	% Hz

t = [Tpad Trf Tpad];
dt = Trf;
% random b1
b1 = [0 pi/180*alpha/Trf/gamma/2/pi, 0]';	% Gauss.
b1 = b1*exp(1i*0.1*pi); % single coil 
b2 = [b1,randn(size(b1))]; % 2 coils
b3 = [b2,randn(size(b2))]; % 4 coils
b4 = [b3,b3]; % 8 coils

% random gradient
gx = randn(size(b1)); 
gy = randn(size(b1)); 
gz = randn(size(b1)); 

% random sensitivity
sens = randn(1,length(xx(:))); 
sens = sens + 1i*randn(size(sens)); 

%% compare blochCim with old matlab simulator
[mxcc,mycc,mzcc] = blochCim(b1,[gx,gy,gz],dt,T1,T2,freq',[xx(:),yy(:),zz(:)],0,sens);
mcc = mxcc+1i*mycc;
roi_simu = ones(size(xx)); 
[m1,mz1] = parallel_blochsim_field_roi(0,b1/10000,gx,gy,gz,sens, x(:),y(:),z(:),dt,freq,0,roi_simu);
figure,plot(abs(m1(:)-mcc)); title('compare c simu with matlab(plot difference)');  



%% compare steady state simulation using loop of blochCim.c blochsim2.m and
% steady state mode of blochCim.c 
N = round(T1/TR*5); 
m1_all = zeros(N,1);  
mcc_all = zeros(N,1);  
% steady state using for loop of matlab bloch simulation code
tic
for n = 1:N
  minit1 = [real(m1(:)), imag(m1(:)), mz1'];  
  [m1,mz1] = parallel_blochsim_field_roi(minit1,b1/10000,gx,gy,gz,sens, x(:),y(:),z(:),dt,freq,0,roi_simu);
  m1_all(n) = m1(round(length(m1)/2)); 
end
matLoopTime = toc; 

% steady state using for loop of blochCim
tic
for n = 1:N
  [mxcc,mycc,mzcc] = blochCim(b1,[gx,gy,gz],dt,T1,T2,freq',[xx(:),yy(:),zz(:)],100,sens,mxcc,mycc,mzcc);
  mcc_all(n) = mxcc(round(length(m1)/2)) + 1i*mycc(round(length(m1)/2)); 
end
cLoopTime = toc; 
mcc = mxcc + 1i*mycc; 
figure,plot(abs(m1_all)), hold on, plot(abs(mcc_all),'r'); title('time curve'); legend('Matlab', 'C'); 

% steady state using steady state simulation mode of blochCim
tic
[mxss,myss,mzss] = blochCim(b1,[gx,gy,gz],dt,T1,T2,freq',[xx(:),yy(:),zz(:)],1,sens);
ssTime = toc
mss = mxss+1i*myss;
figure,plot(abs(mcc(:))), hold on, plot(abs(mss(:)),'r'); legend('loop', 'mode 1'); 

% compare the run time
matLoopTime
cLoopTime
ssTime



