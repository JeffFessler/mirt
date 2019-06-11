% demo program: 
% This code demos:
% (1) How to use the blochCim.c code.
% (2) Compare the steady state mode of blochCim.c with for loop steady
% state simulation
clear 
% time parameters
dt = 0.1/1000;     %sec
TR = 0.3/1000;		% Sec.
Trf = dt;		% 100 us "hard" RF pulse.
T1 = 520/1000;			% Sec.
T2 = 50/1000;		% Sec.
Tpad = (TR-Trf)/2;	% Sec.
t = [Tpad Trf Tpad]; % time array

% position paramters
x = [-20:1:20];   %cm
y = [-10:10]; 
z = [2];
[xx,yy,zz] = ndgrid(x,y,z); 
freq = linspace(-100,100,length(xx(:))); 	% Hz

% generate b1
alpha = 20;		% Tip angle (Degrees).
gamma = 4258;		% Hz/G.
b1amp = pi/180*alpha/Trf/gamma/2/pi; 
b1 = [0 b1amp, 0]';	% Gauss.
b1 = b1*exp(1i*0.1*pi); % add some phase to b1
b2 = [b1,b1]; % 2 coils
b3 = [b2,b2]; % 4 coils 
b4 = [b3,b3]; % 8 coils 

% generate gradient wave
gx = zeros(size(b1)); 
gy = zeros(size(b1)); 
gz = zeros(size(b1)); 

% sensitivities
b1 = b3; % choose to simulate 4 coil case
sens = ones(4,length(xx(:))); 


%% do bloch simulation
[mxcc,mycc,mzcc] = blochCim(b1,[gx,gy,gz],dt,T1,T2,freq',[xx(:),yy(:),zz(:)],0,sens);
mcb = mxcc+1i*mycc;

%% do steady state simulation using for loop and steady state simu mode
N = T1/TR*5; 
tic
for n = 1:N
  [mxcc,mycc,mzcc] = blochCim(b1,[gx,gy,gz],dt,T1,T2,freq',[xx(:),yy(:),zz(:)],100,sens,mxcc,mycc,mzcc);
end
mtime = toc; 
mcc = mxcc + 1i*mycc;


% simu steady state using steady state mode and compare it with for loop
[mxss,myss,mzss] = blochCim(b1,[gx,gy,gz],dt,T1,T2,freq',[xx(:),yy(:),zz(:)],1,sens);
mss = mxss+1i*myss;
figure,plot(abs(mcc(:)),'o'), hold on, plot(abs(mss(:)),'r'); legend('loop', 'mode 1'); 



