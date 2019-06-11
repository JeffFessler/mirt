   function [mean_pos,ma_pos] = mean_posit_track(t,s,tau,train_len)
%| function [mean_pos,ma_pos] = mean_posit_track(t,s,tau,train_len)
%| INPUT
%|	t: time indices
%|	s: position sequence (should be same length as t)
%|	tau: delay parameter,
%|		set to be the equivalence index to about 0.5sec for breathing.
%|	train_Len: training length (starting up stage),
%|		approximately one period's equivalence
%|
%| OUTPUT
%|	mean_pos: a N*2 matrix,
%|		for practical purposes, take the first column mean_pos(:,1)
%|	ma_pos: naive moving average.
%|
%| Dan Ruan, University of Michigan, 2007
%|
%| Based on the paper by Dan Ruan, J A Fessler, James M Balter
%| Mean position tracking of respiratory motion
%| Med. Phys. 35(2):782-92, Feb. 2008
%| doi 10.1118/1.2825616

if nargin < 1, help(mfilename), error(mfilename), end
if streq(t, 'test'), self_test, return, end

s_complete = s;
s_d = s(1:end-tau);
s = s(tau+1:end);
for tt = 1:length(t)-train_len-tau;
        s_temp = s(tt:tt+train_len);
        sd_temp = s_d(tt:tt+train_len);
        params = fitellipse(s_temp,sd_temp);
        mean_pos(tt,:) =  [params(1),params(2)];
        ma_pos(tt) = mean(s(tt:tt+train_len));
end

plot(t,s_complete,'linewidth',2);
hold on
plot(t(1+tau+train_len:end),ma_pos,'g:','linewidth',2);
plot(t(1+tau+train_len:end),mean_pos(:,1),'r--','linewidth',2);
hold off
set(gca,'DataAspectRatio',[5 5 1]);
set(gca,'fontsize',14,'FontWeight','demi');
h = axis;
axis(h + [0 0 -1 1]);
xlabel('time (second)');
ylabel('displacement');
legend('Observation','Estimation with MA','Estimate with Ellipse Center','Location','SouthWest');


%
function a = fitellipse(X,Y)
% fit ellipse with sample locations in X,Y
% normalize data
mx = mean(X);
my = mean(Y);
sx = (max(X)-min(X))/2;
sy = (max(Y)-min(Y))/2; 

x = (X-mx)/sx;
y = (Y-my)/sy;

% Force to column vectors
x = x(:);
y = y(:);

% Build design matrix
D = [ x.*x  x.*y  y.*y  x  y  ones(size(x)) ];

% Build scatter matrix
S = D'*D;

% Build 6x6 constraint matrix
C(6,6) = 0; C(1,3) = -2; C(2,2) = 1; C(3,1) = -2;

  
  % Break into blocks
  tmpA = S(1:3,1:3); 
  tmpB = S(1:3,4:6); 
  tmpC = S(4:6,4:6); 
  tmpD = C(1:3,1:3);
  tmpE = inv(tmpC)*tmpB';
  [evec_x, eval_x] = eig(inv(tmpD) * (tmpA - tmpB*tmpE));
  
  % Find the positive (as det(tmpD) < 0) eigenvalue
  I = find(real(diag(eval_x)) < 1e-8 & ~isinf(diag(eval_x)));
  
  % Extract eigenvector corresponding to negative eigenvalue
  A = real(evec_x(:,I));
  
  % Recover the bottom half...
  evec_y = -tmpE * A;
  A = [A; evec_y];



% unnormalize
par = [
  A(1)*sy*sy,   ...
      A(2)*sx*sy,   ...
      A(3)*sx*sx,   ...
      -2*A(1)*sy*sy*mx - A(2)*sx*sy*my + A(4)*sx*sy*sy,   ...
      -A(2)*sx*sy*mx - 2*A(3)*sx*sx*my + A(5)*sx*sx*sy,   ...
      A(1)*sy*sy*mx*mx + A(2)*sx*sy*mx*my + A(3)*sx*sx*my*my   ...
      - A(4)*sx*sy*sy*mx - A(5)*sx*sx*sy*my   ...
      + A(6)*sx*sx*sy*sy   ...
      ]';

% Convert to geometric radii, and centers

thetarad = 0.5*atan2(par(2),par(1) - par(3));
cost = cos(thetarad);
sint = sin(thetarad);
sin_squared = sint.*sint;
cos_squared = cost.*cost;
cos_sin = sint .* cost;

Ao = par(6);
Au =   par(4) .* cost + par(5) .* sint;
Av = - par(4) .* sint + par(5) .* cost;
Auu = par(1) .* cos_squared + par(3) .* sin_squared + par(2) .* cos_sin;
Avv = par(1) .* sin_squared + par(3) .* cos_squared - par(2) .* cos_sin;

% ROTATED = [Ao Au Av Auu Avv]

tuCentre = - Au./(2.*Auu);
tvCentre = - Av./(2.*Avv);
wCentre = Ao - Auu.*tuCentre.*tuCentre - Avv.*tvCentre.*tvCentre;

uCentre = tuCentre .* cost - tvCentre .* sint;
vCentre = tuCentre .* sint + tvCentre .* cost;

Ru = -wCentre./Auu;
Rv = -wCentre./Avv;

Ru = sqrt(abs(Ru)).*sign(Ru);
Rv = sqrt(abs(Rv)).*sign(Rv);

a = [uCentre, vCentre, Ru, Rv, thetarad];


% self_test
function self_test
t = linspace(0, 60, 601);
dt = t(2) - t(1);
s = (8 + cos(2*pi/6*t)).^2 - 8^2;
tau = round(0.5 / dt);
train_len = round(5.1 / dt);

[mean_pos ma_pos] = mean_posit_track(t, s, tau, train_len);

prompt
tt = t(tau+train_len+1:end);
clf, plot(t, s, '-', tt, mean_pos(:,1), '--', tt, ma_pos, ':', 'linewidth', 2)
legend('signal', 'proposed', 'moving avg')
xlabel t
