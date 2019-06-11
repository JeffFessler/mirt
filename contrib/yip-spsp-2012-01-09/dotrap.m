function trap=dotrap(area,gmax,dgdt,dt);
% dotrap(area,gmax,dgdt,dt);
%   area = pulse area in (g sec)/cm
%   gmax = max gradient in g/cm
%   dgdt = max slew in g/cm/sec
%   dt   = sample time in sec

ramppts=ceil(gmax/dgdt/dt);
triareamax = ramppts*dt*gmax;
if triareamax > area
  % triangle pulse
  newgmax = sqrt(area*dgdt);
  ramppts=ceil(newgmax/dgdt/dt);
  pulse = [(0:ramppts)./ramppts (ramppts:-1:0)./ramppts];
else
  % trap pulse
  nflat = ceil((area-triareamax)/gmax/dt/2)*2;
  pulse = [(0:ramppts)./ramppts ones([1 nflat]) (ramppts:-1:0)./ramppts];
end


trap = pulse*(area/(sum(pulse)*dt));
%peakslew = max(diff(trap)./dt/100)

