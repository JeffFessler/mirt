function temp = lpass(E_map)
n = size(E_map,1);
t = E_map;
sn = 33;




vec = linspace(-pi,pi,sn)'*ones(1,sn);
gaus = exp((-vec.^2- (vec').^2)/2)/sqrt(2*pi);

span = [fliplr(fliplr(t).') flipud(t) t.';fliplr(t) t fliplr(t);t.' flipud(t) fliplr(fliplr(t).')];
convdata = conv2(span,gaus);
l = (size(convdata,2)-n)/2;
temp = convdata((l+1):(end-l) , (l+1):(end-l));
temp = temp./sum(sum(temp)).*sum(sum(E_map));