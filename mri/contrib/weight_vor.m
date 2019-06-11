function area=weight_vor(kx,ky,nl)
% function area=weight_vor(kx,ky,nl)
% nl is the number of interleaves
% For calculating sampling density function
%   for spiral trajectory
%  Sangwoo Lee and Brad Sutton
%   University of Michigan
%
% Updated for Matlab R14 9/13/05 BPS

  
[V,C]=voronoin([kx,ky],{'Qbb','QJ'});

%offset = 20;
offset = 10;

  area=zeros(1,length(kx));
  for ii=1:length(kx)
    xx=[];yy=[];s=[];
    s=C{ii}(1:end);
    xx=V(s,1);
    yy=V(s,2);
    area(ii)=polyarea(xx,yy);
  end;



%Now fix edges
%keyboard
for jj = 1:nl

end_j = jj*(length(kx))/nl;

   ii = end_j;
    %while (isnan(mean(area(ii-offset:ii))) | (sum(abs(area(ii-offset:ii)-mean(area(ii-offset:ii))))>(1/nl)))
   while (isnan(mean(area(ii-offset:ii))) | (sum((abs(area(ii-offset:ii))>(offset*3)))))
     ii = ii-10;
   end

   area(ii+1:end_j) = mean(area(ii-offset:ii));

end


area = area.';


