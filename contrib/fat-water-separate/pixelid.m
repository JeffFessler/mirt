function [E_map,w,f,costf,iter,value] = pixelid(data,F_fw,te,varargin)

% in
%   data[n,n,nsets]     source data
%   F_fw = fat water shift (omega)
%   te [nsets]   time
%options
%   niter  # of iterations
%   x       [array or one value]
%   y       [array or one value]
%           x,y are array or one value of object field map 
%           If x = 100, y = 200 then we choose just (100,200)th pixel of field map 
%           If x = 1:100, y = 200:202 then we choose (1:100,200:202) of
%           field map ( = estimate 100*103 pixels)
%   init  [n,n] Initial guess of field map
%   cost    If cost = 0 : Do not calculate cost function (Default)
%           If cost = 1 : Calculate cost function
%% initialize


n = size(data,1); % length of the image

%defaults
arg.niter = 30;
arg.x = 1:512;
arg.y = 1:512;  % Just choose (100,100)th pixel of field map
arg.init = zeros(n,n);
arg.cost = 0;
arg = vararg_pair(arg, varargin);

c = cos(F_fw*te).'; 
d = sin(F_fw*te).'; % for Appendix A and B
A1 =  [1 0 c(1) -d(1); 1 0 c(2) -d(2); 1 0 c(3) -d(3); 0 1 d(1) c(1); 0 1 d(2) c(2); 0 1 d(3) c(3)];
A2 = inv(A1'*A1)*A1';  % for Appendix A
map = arg.init; % initial field map
rgmap = arg.init;
A3(:,1) = ones(length(te),1);
A3(:,2) = exp(j*F_fw.*te); % for calculate cost function
if arg.cost == 1
    length_cost = size(arg.x,2)*size(arg.y,2);
    costf = zeros(length_cost*2,arg.niter*2);
end
iteration = 0;

%% pixel independent method

for arr = arg.x
    if (rem(arr,100) == 0)
         disp(['# of iteration of RE = ' num2str(arr)])
    end

    for col = arg.y
        
        iteration = iteration+1;
        init_guess = arg.init(arr,col);
        [est_map,cost,iter] = regionpixel(data,arr,col,init_guess,te,c,d,A1,A2,A3,arg.niter,arg.cost);
        if arg.cost == 1
                costf((iteration*2-1):iteration*2,:) = cost; 
                % cost(2*k-1,:) = field map values(record of estimated field map values in each iteration)
                % cost(2*k,:) = cost values 
        end
        
        % My constraint
        if est_map > 90
            est_map = 90;
        elseif est_map <-80
            est_map = -80;
        end
        map((col-1)*n+arr) = est_map;
    end
end
if arg.cost == 1
    %% plot the cost functions (plot only the first pixel's cost function)
    % just plot costfunction(field map,water^,fat^),water^,fat^ is estimated
    % water, fat based on field map
    figure(5);
    plot(costf(1,(1:arg.niter)*2-1),costf(2,(1:arg.niter)*2-1),'o');title('cost function');
    hold on
    plot(costf(1,iter*2-1),costf(2,iter*2-1),'o','color','red'); % final estimated field map
    hold off

    %% plot the whole cost function (from -150Hz - 150Hz)
    a = linspace(-1200,1200,10000);
    orig_s = reshape(data(arg.x(1),arg.y(1),:),3,1);
    for k = 1:10000
        s =  orig_s.*(exp(-j*2*pi*a(k)*te).'); % shat (m is freq)
        shat = [real(s);imag(s)];
        rowhat = A2*shat;
        water = rowhat(1) + j*rowhat(2);
        fat = rowhat(3) + j*rowhat(4);
        temp(1,k) = norm(orig_s - A3*[water;fat].*(exp(j*2*pi*a(k)*te).'));
    end
    figure(6);
    plot(a,temp,'o');
    w = 0;
    f = 0;
    E_map = map;
    value = est_map;
else
    costf = 0;
    value = 0;
    E_map = zeros(n,n);
    E_map(find(map>0)) = map(find(map>0));
    E_map(find(map<0)) = map(find(map<0));
    E_map = lpass(E_map);
    E_map = E_map*2*pi;
    rgmap = rgmap*2*pi;
   [w f] = restore(data,E_map,F_fw,te);
end
 