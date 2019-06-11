function [map w f,temp] = region(data,map,F_fw,te)
%% Initialize

map = map/2/pi;
n = size(data,1);
temp = zeros(512,512);
% A is a sorted matrix [1 arr col]
A = zeros(n^2,3);
A(:,1) = ones(n^2,1);
B = (ones(n,1))*(1:n);
A(:,2) = reshape(B',n^2,1);
A(:,3) = reshape(B,n^2,1); 

afterest = zeros(n,n); 
% afterest(arr,col) = 1 means we finish estimated fieldmap(arr,col)

% sarr = 164;
% scol = 96;
% earr = 195;
% ecol = 127;
sarr = 94;
scol = 108;
earr = 94+31;
ecol = 108+31;

afterest(sarr:earr,scol:ecol) = 1;
map = map.*afterest;
check = 0;
check1 = 0;
check2 = 0;
check3 = 0;
check4 = 0;
iter = 0;
map = map(:);

% For pixel independent estimation
c = cos(F_fw*te).'; 
d = sin(F_fw*te).';
A1 =  [1 0 c(1) -d(1); 1 0 c(2) -d(2); 1 0 c(3) -d(3); 0 1 d(1) c(1); 0 1 d(2) c(2); 0 1 d(3) c(3)];
A2 = inv(A1'*A1)*A1'; 
niter = 30; % # of iterstion
cal = 0;

%% Region growing method

while (check == 0)

    iter = iter +1;
    if rem(iter,50) == 0
        sprintf('# of iteration of RE = %d',iter)
    end
    % map   [n,n]   initial field map 
    
    if sarr>1 % direction of region growing is down
%         sprintf('The number of iteration = %d - %d',iter,1)
        sarr = sarr-1;  
        col = scol;
        while(col < ecol+1) 
           
            [init_map] = region_sub(map,A,afterest,sarr,col);
            if init_map > 100
                init_map = 100;
            elseif init_map<-100
                init_map = -100;
            end
            temp((col-1)*n+sarr) = init_map;    
        
            % init_map is new initial guess for field map(sarr,col) by region
            % growing method
            afterest((col-1)*n+sarr) = 1;
            [new_map] = regionpixel(data,sarr,col,init_map,te,c,d,A1,A2,0,niter,cal);
            % new_map is final estimated field map(sarr,col) by pixel 
            % independent method
            if new_map > 100
                new_map = 100;
            elseif new_map<-100
                new_map = -100;
            end
            map((col-1)*n+sarr) = new_map;
            
            col = col+1;
        end
    else
        check1 = 1;
%         sprintf('Finish %d',1)
    end
    
    if ecol<n % direction of region growing is right
%         sprintf('The number of iteration = %d - %d',iter,2)
        ecol = ecol+1;
        arr = sarr;
        while(arr < earr+1)
            [init_map] = region_sub(map,A,afterest,arr,ecol);
            if init_map > 100
                init_map = 100;
            elseif init_map<-100
                init_map = -100;
            end
            temp((ecol-1)*n+arr) = init_map;
    
            afterest((ecol-1)*n+arr) = 1;
            
            [new_map] = regionpixel(data,arr,ecol,init_map,te,c,d,A1,A2,0,niter,cal);
            if new_map > 100
                new_map = 100;
            elseif new_map<-100
                new_map = -100;
            end
            map((ecol-1)*n+arr) = new_map;    
          
            arr= arr+1;
        end
    else
        check2 = 1; 
%         sprintf('Finish %d',2)
    end

    if earr<n % direction of region growing is up
%         sprintf('The number of iteration = %d - %d',iter,3)
        earr = earr+1;
        col = ecol;
        while(col > scol-1)
            [init_map] = region_sub(map,A,afterest,earr,col);
            if init_map > 100
                init_map = 100;
            elseif init_map<-100
                init_map = -100;
            end
            temp((col-1)*n+earr) = init_map;

            afterest((col-1)*n+earr) = 1;
            [new_map] = regionpixel(data,earr,col,init_map,te,c,d,A1,A2,0,niter,cal);
            if new_map > 100
                new_map = 100;
            elseif new_map<-100
                new_map = -100;
            end
            map((col-1)*n+earr) = new_map;    
            
            col = col-1;
        end
    else
        check3 = 1;
%         sprintf('Finish %d',3)
    end
    
    if scol>1   % direction of region growing is left
%         sprintf('The number of iteration = %d - %d',iter,4)
        scol = scol-1;
        arr = earr;
        while(arr > sarr-1)
            [init_map] = region_sub(map,A,afterest,arr,scol);
            if init_map > 100
                init_map = 100;
            elseif init_map<-100
                init_map = -100;
            end
            temp((scol-1)*n+arr) = init_map;
   
            afterest((scol-1)*n+arr) = 1;
            [new_map] = regionpixel(data,arr,scol,init_map,te,c,d,A1,A2,0,niter,cal);
            if new_map > 100
                new_map = 100;
            elseif new_map<-100
                new_map = -100;
            end
            map((scol-1)*n+arr) = new_map;    
            
            arr = arr-1;
        end
    else
        check4 = 1;
%         sprintf('Finish %d',4)
    end
    check = check1*check2*check3*check4;
    
end

%% Plot the result

tempm = reshape(map,n,n)*2*pi;

map = zeros(512,512);
map(find(tempm>0)) = tempm(find(tempm>0));
map(find(tempm<0)) = tempm(find(tempm<0));
map(find(map>110*2*pi)) = 110*2*pi;
map(find(map<-110*2*pi)) = -110*2*pi;
[w f] = restore(data,map,F_fw,te); 

% % w = water image, f = fat image based on the estimated field map
% figure(1)
% im(map);title('field map');
% figure(2)
% im(w);title('water');
% figure(3)
% im(f);title('fat');






