% interactive tps - streamlined version
clear;
close all;

delta = 0.05; %grid size
bounds=[-0.5 1.5 -0.5 1.5];
disp_force = 1;
maxpnts = 30; %Maximum number of control points
[gridrefx, gridrefy] = meshgrid(0:delta:1, 0:delta:1); %define inital grids
[gridholx, gridholy] = meshgrid(0:delta:1, 0:delta:1); 

figure(1); 
set(1,'position',[100 300 400 400]); % position figure 1
plot(gridrefx, gridrefy, 'k-', gridrefx', gridrefy', 'k-');
title('Reference frame coordinates');xlabel('x');ylabel('y');
axis square;hold on; % plot ref intial grid setup.

figure(2);
set(2,'position',[510 300 400 400]); % position figure 2 
plot(gridholx, gridholy, 'k-', gridholx', gridholy', 'k-');
title('Homologous frame coordinates');xlabel('x');ylabel('y');
axis square; hold on; % plot hol intial grid setup.
refx=0;refy=0;
figure(1); xlabel('Point outside when done adding points');
for i=1:maxpnts,
    figure(1);[x, y, button] = ginput(1);
    if  x<0 | x>1 | y<0 | y>1 break; end;
    refx = [refx; x];
    refy = [refy; y];
    figure(1);plot(x,y,'o');
    figure(2);plot(x,y,'o');
end;

n = length(refx) - 1; %discard first zero , number of control points
refx = refx(2:n+1);
refy = refy(2:n+1);
    
holx = refx;
holy = refy; % initial set up
    
h = 2; %intial frame to work on = hol frame;
%now try to interact with points and calculate warp with the given point
while 1,
    n = length(refx);
    [L, W, A, BE] = calc_tps(refx, refy, holx, holy,1);
   
    W
    
    % now calculate deformed holmogous grids
    [size_n, size_m] = size(gridrefx);
    U = zeros(n,1);
    for i=1:size_n,
       for j=1:size_m,
           [gridholx(i,j) gridholy(i,j)] = map_points(A, W, refx, refy, gridrefx(i,j), gridrefy(i,j));
       end;
    end;
    
    %calculate forces at the control points -- forces are W coef's
    F = zeros(n, 2); % force vector in 2D
    for i=1:n,
        F(i,1) = W(1,i);
        F(i,2) = W(2,i);
    end;
    
    %calculate effective of control points 
    %how much additional displacements are possible when adding a control point
    %at specific location
    
    %calculate total displacement
    DISx=0;DISy=0;
    [DISx, DISy] = calc_displacement(A, W, refx, refy, gridrefx, gridrefy);
    
    figure(3);clf; % Display information
    set(3,'position',[200 50 600 200]);
    text(0.1, 0.9, 'Bending Energy = ' );text(0.4, 0.9, num2str(BE) );
    text(0.1, 0.8, 'Displacement x = ' );text(0.4, 0.8, num2str(DISx) );
    text(0.1, 0.7, 'Displacement y = ' );text(0.4, 0.7, num2str(DISy) );
    axis off; title('Status');
      
    %update grids and points
    figure(1);hold off;
    plot(gridrefx, gridrefy, 'k-', gridrefx', gridrefy', 'k-');
    title('Reference frame coordinates');xlabel('x');ylabel('y');
    hold on; axis square;
    for i=1:n,
        plot(refx(i),refy(i),'o');
    end; 
    figure(2);hold off;
    plot(gridholx, gridholy, 'k-', gridholx', gridholy', 'k-');
    title('Homologous frame coordinates');xlabel('x');ylabel('y');
    hold on; axis square;
    for i=1:n,
        plot(holx(i),holy(i),'o');
    end; 
    
    % update forces on control points if necessary
    force_scale = 0.5;
    figure(2);hold on;
    for i=1:n,
        fvector=[holx(i) holy(i);holx(i)+F(i,1)*force_scale  holy(i)+F(i,2)*force_scale];
        if disp_force == 1,
            plot(fvector(:,1),fvector(:,2),'r');
        end;
    end;
          
    
    if h == 1, 
        coordx = refx; coordy = refy; gridx = gridrefx; gridy = gridrefy; end;
    if h == 2, 
        coordx = holx; coordy = holy; gridx = gridholx; gridy = gridholy; end; 
    % choose which frame to work with
    
    figure(h);xlabel('Move/ Add/ Delete/ Relax/ 1-ref/ 2-hol / Force');
    [x, y, button] = ginput(1); %mouse input
    switch button
    case {'m', 'M'} %move
        best=1000000;foundidx=0;
        for i=1:n,
            r = (x - coordx(i))^2 + (y - coordy(i))^2;
            if r < best, best = r; foundidx = i; end;
        end;
        figure(h);hold off;
        plot(gridx, gridy, 'k-', gridx', gridy', 'k-');
        if h == 1, title('Reference frame coordinates'); end;
        if h == 2, title('Homologous frame coordinates'); end;
        xlabel('x');ylabel('y');hold on; axis square;
        for i=1:n,
            plot(coordx(i),coordy(i),'o');
        end;
        plot(coordx(foundidx),coordy(foundidx),'x');
        
        %now move the point
        xlabel('Move the chosen point');    
        [x, y, button] = ginput(1);
        if h==1, refx(foundidx)=x;refy(foundidx)=y; end;
        if h==2, holx(foundidx)=x;holy(foundidx)=y; end;
        
    case '1', % change into fig1, reference frame
        h = 1;
    case '2', % change into fig2, homologous frame
        h = 2;
    case {'r', 'R'} % relax
        holx = refx; holy = refy;
    case {'f', 'F'} % display forces
        if disp_force == 0, disp_force = 1; else disp_force = 0; end;
        
    case {'d', 'D'} % delete
        best=1000000;foundidx=0;
        for i=1:n,
            r = (x - coordx(i))^2 + (y - coordy(i))^2;
            if r < best, best = r; foundidx = i; end;
        end;
        
        refx(foundidx)= refx(n); refx = refx(1:n-1); 
        refy(foundidx)= refy(n); refy = refy(1:n-1); 
        holx(foundidx)= holx(n); holx = holx(1:n-1); 
        holy(foundidx)= holy(n); holy = holy(1:n-1); %move last one to deleted spot and delete last one
           
    case {'a', 'A'} % add 
        if h == 1,
            [newx, newy] = map_points(A, W, refx, refy, x, y);
            disp('add at ref side'); W
            
            refx=[refx; x]; refy=[refy; y];
            holx=[holx; newx]; holy=[holy; newy];
        end;
        if h == 2, % need to calculate inverse TPS 
            [L, W, A, BE]=calc_tps(holx, holy, refx, refy, 0);
            disp('add at hol side'); W
            
            [newx, newy] = map_points(A, W, holx, holy, x, y);
            holx=[holx; x]; holy=[holy; y];
            refx=[refx; newx]; refy=[refy; newy];
        end;
        
    end; %end of switch
    
end; %end of while
