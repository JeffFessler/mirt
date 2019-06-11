% Execute a toy pulse design problem to illustrate fast small- and
% large-tip-angle parallel RF pulse design.
% The key objects used here are Gmri_SENSE (small-tip) and
% Gocrf (large-tip).
%
% Copyright 2009-3-17, Will Grissom, Stanford University
% (based on work he initiated at University of Michigan)
% modified slightly by JF

tipangle = 180; % pulse flip angle
traj = 'ep';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load sensitivities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensitivities Courtesy of Steve Wright: 10th ISMRM, 2002, p. 854
if ~isvar('sens')
	load fdtdsens
	mask = sum(abs(sens),3) > 0.1*max(col(sum(abs(sens),3)));
	dim = size(sens,1); % dimension of square x-y grid
	Nc = size(sens,3); % number of Tx coils
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the desired flip angle pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvar('md')
	FOV = 24; % cm
	dimrf = 32;
	[xx, yy]=ndgrid(-FOV/2:FOV/dim:FOV/2-FOV/dim);
	d = double((abs(xx)<=10/2)&(abs(yy)<=5/2)); % the desired pattern

	% blur the desired pattern to reduce gibbs ringing
	d = conv2(d,exp(-(xx.^2+yy.^2)/0.5),'same');
	d = d./max(abs(d(:)));

	% magnetization-domain desired patterns
	md = sin(d*tipangle*pi/180);
	mzd = cos(d*tipangle*pi/180);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a k-space trajectory and its gradients
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load k-space trajectories
if ~isvar('k')
	load exckspace

	if strcmp(traj,'spiral')
	k = ksp;g = gsp;NN = NNsp;
	else
	k = kep;g = gep;NN = NNep;
	end

	Nt = size(k,1); % number of samples in pulses
	dt = 4e-6; % seconds, RF and gradient sampling period
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the small-tip object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvar('Gsml')
	% field-map params
	tt = 0:dt:(Nt-1)*dt;tt = tt-(Nt-1)*dt;L = 4;fmap = []; % Hz
	% nufft params
	J = 6;K = 2*dim;
	nufft_args = {[dim dim],[J J],[K K],[dim dim]/2,'minmax:kb'};
	gambar = 4257;	% gamma/2pi in Hz/g
	gam = gambar*2*pi; % gamma in radians/g
	% trick: the system matrix is just the transpose of a SENSE image recon matrix!
	Gsml = Gmri_SENSE(k,true(dim),'fov',[FOV FOV],'basis',{'dirac'}, ...
		'nufft',nufft_args,'exact',0, ...
		'sens',conj(reshape(sens,[dim*dim Nc]))*(-1i*gam*dt), ...
		'ti',-tt,'L',L,'zmap',2*pi*1i*fmap)';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design the initial pulse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvar('B1')
	% get the tikhonov penalty vector
	beta = 0.01;
	betavec = ones(Nc*Nt,1)*sqrt(beta);
	% penalize RF during the rewinder to force it to zero
	betavec = betavec+1000000*kron(ones(Nc,1),[zeros(NN(1),1);ones(NN(2),1)]);
%	R = diag_sp(betavec);
	R = Gdiag(betavec); % jf

	niters = floor(Nt/4);
	printm 'Designing initial pulse'
	[xS, info] = qpwls_pcg(zeros(Nc*Nt,1),Gsml,1,d*tipangle*pi/180,0, ...
			R,1,niters,ones(size(d)));
	B1(:,1) = xS(:,end);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Large-tip pulse design
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isvar('B1_final')
	% build the fast large-tip object
	Nsvd = 6;
	mdom = 1; % this is a magnetization-domain design
	[Glrg, ms0, mzs0] = Gocrf(k, [FOV FOV], [dim dim], mdom, dt, ...
		'd', d, 'indmask', mask, ...
		'sens', reshape(sens,[dim*dim Nc]), ...
		'baseB1', B1(:,1), 'g', g, 'Nsvd', Nsvd, ...
		'we', fmap, 'ti', 0:dt:(Nt-1)*dt);

	nrmse = @(x,y,mask) sqrt(mean(abs(x(mask) - y(mask)).^2) / mean(abs(x(mask).^2)));
	err(1) = 1;
	err(2) = nrmse([md;mzd],[ms0;mzs0],[mask;mask])

	printm 'Running fast large-tip design'
	cgiter = 25;
	db1beta = 0.01;
	ii = 2;
	while ii == 2 | (err(end) <= 0.9999*err(end-1) & ii <= 50)

		% error weighting includes a 'total pulse' regularization (the betavec's)
	%	W = spdiag([ones(2*dim*dim,1);betavec;betavec]);
		W = Gdiag([ones(2*dim*dim,1);betavec;betavec]); % jf

		% pre-subtract the baseline excitation from the desired excitation
		if ii == 2
			mdcg = -ms0 + md;
			mzdcg = -mzs0 + mzd;
		else
			mdcg = -ms + md;
			mzdcg = -mzs + mzd;
		end
	
		printm 'Designing Perturbations'
		[xS,info] = qpwls_pcg(zeros(2*Nc*Nt,1),Glrg,W,[mdcg(:);mzdcg(:);-real(B1(:,end));-imag(B1(:,end))],0,db1beta,1,cgiter);
	
		% get complex pulse back
		db1 = xS(1:Nc*Nt,end) + 1i*xS(Nc*Nt+1:end,end);

		% simulate the total pulse, decreasing the perturbation magnitude
		% if error isn't decreased
		keeptrying = 1;ntrys = 0;
		while keeptrying
			B1tmp = B1(:,end) + db1;

			% simulate the designed pulse
			printm 'Simulating Total Pulse'
			[Glrg, ms, mzs] = feval(Glrg.arg.update, Glrg, B1tmp);

			errtmp = nrmse([md;mzd],[ms;mzs],[mask;mask]);

			if errtmp > 0.9999 * err(end)
				% reduce db1 magnitude, re-simulate pulse
				db1 = db1/2;
				ntrys = ntrys + 1;
				if ntrys > 4 % looks like we have converged, stop trying
					keeptrying = 0;
				end
			else
				% we got a good improvement, proceed
				keeptrying = 0;
				ntrys = 0;
			end
		end

		B1(:,end+1) = B1tmp;
		err(end+1) = errtmp

		ii = ii + 1;
	end

	B1_init = reshape(B1(:,1),[Nt Nc]);
	B1_final = reshape(B1(:,end),[Nt Nc]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show the results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% plot channel 1's pulses
subplot(221)
plot(0:dt*1000:(NN(1)-1)*dt*1000,abs(B1_init(1:NN(1),1)));
hold on
plot(0:dt*1000:(NN(1)-1)*dt*1000,abs(B1_final(1:NN(1),1)),'r');
xlabel 'Time (ms)'
ylabel '|b_1(t)| (a.u.)'
axis([0 (NN(1)-1)*dt*1000 0 max(abs([B1_init(:,1);B1_final(:,1)]))]);
subplot(222)
plot(0:dt*1000:(NN(1)-1)*dt*1000,angle(B1_init(1:NN(1),1)));
hold on
plot(0:dt*1000:(NN(1)-1)*dt*1000,angle(B1_final(1:NN(1),1)),'r');
xlabel 'Time (ms)'
ylabel '\angle b_1(t) (Radians)'
axis([0 (NN(1)-1)*dt*1000 -pi pi]);
legend('Small-tip','Fast large-tip');

% meshplot Mz
subplot(223)
mesh(xx,yy,mzs0)
title 'Small-tip'
axis([-FOV/2 FOV/2 -FOV/2 FOV/2 -1 1]);
xlabel 'x (cm)'
ylabel 'y (cm)'
zlabel 'M_z'
subplot(224)
mesh(xx,yy,mzs)
title 'Fast large-tip'
axis([-FOV/2 FOV/2 -FOV/2 FOV/2 -1 1]);
xlabel 'x (cm)'
ylabel 'y (cm)'
zlabel 'M_z'
