%
% convergence_test.m
% Compare convergence rate between WLS-GCD, WLS-PPCD, and WLS-CD algorithms
%

clear all,close all
load testdata
tol = 4.4522974246e4;         
[ma na] = size(x0);

niter = 50;
% 4x4 groups of WLS-GCD algorithm
disp('running gcd')
[xg,phig,dg,timeg,niterg] = wls_gcd(A, W, yy, x0, [4 4], niter, tol);

% 4 sets of WLS-PPCD algorithm
disp('running ppcd')
[xp,phip,dp,timep,niterp] = wls_ppcd(A, W, yy, x0, 2, 2, niter, tol);

% WLS-CD algorithm
disp('running cd')
[xc,phic,dc,timec,niterc] = wls_cd(A, W, yy, x0, niter, tol);

% Calculate normalized norms of estimates at each iterations
for i=1:niterg
  normxg(i) =norm(xg(:,i)-xhat)/norm(xhat);
end;
for i=1:niterp
  normxp(i) = norm(xp(:,i)-xhat)/norm(xhat);
end;
for i=1:niterc
   normxc(i) = norm(xc(:,i)-xhat)/norm(xhat);
end;

figure(1),orient tall
subplot(211),plot(1:niterg,phig,'go-.'),hold on
plot(1:niterp,phip,'m*:'),hold on
plot(1:niterc,phic,'bd--'),hold off,
,xlabel('#iterations'),ylabel('\Phi decrease'),
title('\Phi Decrease vs. Number of Iterations'),
legend('4x4 GCD','4-PPCD','CD')

subplot(212),semilogy(1:niterg,normxg,'g-.'),hold on,
semilogy(1:niterp,normxp,'m-'),hold on,
semilogy(1:niterc,normxc,'b--'),hold off,
xlabel('#iterations'),ylabel('Normalized norm'),
title('Normalized Norm vs. Number of Iterations'),
legend('4x4 GCD','4-PPCD','CD')

