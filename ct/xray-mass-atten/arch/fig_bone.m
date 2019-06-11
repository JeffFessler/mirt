% fig_bone.m
% compare NIST and Sukovic bone data

b1 = load('../bone.dat');
b1(:,1) = b1(:,1) * 1000; % put in keV

bs = load('bone-sukovic.dat');
bs(:,1) = bs(:,1) / 1000; % put in keV
bs(:,2) = bs(:,2) * 10; % why this factor of 10?

clf
subplot(131)
i1 = b1(:,1) >= 10 & b1(:,1) <= 140;
loglog(b1(:,1), b1(:,2), 'y.-', b1(i1,1), b1(i1,2), 'co-')
set(gca, 'xtick', 10 .^ [0:2:5]), xlabel 'E'
subplot(132)
i2 = bs(:,1) >= 10 & bs(:,1) <= 140;
loglog(bs(:,1), bs(:,2), 'g.-', bs(i2,1), bs(i2,2), 'ro-')
set(gca, 'xtick', 10 .^ [0:2:5])
subplot(133)
plot(b1(:,1), b1(:,2), 'c-o', bs(:,1), bs(:,2), 'y-x')
%plot(b1(i1,1), b1(i1,2), 'c-o', bs(i2,1), bs(i2,2), 'y-x')
axisx([10 200]), axisy([0 6]), legend('NIST', 'Sukovic')
