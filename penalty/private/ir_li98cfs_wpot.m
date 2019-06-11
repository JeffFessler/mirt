function w = ir_li98cfs_wpot(z, d)
% weighting for potential function from li:98:cfs
z = z ./ d;
w = atan(z);
w(z == 0) = 1;
z(z == 0) = 1;
w = w ./ z;
