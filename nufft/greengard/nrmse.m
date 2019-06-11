function nn = nrmse(xtrue, xhat, dummy)

nn = norm(xhat(:) - xtrue(:)) / norm(xtrue(:));
