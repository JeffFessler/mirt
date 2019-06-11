function p = jf_normcdf(pn, x, mu, sig)
% avoid matlab's stat toolbox...
z = (x - mu) / sig; 
p = 0.5 * erfc(-z ./ sqrt(2));
