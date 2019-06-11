 function [kbli, alpha, kb_m] = kaiser_bessel_li(x, J, alpha, kb_m)
%function [kbli, alpha, kb_m] = kaiser_bessel_li(x, J, alpha, kb_m)
%
%	Line integrals of generalized Kaiser-Bessel function 
%	at radial positions x in support [-J/2,J/2]
%	shape parameter "alpha" (default 2.34 J)
%	order parameter "kb_m" (default 0)
%	see (A1) in lewitt:90:mdi, JOSA-A, Oct. 1990
%	in
%		x	[M,1]	arguments
%	out
%		kbli	[M,1]	KB line integrals

if ~isvar('J'), J = 4; end
if ~isvar('alpha') | isempty('alpha'), alpha = 2.34 * J; end
if ~isvar('kb_m') | isempty('kb_m'), kb_m = 1; end

if alpha == 0		% Uniform circle
	kbli = zeros(size(x));
	ii = abs(x) < J/2;
	kbli(ii) = 2 * sqrt((J/2)^2 - x(ii).^2);
else
	kbli = kaiser_bessel(x, J, alpha, kb_m+0.5);
	kbli = J/2 * sqrt(2*pi/alpha) * besseli(kb_m+0.5, alpha) / besseli(kb_m, alpha) * kbli;
end

