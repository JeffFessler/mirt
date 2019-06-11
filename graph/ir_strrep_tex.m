 function out = ir_strrep_tex(in)
%function out = ir_strrep_tex(in)
%|
%| tricks for image processing book that allow me to use macros in commands:
%| titlef, xlabelf, ylabelf, ir_text
%| replace \kx with code for k_X
%| replace \yx with code for k_Y
%|
%| 2016-05-10 Jeff Fessler

if nargin < 1, help(mfilename), fail(mfilename), end

kx = '\ensuremath{\nu_{\mathrm{\scriptscriptstyle X}}}';
ky = '\ensuremath{\nu_{\mathrm{\scriptscriptstyle Y}}}';
omf = '\ensuremath{\Omega}';
omx = '\ensuremath{\Omega_1}';
omy = '\ensuremath{\Omega_2}';
omz = '\ensuremath{\Omega_3}';

bfun = @(c) ['\ensuremath{\mathbf{' c '}}'];

pairs = {
'\x', bfun('x');
'\y', bfun('y');
'\A', bfun('A');
'\U', bfun('U');
%
'\Ex', '\ensuremath{\mathrm{E}}'; % \mathsf fails
'\usym', '\ensuremath{\nu}';
'\Dx', '\ensuremath{\Delta_{\mathrm{\scriptscriptstyle X}}}';
'\Dy', '\ensuremath{\Delta_{\mathrm{\scriptscriptstyle Y}}}';
'\gxy', '\ensuremath{g(x,y)}';
%
'\bmn', '\ensuremath{b[m,n]}';
'\cmn', '\ensuremath{c[m,n]}';
'\fmn', '\ensuremath{f[m,n]}';
'\fhmn', '\ensuremath{\hat{f}[m,n]}';
'\gmn_1', '\ensuremath{g_1[m,n]}';
'\gmn_2', '\ensuremath{g_2[m,n]}';
'\gdmn', '\ensuremath{g_d[m,n]}';
'\gmn', '\ensuremath{g[m,n]}';
'\hmn', '\ensuremath{h[m,n]}';
'\wmn', '\ensuremath{w[m,n]}';
'\smn', '\ensuremath{s[m,n]}';
'\xmn', '\ensuremath{x[m,n]}';
'\vmn', '\ensuremath{v[m,n]}';
'\ymn', '\ensuremath{y[m,n]}';
'\zmn', '\ensuremath{z[m,n]}';
'\Rfmn', '\ensuremath{R_f[m,n]}';
'\Rvmn', '\ensuremath{R_v[m,n]}';
'\Rxmn', '\ensuremath{R_x[m,n]}';
'\Rymn', '\ensuremath{R_y[m,n]}';
'\Rzmn', '\ensuremath{R_z[m,n]}';
'\hRsmn', '\ensuremath{\hat{R}_s[m,n]}';
'\hRxmn', '\ensuremath{\hat{R}_x[m,n]}';
'\hRymn', '\ensuremath{\hat{R}_y[m,n]}';
'\hRzmn', '\ensuremath{\hat{R}_z[m,n]}';
'\Om', '\ensuremath{\Omega}';
'\kx', kx;
'\ky', ky;
'\omf', omf;
'\omx', omx;
'\omy', omy;
'\omz', omz;
'\rect', '\mathrm{rect}';
'\jinc', '\mathrm{jinc}';
'\sinc', '\mathrm{sinc}';
'\tri', '\mathrm{tri}';
'\cov', '\mathrm{cov}';
'\mod', '\mathrm{mod}';
'\var', '\mathrm{var}';
'\cconv', '\; {**} \;';
'\cos', '\mathrm{cos}';
'\sin', '\mathrm{sin}';
};

out = in;
for ip=1:size(pairs,1)
	expr = pairs{ip,1};
	expr = strrep(expr, '\', '\\');
	expr = [expr '(?![a-zA-Z])']; % not followed by a letter
	repl = pairs{ip,2};
	repl = strrep(repl, '\', '\\');
%	out = strrep(out, pairs{ip,1}, pairs{ip,2});
	tmp = regexprep(out, expr, repl);
	if 0 && ~streq(tmp, out)
		pr expr
		pr out
		pr tmp
	end
	out = tmp;
end

amps = strfind(out, '&');
switch numel(amps)
case 0
	% ok;
case 1
	amp = amps(1);
	if amp == 1
		out = ['\' out];
	elseif out(amp-1) ~= '\'
		out = [out(1:amp-1) '\' out(amp:end)]; % insert ampersand
	end
otherwise
	fail 'not done'
end

%{
for ii=1:numel(amps)
	amp = amps(ii);
end
%}
