 function subplo(varargin)
% my version of subplot

persistent dim
persistent top		% top percentage
persistent bottom	% bottom percentage
persistent bottom0	% bottom percentage
persistent hshrink	% shrink horizontal (to squeeze plots)

if length(varargin) < 1
	printf('top=%g bottom=%g bottom0=%g dim=', top, bottom, bottom0)
	disp(dim)
	return
end

if ischar(varargin{1})
	while length(varargin)
		arg1 = varargin{1};
		if streq(arg1, 'new')
			if length(varargin) < 2, error 'need dim', end
			dim = varargin{2};
			varargin = {varargin{3:end}};

		elseif streq(arg1, 'bottom')
			if length(varargin) < 2, error 'need bottom', end
			bottom = varargin{2} / 100;
			varargin = {varargin{3:end}};

		elseif streq(arg1, 'bottom0')
			if length(varargin) < 2, error 'need bottom0', end
			bottom0 = varargin{2} / 100;
			varargin = {varargin{3:end}};

		elseif streq(arg1, 'top')
			if length(varargin) < 2, error 'need top', end
			top = varargin{2} / 100;
			varargin = {varargin{3:end}};

		elseif streq(arg1, 'hshrink')
			if length(varargin) < 2, error 'need hshrink', end
			hshrink = varargin{2};
			varargin = {varargin{3:end}};

		else
			error 'unknown arg'
		end
	end
	return
end

if isempty(dim)
	error 'subplot set first'
end

nr = dim(1);
nc = dim(2);

if isempty(top)
	top = 0;
end
if isempty(bottom)
	bottom = 0;
end
if isempty(bottom0)
	bottom0 = 0;
end
if isempty(hshrink)
	hshrink = 1;
end

if length(varargin) == 1
	arg = varargin{1};
	if arg < 1 || arg > nr * nc, error 'bad plot index', end
	row = 1 + floor((arg-1) / nc);
	col = 1 + mod(arg-1, nc);
elseif length(varargin) == 2
	row = varargin{1};
	col = varargin{2};
else
	error 'too many args'
end

s.width = 1 / nc * hshrink;
s.height = (1 - bottom0) / (nr * (1+top+bottom));
s.left = (col-1) / nc * hshrink;
s.bottom = (1 - row / nr)*(1-bottom0) + bottom/nr + bottom0;
subplot('position', [s.left s.bottom s.width s.height])
