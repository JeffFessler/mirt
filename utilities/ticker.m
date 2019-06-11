 function ticker(varargin)
%|function ticker(varargin)
%|
%| This utility is used in functions that might take a long time
%| to provide a status update so the user knows it is working.
%| "ticker print" to make output print to command window (default)
%| "ticker waitbar" to display using waitbar
%| ticker(waittime) to set time between updates to "waittime"
%| the default is 1 seconds
%| ticker(mfilename, i, n) to report at ith of n iterations
%| trick: i can be 2 elements (iteration, subset) and nn can be (#iter, #subset)
%|
%| Copyright 2003-11-29, Jeff Fessler, University of Michigan

persistent Ticker % stores state

if isempty(Ticker)
	Ticker = ticker_default;
end

% query state
if ~nargin
	printf('Ticker.how=%s Ticker.wait=%g Ticker.who="%s"', ...
		Ticker.how, Ticker.wait, Ticker.who)
return
end


% set mode from command line
if nargin == 1
	arg = varargin{1};
	if ischar(arg)
		if streq(arg, 'help')
			help(mfilename)
		elseif streq(arg, 'form')
			Ticker.form = arg;
		elseif streq(arg, 'reset')
			Ticker = ticker_default;
		elseif streq(arg, 'test')
			ticker_test
		else
			Ticker.how = arg;
		end
	else
		Ticker.wait = arg; % time
	end
return
end


% set who and mode (or give help)
if nargin == 2
	if ~ischar(varargin{1})
		fail 'first arg must be mfilename'
	end

	if isempty(Ticker.who)
		Ticker.who = varargin{1};
	else
		return % some other mfile has control!
	end

	if ischar(varargin{2})
		Ticker.how = varargin{2};
	else
		Ticker.wait = varargin{2};
	end
return
end


% within a loop
if nargin == 3
	ww = varargin{1};
	ii = varargin{2};
	nn = varargin{3};

	if isempty(Ticker.who)
		Ticker.who = ww;
	end

	if ~streq(ww, Ticker.who, length(Ticker.who))
		% fix: what about nested uses of ticker?
%		printf(['ticker who mismatch? ' Ticker.who ' ' ww])
		return
	end

	if streq(Ticker.how, 'waitbar')
		if isempty(Ticker.h)
			Ticker.h = waitbar(ii(1)/nn(1));
		else
			waitbar(ii(1)/nn(1), Ticker.h)
		end
		if ii(1)==nn(1)
			close(Ticker.h)
			Ticker = ticker_default;
		end

	elseif streq(Ticker.how, 'print')
		if isempty(Ticker.t0)
			Ticker.t0 = clock;
		end
		if ii(1) > 0
			t1 = clock;
			te = etime(t1, Ticker.t0);
			if te >= Ticker.wait
				Ticker.t0 = t1;
				switch Ticker.form
				case 'hour'
%					form = sprintf('%d:%2d:%2d', t1(4:6));
%					form = datestr(t1-Ticker.t0, 'HH:MM:SS');
					te = etime(t1, Ticker.tinit);
					form = sprintf('%.1f', te/60);
				otherwise
					form = '';
				end
				if length(nn) == 1
					printf('%s: %d of %d %s', ww, ii, nn, form)
				else
					printf('%s: %d,%d of %d,%d %s', ...
						ww, ii(1), ii(2), nn(1), nn(2), form)
				end
			end
		end
		if ii(1) == nn(1) % reset
			Ticker = ticker_default;
		end

	else
		error('unknown Ticker mode %s', Ticker.how)
	end
end


function t = ticker_default
t.who = '';		% which mfile controls
t.how = 'print';	% how to display
t.h = [];		% waitbar handle
t.t0 = [];		% start of recent elapsed time
t.tinit = clock;	% start of total elapsed time
t.wait = 1;		% (seconds) how often to update
t.form = 'hour';	% format for displaying time


% self test
function t = ticker_test
ticker reset
nloop=100;
for ii=1:nloop
	ticker(mfilename, ii, nloop)
	pause(0.1)
end
