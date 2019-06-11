 function argpair = fatrix2_construct_args(ob, varargin)
%function argpair = fatrix2_construct_args(ob, options)
%|
%| Given a fatrix2_block object, return a cell array of name, value pairs
%| such that the object could be re-constructed by
%|	ob = fatrix2(argpair{:})
%| This is useful for constructing variations of the object.
%|
%| out
%|	argpair {2n}	{name1, value1, ...}
%|
%| Copyright 2010-12-02, Jeff Fessler, University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end

% user-definable properties
name_user = {'arg', 'dim', 'caller', 'odim', 'mask', 'meth'};

% private/internal properties
name_priv = {'dim', 'is_transpose', 'nblock', 'iblock'}

% handles
name_hand = {'back', 'forw', 'power', 'abs', 'free', 'gram', ...
	'mtimes_block', 'block_setup', 'blockify_data'};

argpair = {};

for ii=1:numel(name_user)
	name = name_user{ii};
	argpair{end+[1:2]} = {name, ob.(name)};
end

for ii=1:numel(name_priv)
	name = name_priv{ii};
	argpair{end+[1:2]} = {name, ob.(name)};
end

for ii=1:numel(name_hand)
	name = name_hand{ii};
	argpair{end+[1:2]} = {name, ob.(['handle_' name])};
end
