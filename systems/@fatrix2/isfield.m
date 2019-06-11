function tf = isfield(ob, varargin)
tf = isfield(struct(ob), varargin{:});
