%function h = jim(varargin) - alias to "im" for julia compatibility
 function h = jim(varargin)
   out = im(varargin{:});
if nargout
   h = out;
end
end
