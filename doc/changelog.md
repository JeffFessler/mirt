This file is an extremely incomplete list of changes,
in reverse chronological order.

2016-10-24
added odiag and idiag features to fatrix2.
added Gblock(., ., 'odiag', d) feature
removed chat option Gblock(.,.,chat) 

2012-09-23
modified strum so that methods must have one and only output,
because octave (and freemat) do not do "nargout" well.
even matlab does not do nargout with implicit functions properly.

2007-03-05
moved my assert.m to jf_assert.m because matlab 7.4 has its own assert.m
