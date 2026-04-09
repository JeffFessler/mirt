mex/readme

This directory contains compiled MEX files for some subroutines.
Also provided are some of the source code files for compiling MEX files.

Windows users: mathworks provides help on their web pages for compiling mex files.


2006-6-20

getting mex files to run on linux has proved to be somewhat challenging.
at issue seems to be to use the "correct" (matlab endorsed?) gcc version.
see:
	http://www.mathworks.com/support/tech-notes/1600/1601.html

today, i switched to gcc 3.4.4 for linux compiles.  locally, i could
not get the mex files to run unless i following this tech note:
	http://www.mathworks.com/support/solutions/data/1-2H64MF.html?1-2H64MF
which required this addition to my .cshrc:
	setenv LD_PRELOAD /lib/libgcc_s.so.1


2005-3-23 (circa when matlab7 first came out)

unfortunately, it seems that v6 and v7 matlab mex files
are not quite compatible.  in particular, compiling mex files under v7
leads to mex files that v6 will not execute.
however, mex files compiled under v6 *will* run under v7, but it seems
that if those mex files call mexErrMsgTxt(), then v7 will abort,
due to exception handling issues.

mathworks claims the opposite here:
http://www.mathworks.com/access/helpdesk/help/techdoc/rn/externab.html

well, i have purged almost all of my mexErrMsgTxt() calls from my mex
files, so now the mex files should work under v7, although as of
2005-3-23 i have not tested them extensively under v7.   but i will
be moving in that direction soon...
