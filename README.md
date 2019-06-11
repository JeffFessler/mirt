Michigan Image Reconstruction Toolbox (MIRT) - Matlab Version
https://github.com/JeffFessler/mirt


This directory contains various algorithms for image reconstruction and
other inverse problems such as image restoration and image registration.
It also contains code for related applications including MRI pulse design.

For many years this toolbox was posted as a (too large) .tgz file here:
http://web.eecs.umich.edu/~fessler/code/index.html

This github version is a work-in-progress.
All the code is here,
but to keep the repo size small,
the data and the compiled mex files are not included.
Currently one must still use the .tgz file to get all of that. 
I will eventually work on some other way to distribute the data and mex files.

This Matlab code is no longer being very actively maintained
because I am switching to using Julia instead.  See:
https://github.com/JeffFessler/MIRT.jl

## GETTING STARTED

After installing the toolbox, use matlab's `path` functionality to put
the top level directory in its path (or launch matlab from that directory).
Then run the file `setup.m` that will add all the appropriate subdirectories
to the path.  You may find it convenient to read `setup.m` and customize it.

I recommend running and examining some of the files in the `example/` directory
or any of the many `..._example.m` files around, such as
	`emission/eml_osem_example.m`

Many example files prompt you to hit enter before continuing, so you (and I)
can see the output of each stage before proceeding.
To change this behavior, execute `prompt run`
(See `utilities/prompt.m` for help.)

A few examples may require the binary programs `wt` or `op` that are part
of Aspire.  You can also get Aspire for free by following the instructions at
http://web.eecs.umich.edu/~fessler/aspire/index.html
There are also mex files that you may need:
`wtfmex` and `f3dmex` for example.
I distribute these only in linux/mac formats;
see `mex/` directory.

Part of my motivation for creating these files is to accompany a book on
image reconstruction that I am currently writing.  If you have any problems
with these m-files, or any suggestions, I welcome your input!


Jeff Fessler, http://web.eecs.umich.edu/~fessler/


## Subdirectories (in alphabetical order):

* `align`
	image registration tools

* `blob`
	SPECT reconstruction with blob basis functions (not recommended)

* `contrib`
	algorithms contributed by others.  these directories are not added to
	the path by setup.m so the user must modify the path to use them.

* `contrib/ppcd`
	test routines comparing WLS-CD, WLS-GCD, WLS-PPCD
	(these are mostly for internal UM use)

* `ct`
	polyenergetic CT routines (beam hardening, dual energy, etc.)

* `data`
	data for examples (not in github version)

* `doc`
	see the pdf file within (not in github version)
	for some introductory documentation.

* `emission`
	algorithms for Poisson emission tomography PET/SPECT/ Poisson regression

	+ `eml_` emission maximum likelihood
	+ `eql_` emission quadratically penalized likelihood
	+ `epl_` emission penalized likelihood

* `example`
	example(s) of usage.  there are more examples in other directories too.
	running any of these examples is a good place to start!

* `fbp`
	filter-backproject reconstruction, including 2D parallel and fan-beam
	and 3D Feldkamp (FDK) cone beam reconstruction

* `freemat` (not in github version)
	work in progress, towards making the code run with freemat
	(obsolete: use octave instead of freemat)
	(better yet, just use Julia)

* `general`
	some algorithms that work for generic image reconstruction problems

* `graph`
	graphics functions

* `mex`
	MEX (matlab executables), including some C99 source code

* `mri`
	MR image reconstruction

* `mri-rf`
	MR pulse design tools, including Spectral-spatial pulse design
	for phase precompensatory slice selection.  (more to come)

* `nufft`
	non-uniform FFT (NUFFT) toolbox

* `octave` (not in github version)
	work in progress, towards making the code run with octave
	(this is a moving target as both octave and matlab evolve)

* `penalty`
	functions associated with regularization

* `systems`
	system matrices and system matrix object classes
	If you are interested in edge-preserving image restoration
	for a shift-invariant blur model with additive gaussian noise,
	then start with `systems/Gblur_test.m` and `example/restore_example.m`
	For 2D tomography, consider starting with `systems/Gomo2_strip.m`,
	which is used in many of the examples.

* `transmission`
	algorithms for Poisson transmission tomography

	+ `tml_` transmission ML
	+ `tql_` transmission quadratically penalized likelihood
	+ `tpl_` transmission penalized likelihood

* `utilities`
	useful functions for image reconstruction algorithms.

* `wls`
	algorithms associated with the weighted least squares (WLS)
	cost function and penalized versions thereof

	+ `pwls_`	penalized weighted least squares
	+ `qpwls_`	quadraticaly penalized weighted least squares

Most algorithms also include a test routine...


## Additional notes:

Raymod Muzic has matlab routines for reading ECAT files available:
	http://www.nuclear.uhrad.com/comkat
(I have not yet tried them myself.)

One of many annoying issues with Matlab is that it can store sparse matrices
only as doubles, wasting memory, and if you do `S * x`, where `S` is a
sparse matrix and `x` is a vector of class single, Matlab (as of 2016a)
gives an error message rather than politely upgrading `x` to a double.
The object `Gsparse.m` provides a work around for this.
Complain to Mathworks that they should fix this annoyance...

Windows users:

Some of the subdirectories of the `systems` directory contain "links"
to m-files in other directories.  (These are created using `ln -s` in
unix.)  These links are also in the `.tar` file as soft links.  But these
links may not be recognized by Windoze, resulting in various error messages.
They work fine on Mac OSX since it is unix "under the hood."
I recommend that you avoid using Windows.
But if you insist, then you will have to figure out how to fix those links
or copy the appropriate m-files into the appropriate directories.

Another problem is that apparently windoze is case insensitive.
I believe I have purged most of the m-files that had capitalized names now.
Nevertheless, at this point it would be better to just install linux instead.
Or switch to Julia.
