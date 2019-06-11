#!/bin/csh
# do.sh osemn {post|iter} [mean|resol|bias|noise|std]
# do.sh ospsn {beta|post|iter} [mean|resol|bias|noise|std]

set arg = ""
if ($#argv >= 1) then
        set arg = $1
endif

alias printf /usr/bin/printf

alias i $HOME/blob/c/i
alias op $HOME/blob/c/op

# Choose image display type
#alias ji 'echo DISPLAY \!*; j --red -b -1 -nrow 5 \!*'	# Show image
#alias js 'echo DISPLAY \!*; j --red -b -1 -nrow 10 \!*' # Show projection views
alias ji 'op range'					# Just print range
alias js 'op range'					# Just print range


if 0 then			# Which phantom to do?
	set case = sph		# Small sphere phantom
else if 0 then
	set case = unif		# Uniform cylinder phantom
else
	set case = chest	# Chest phantom
	set backg = 2		# True background (soft tissue) value
endif

if 0 then			# To downsample or not?
	set down = no
else
	set down = yes
endif

if ($case == sph) then
	set nx = 64		# Matrix size
	set ny = $nx
	set nz = 25		# Number of slices

	set rad = `echo "2 * $nx/64" | bc -l`	# Sphere radius
	set act = 2		# Sphere value
else if ($case == unif) then
	set nx = 64		# Matrix size
	set ny = $nx
	set nz = 25		# Number of slices

	set rad = 28		# Cylinder radius
	set act = 2		# Cylinder value
else if ($case == chest) then
	set nx = 64		# Matrix size
	set ny = $nx
	set nz = 47		# Number of slices
	set nz = 23
endif

set nu = $nx		# Projection size
set nv = $nz
set na = 60		# Number of projections

set orbit = 360		# Orbit range
set ostart = 0		# Orbit start

if 0 then
	set bl = uhe
	set sx = `echo "64 /$nx * 7.196" | bc -l`	# Pixel size for I-131
else
	set bl = tc
	set sx = `echo "64 /$nx * 4" | bc -l`	# Pixel size for Tc-99
endif
set sx = `printf '%g' $sx`
set su = $sx					# Ray spacing
set sv = $sx					# Ray spacing

set sfilter = 1		# Interpolation filter used for rotations

# Set directory for phantom image
#set phdir	= /Users/nastazia
set phdir	= /n/ir7/y/fessler/data/phantom

# Set directory for results
#set dir = /Users/nastazia/blob/data
set dir = /n/ir22/y/nastazia/blobdata

set cdir = $dir/$case
set indir = $cdir/phantom
set fbpdir = $cdir/fbp,hanning,0,1

if 0 then		# Rotation- or blob-based system model?
	set btype = conv1sym
else
	#set R		= `echo "1.25 * $sx" | bc -l`	# Blobnn
	#set alpha	= 3.6
	#set R		= `echo "1.50 * $sx" | bc -l`	# Blobn
	#set alpha	= 6.4
	set R		= `echo "2.00 * $sx" | bc -l`	# Blob
	set alpha	= 10.4
	#set R		= `echo "2.25 * $sx" | bc -l`	# Blobw
	#set alpha	= 3.5
	#set R		= `echo "2.00 * $sx" | bc -l`	# Blob0w
	#set alpha	= 0
	#set R		= `echo "2.60 * $sx" | bc -l`	# 
	#set alpha	= 14.5
	set R		= `printf '%g' $R`
#set alpha	= `echo "2.34 * 2 * $R" | bc -l`
#set alpha	= `printf '%1.4f' $alpha`
	set kbm	= 2
#set kbm		= 1
	set dR		= `echo "$sx / 10" | bc -l`
#set dR		= `echo "$sx / 50" | bc -l`
	set dR		= `printf '%1.4f' $dR`
	set dR		= `printf '%g' $dR`
	#set dR		= $sx
	set btype	= blob,$dR,$su,$sv
endif
if ($btype == conv1sym) then
	set outdir = $cdir/conv1sym,b$sfilter
else
	set outdir = $cdir/blob,R$R,dR$dR,a$alpha,m$kbm
endif

# Create directories
if !(-d $dir)		mkdir $dir
if !(-d $cdir)		mkdir $cdir
if !(-d $indir)		mkdir $indir
if !(-d $fbpdir)	mkdir $fbpdir
if !(-d $outdir)	mkdir $outdir

# Dimensions of original (high-res) projection set
if ($down == yes) then
	set m = 2
	set m = 4
else
	set m = 1
endif
set nxP	=	`echo "$m * $nx" | bc`
set nyP	=	`echo "$m * $ny" | bc`
set nzP	=	`echo "$m * $nz" | bc`
set sxP =	`echo "$sx / $m" | bc -l`
set sxP =	`printf '%g' $sxP`
set suP =	$sxP
set svP =	$sxP

if (($case == sph) || ($case == unif)) then
	set radP =	`echo "$m * $rad" | bc`
endif

# Set data filenames
set tmp = /tmp						# Directory for temp. files

set phantom	= $indir/xtrue,$nx,$ny,$nz.fld		# Phantom
set phantomP	= $indir/xtrue,$nxP,$nyP,$nzP.fld
set phantomP_imp = $indir/xtrue,imp,$nxP,$nyP,$nzP.fld	# Phantom plus an impulse
set mask	= $indir/mask,$nx,$ny,$nz.fld		# System mask
set maskP	= $indir/mask,$nxP,$nyP,$nzP.fld
set mumap	= $indir/mumap,$nx,$ny,$nz.fld		# Attenuation map
set mumapP	= $indir/mumap,$nxP,$nyP,$nzP.fld
set proj	= $indir/proj,$nx,$nz,$na.fld		# Sinogram
set proj_imp	= $indir/proj,imp,$nx,$nz,$na.fld
set back	= $outdir/back.fld			# Plain back-projection
set chang	= $fbpdir/chang.fld			# Chang atten. correction factors

if ($btype == conv1sym) then
	set psfs = $dir/psfs230,$bl,$nx.fld		# Depth-dependent blurring PSFs
	set psfsP = $dir/psfs230,$bl,$nxP.fld
else
	set blobli	= $dir/blobli,$bl,R$R,dR$dR,alpha$alpha,m$kbm.fld
	set radii	= $dir/radii,$bl,R$R,dR$dR,alpha$alpha,m$kbm.fld
	set psfs	= $blobli\:$radii

	if (($arg == blobli) || !(-e "$blobli") || !(-e "$radii")) then
		echo "cd $dir; blob_li_blur($R, $dR, $alpha, $kbm, '$bl', $su, $sx, $ny, 1);" \
			| matlab -nodesktop -nosplash

		op range $blobli
		op range $radii
		exit
	endif

	set RP		= `echo "2.00 * $sxP" | bc -l`	# Blob
	set RP =	`printf '%g' $RP`
	set dRP		= `echo "$sxP / 10" | bc -l`
	set dRP =	`printf '%g' $dRP`

	set blobliP	= $dir/blobli,$bl,R$RP,dR$dRP,alpha$alpha,m$kbm.fld
	set radiiP	= $dir/radii,$bl,R$RP,dR$dRP,alpha$alpha,m$kbm.fld
	set psfsP	= $blobliP\:$radiiP

	if (!(-e $proj) && (!(-e "$blobliP") || !(-e "$radiiP"))) then
		echo "cd $dir; blob_li_blur($RP, $dRP, $alpha, $kbm, '$bl', $suP, $sxP, $nyP, 1);" \
			| matlab -nodesktop -nosplash

		op range $blobliP
		op range $radiiP
		exit
	endif
endif

# Set parameters
set ftype	= 3s@$sx,$sx,$sx,$orbit,$ostart,$sfilter,$btype@$mumap@$psfs@-$nx,$nz,$na
			# 3D SPECT model
set ftypeP	= 3s@$sxP,$sxP,$sxP,$orbit,$ostart,$sfilter,$btype@$mumapP@$psfsP@-$nxP,$nzP,$na
			# 3D SPECT model for generating original (high-res) projection set
set ftype0	= 3s@$sx,$sx,$sx,$orbit,$ostart,$sfilter,none@$mumap@-@-$nx,$nz,$na 
			# 3D SPECT model for Chang calculation

# Create phantom
if (($arg == phantom) || !(-e "$phantom")) then
	# High-res version to generate projection set
	if ($case == sph) then		# Uniform sphere
		op ellipsoid $phantomP $nxP $nyP $nzP 0 0 0 $radP $radP $radP 0 0 $act 3
	else if ($case == unif) then	# Uniform cylinder
		op ellipse t0 $nxP $nyP 0 0 $radP $radP 0 $act 3
		op rep $phantomP t0 $nzP
	else if ($case == chest) then
		set nx0 = 128			# Original chest phantom size
		if ($nxP <= $nx0) then
			set d = `echo "$nx0/$nxP" | bc`

			op slice t0 $phdir/chest/t0 - 0 `echo "$nzP*$d-1" | bc` 0 0 0 0
			op sample mean t0 t0 $d $d	# Down. along x, y
			op transpose t0 t0 1,2		# Down. along z
			op sample mean t0 t0 1 $d
			op transpose $phantomP t0 1,2
		else
			set d = `echo "$nxP/$nx0" | bc`

			op slice t0 $phdir/chest/t0 - 0 `echo "$nzP/$d-1" | bc` 0 0 0 0
			echo "in = fld_read('t0'); out = upsamp(in, $d); mat_write('t1', out, '-nocheck');" \
				| matlab -nodesktop -nosplash
			op float t1 t1
			op blur $phantomP t1 2,7,zer 2,7,zer 2,7,renorm
		endif
	endif

	ji $phantomP

	if ($down == yes) then
	# Low-res version to compare to reconstructions
		op sample mean t0 $phantomP $m $m	# Down. along x, y
		op transpose t0 t0 1,2			# Down. along z
		op sample mean t0 t0 1 $m
		op transpose $phantom t0 1,2

		ji $phantom
	endif

	exit
endif

# Specify system mask
if ($arg == mask || (("$mask" != "-") && !(-e "$mask"))) then
	if ($case == sph) then	
		set mask_rad = `echo "$rad + 2" | bc -l`	# Small mask
	else
		set mask_rad = `echo "$nx / 2 - 2" | bc -l`
	endif
	op ellipse t0 $nx $ny 0 0 $mask_rad $mask_rad 0 1 1	# A circle
	op conv $mask t0 byte

	# Display mask with true image superimposed
	op add t0 $phantom $mask
	ji t0

	if ($down == yes) then
		if ($case == sph) then	
			set mask_rad = `echo "$radP + 2" | bc -l`	# Small mask
		else
			set mask_rad = `echo "$nxP / 2 - 2" | bc -l`
		endif
		op ellipse t0 $nxP $nyP 0 0 $mask_rad $mask_rad 0 1 1	# A circle
		op conv $maskP t0 byte

		# Display mask with true image superimposed
		op add t0 $phantomP $maskP
		ji t0
		endif
	exit
endif

# Specify attenuation map
if ($arg == mumap || (("$mumap" != "-") && !(-e "$mumap"))) then
	if ($bl == uhe) then
		set mu = 0.011				# 0.011/mm for I-131 360KeV in water
	else if ($bl == tc) then
		set mu = 0.015				# 0.015/mm for Tc-99 in water
	endif
	if (($case == sph) || ($case == unif)) then
		set mu_rad_pix = `echo "$nxP / 2 - 3" | bc -l`	# Radius of attenuation disc
		# High-res version to generate projection set
		if 1 then				# A smooth circle
			op ellipse $mumapP $nxP $nyP 0 0 $mu_rad_pix $mu_rad_pix 0 $mu 3
		else					# Non-uniform
			op ellipse t0 $nxP $nyP \
				-8 -8 18 18 0 1 3 \
				19 0 9 9 0 0.5 3 \
				6 18 8 8 0 1 3
			op mul $mumapP t0 - $mu
		endif

		ji $mumapP

		if ($down == yes) then
		# Low-res version to use in reconstruction
			op sample mean $mumap $mumapP $m $m	# Down. along x, y

			ji $mumap
		endif
	else if ($case == chest) then		# Chest phantom atten. map
		# Scale mu values
		op mul t0 $phdir/chest/t1 - `echo $mu / 0.011 | bc -l`

		# High-res version to generate projection set
		set nx0 = 128			# Original chest phantom size
		if ($nxP <= $nx0) then
			set d = `echo "$nx0/$nxP" | bc`

			op slice t0 t0 - 0 `echo "$nzP*$d-1" | bc` 0 0 0 0
			if ($bl == tc) then
				op mul t0 t0 - `echo 0.015/0.011|bc -l`
			endif
			op sample mean t0 t0 $d $d	# Down. along x, y
			op transpose t0 t0 1,2		# Down. along z
			op sample mean t0 t0 1 $d
			op transpose $mumapP t0 1,2
		else
			set d = `echo "$nxP/$nx0" | bc`

			op slice t0 t0 - 0 `echo "$nzP/$d-1" | bc` 0 0 0 0
			echo "in = fld_read('t0'); out = upsamp(in, $d); mat_write('t1', out, '-nocheck');" \
				| matlab -nodesktop -nosplash
			op float t1 t1
			op blur $mumapP t1 2,7,zer 2,7,zer 2,7,renorm
		endif

		ji $mumapP

		if ($down == yes) then
		# Low-res version to use in reconstruction
			op sample mean t0 $mumapP $m $m	# Down. along x, y
			op transpose t0 t0 1,2		# Down. along z
			op sample mean t0 t0 1 $m
			op transpose $mumap t0 1,2

			ji $mumap
		endif
	endif

	exit
endif

# Display depth-dependent blurring PSFs (generated by gen_psf.m)
if ($btype == conv1sym) then
if ($arg == psfs || (("$psfs" != "-") && !(-e "$psfs"))) then
	echo RUN psfs_$bl.m to produce PSF file

	ji $psfs
	op range $psfs
	ji $psfsP
	op range $psfsP
	exit
endif
endif

# Generate projections using 3D SPECT system model
if ($arg == proj || !(-e $proj)) then
	i proj3 t0 $phantomP $ftypeP $maskP
	op nonlin max t0 t0 0 0		# Set negatives to zero 

	if ($down == yes) then		# Downsample projection views
		js t0
		op div t0 t0 - $m	# !!!
		op sample mean $proj t0 $m $m
	else
		cp t0 $proj
	endif

	js $proj
	exit
endif

# Plain back-projection
if ($arg == back || !(-e $back)) then
	i back3 $back $nx $ny $nz $proj - $ftype $mask

	js $back
	exit
endif

# Calculate Chang attenuation correction factors
if ($arg == chang || !(-e $chang)) then
	# To get Chang correction factors, backproject a uniform sinogram
	i back3 t0 $nx $ny $nz - - $ftype0 $mask
	op div0 t1 - t0 $na

	op index t11 $mask - 1 0		# Make it 1 instead of 0 outside mask
	op add $chang t1 t11

	op range $chang
	ji $chang
	exit
endif

# FBP reconstruction
set fbpwin = hanning,0.6,0,1
set fbpout = $fbpdir/fbp,$fbpwin.fld

if ($arg == fbp || !(-e $fbpout)) then
	op transpose t0 $proj 1,2		# Make sinograms nx x na x nz

	set cut = 0.6
	set ic = 0
	set nc = 1

	while ($ic < $nc)
		set fbpwin = hanning,$cut,0,1
		#set fbpwin = boxcar,1,1
		set fbpout = $fbpdir/fbp,$fbpwin.fld

		op fbp t1 t0 user,t $nx $ny $fbpwin 1 1 0  0 0 $orbit $ostart	# Pixel size 1
								# for consistency w/ 3s projector
		op mul $fbpout t1 $chang

		op range $phantom t1 $fbpout
		ji $fbpout

		set cut = `echo "$cut + 0.05" | bc`
		set cut = `printf '%g' $cut`
		@ ic += 1
	end
	exit
endif

# Iterative reconstruction
if 1 then							# Choose initial image
        set init = $fbpdir/fbp,hanning,0.6,0,1.fld		# FBP-reconstructed image

	if (($arg == em) || ($arg == osem)) then
		op nonlin max t9 $init 0 0.001			# Only positive values for EM
		set init = t9
	endif
else
        set init = -$nx,$ny,$nz,1				# Uniform image
endif

# Scale noiseless projections to have the desired mean number of counts.
set count = 5e6						# Average counts (5e6 ~= 20 * 64*64*60)
set ybi = $indir/ybi,$count.fld				# Scatter-free measurements means
set ci = $indir/ci,$count.fld				# Calibration factors

if (($arg == ci) || !(-e $ci)) then
	op sim calib $ybi $ci $proj - 0.0 $count 1	# (seed)

	op range $ci $ybi
	js $ci
	js $ybi
	exit
endif

# Compute mean scattered counts
set scatpc = 10					# Scatter percentage (%)
set t = "`echo $count | sed 's/e/ * 10^/'`"
set scatmean = `echo "$t * $scatpc / 100 / $nu / $nv / $na" | bc -l`
echo "count = $count, scatpc = $scatpc, scatmean = $scatmean"
set obj = 0

#set saver = disp,100
#set saver = save,1
set saver = -
set nsave = 50
set nsave = 100
set saver = save,$nsave

#  Impulse response (OSPS or OSEM or EM) ##############################
if $arg == imp then
	set niter = 25
		set niter = 1500

	if 0 then		# OSPS
		# Choose beta''s by trial&error to make impulse response spherical
		set l2bx = -8
		set l2bz = -3
		set l2b = $l2bx,$l2bz

		set nsubset = 10
		set alg = ospsc,pc,$nsubset,$na,1,0.1
		set penal = 3d,$l2b,quad,5,-

		set out = $outdir/os,$nsubset,$l2b,%04d.fld
		set out_imp = $outdir/os,imp,$nsubset,$l2b,%04d.fld
		set imp = $outdir/imp,$nsubset,$l2b,%04d.fld
		set piximp = $outdir/imp,$nsubset,$l2b,%04d,pix.fld
		set fwhm = $outdir/fwhm,imp,$nsubset,$l2b.txt
	else if 1 then		# OSEM
		set nsubset = 6
		set alg = osemc,fast,$nsubset,$na,1.0
		set penal = -

		set out = $outdir/osem,%04d.fld
		set out_imp = $outdir/osem,imp,%04d.fld
		set imp = $outdir/imp,osem,%04d.fld
		set piximp = $outdir/imp,osem,%04d,pix.fld
		set fwhm = $outdir/fwhm,imp,osem.txt
	else			# EM
		set alg = em,1
		set penal = -

		set out = $outdir/em,%04d.fld
		set out_imp = $outdir/em,imp,%04d.fld
		set imp = $outdir/imp,em,%04d.fld
		set piximp = $outdir/imp,em,%04d,pix.fld
		set fwhm = $outdir/fwhm,imp,em.txt
	endif

	set method = @$niter@$alg@$penal

	# Add an impulse to phantom
	if !(-e $phantomP_imp) then
		set iximp = `echo "71 * $m/2" | bc`
		set iyimp = `echo "71 * $m/2" | bc`
		set izimp = `echo "15 * $m/2" | bc`
		op point3 t0 $nxP $nyP $nzP 0.1 $iximp $iyimp $izimp
		op add $phantomP_imp $phantomP t0

		ji t0
		ji $phantomP_imp
	endif

	# Generate projections using 3D SPECT system model
	if !(-e $proj_imp) then
		i proj3 t0 $phantomP_imp $ftypeP $maskP
		op nonlin max t0 t0 0 0		# Set negatives to zero 

		if ($down == yes) then		# Downsample projection views
			op div t0 t0 - $m	# !!!
			op sample mean $proj_imp t0 $m $m
		else
			cp t0 $proj_imp
		endif

		js $proj_imp
		op comp $proj_imp $proj
	endif

	set ybi_imp	= $indir/ybi,imp,$count.fld	# Scatter-free measurements means
	set ci_imp	= $indir/ci,imp,$count.fld	# Calibration factors

	if !(-e $ybi_imp) then
		op sim calib $ybi_imp $ci_imp $proj_imp - 0.0 $count 1	# (seed)

		op range $ci_imp $ybi_imp
		js $ci_imp
		js $ybi_imp
		op comp $ci_imp $ci
		op comp $ybi_imp $ybi
	endif

        set init = -$nx,$ny,$nz,1				# Uniform image
#set init = $fbpdir/fbp,hanning,0.6,0,1.fld

	# Iterative reconstruction of phantom w/ and w/o impulse
	if !(-e `printf $out $niter`) then
#		i -chat 20 empl3 $out $init $proj - 1 -shift 1 $ftype $mask \
		i -chat 20 empl3 $out $init $ybi $ci 1 - $scatmean $ftype $mask \
			$method $saver $obj 1 1e9 0 -
	endif
	if !(-e `printf $out_imp $niter`) then
		i -chat 20 empl3 $out_imp $init $ybi_imp $ci_imp 1 - $scatmean $ftype $mask \
			$method $saver $obj 1 1e9 0 -
	endif

#        op range $out $out_imp $phantom

	# Subtract to get impulse response
	set iter = $nsave
	while ($iter <= $niter)
		if !(-e `printf $imp $iter`) then
			op sub t0 `printf $out_imp $iter` `printf $out $iter`
			op transpose `printf $imp $iter` t0 1,2
		else
			op transpose t0 `printf $imp $iter` 1,2
		endif

		if ($btype != conv1sym) then	# Transform output from blob domain to pixel domain
			if !(-e `printf $piximp $iter`) then
				op sub t0 `printf $out_imp $iter` `printf $out $iter`
				echo "in = fld_read('t0'); out = blob_display(in, $sx, $R, $alpha, $kbm); mat_write('t0', out, '-nocheck');" | matlab -nodesktop -nosplash
				op transpose `printf $piximp $iter` t0 1,2
			else
				op transpose t0 `printf $piximp $iter` 1,2
			endif
			set fwhm = $outdir/fwhm,imp,$nsubset,$l2b,pix.txt
		endif

		# Print impulse response FWHM
		#op neg0 t0 t0
		op -chat 0 fwhm3 t0 usemax .5 "$iter" >> $fwhm

		set iter = `echo "$iter + $nsave" | bc`
	end
	#ji t0 				# Display in x-y direction
	ji `printf $imp $niter`		# Display in x-z direction

	exit
endif

#  OSPS  ##############################################################
set nsubset = 10
	set nsubset = 6
echo "Number of views per subset = "`echo "$na / $nsubset" | bc`

set niter = 25
	set niter = 1500
set l2bx = -8
set l2bz = -8
set l2b = $l2bx,$l2bz

set alg = ospsc,pc,$nsubset,$na,1,0.1
set penal = 3d,$l2b,quad,5,-
set method = @$niter@$alg@$penal

set ospsout = $outdir/os,$nsubset,$l2b,%04d.fld
set fwhm = $outdir/fwhm,xhat,$nsubset,$l2b.txt

if $arg == osps then
        set init = -$nx,$ny,$nz,1				# Uniform image
	if !(-e `printf $ospsout $niter`) then
		i -chat 20 empl3 $ospsout $init $ybi $ci 1 - $scatmean $ftype $mask \
			$method $saver $obj 1 1e9 0 -
	endif

	if ($btype != conv1sym) then	# Transform output from blob domain to pixel domain
		set pixout = $outdir/os,$nsubset,$l2b,%04d,pix.fld

		if !(-e `printf $pixout $niter`) then
			echo "for it = [${nsave}:${nsave}:$niter], in = fld_read(sprintf('$ospsout', it)); out = blob_display(in, $sx, $R, $alpha, $kbm); mat_write(sprintf('$pixout', it), out); end" | matlab -nodesktop -nosplash
		endif

		set ospsout = $pixout
		set fwhm = $outdir/fwhm,xhat,$nsubset,$l2b,pix.txt
	endif

	ji `printf $ospsout $niter`

	set iter = $nsave
	while ($iter <= $niter)
		# Print FWHM of reconstructed image
		op -chat 0 fwhm3 `printf $ospsout $iter` - .5 "$iter" >> $fwhm

		set iter = `echo "$iter + $nsave" | bc`
	end

	exit
endif

#  EM  ################################################################
set niter = 10
set alg = em,1
set penal = -
set method = @$niter@$alg@$penal

set emout = $outdir/em,%02d.fld

if $arg == em then
        i -chat 20 empl3 $emout $init $proj - 1 -shift 1e-5 $ftype $mask \
		$method $saver $obj 1 1e9 0 -

        op range $emout $phantom
        ji $emout
	exit
endif

#  OSEM  ##############################################################
set niter = 20
set alg = osemc,fast,$nsubset,$na,1.0
set penal = -
set method = @$niter@$alg@$penal

set osemout = $outdir/osem,$nsubset,%02d.fld

if $arg == osem then
        i -chat 20 empl3 $osemout $init $proj - 1 -shift 1e-5 $ftype $mask \
		$method $saver $obj 1 1e9 0 -

        op range `printf "$osemout" $niter` $phantom
        ji `printf "$osemout" $niter`
	exit
endif


#  Noisy simulations ##################################################
#  Loop over noisy measurement realizations and reconstruct with FBP/OSPS/OSEM
if ($arg == osemn) then		# For OSEM

	if ($#argv < 2) then
		echo "Usage: do.sh osemn {post|iter} [mean|resol|bias|noise|std]"
		exit 1
	endif

	set nsubset	= 6
	set pref	= osem

	if ($argv[2] == post) then	# For varying post-filtering

		set niter	= 1000
		set nsave	= $niter
		set saver	= -

	else if ($argv[2] == iter) then	# For varying iteration

		set niter	= 50
		set nsave	= 5
		#set niter	= 4
		#set nsave	= 1
		set saver	= save,$nsave
	endif

else if ($arg == ospsn) then	# For OSPS

	if ($#argv < 2) then
		echo "Usage: do.sh ospsn {beta|post|iter} [mean|resol|bias|noise|std]"
		exit 1
	endif

	set nsubset	= 10
	set pref	= os

	# The (l2bx,l2bz) pairs were chosen to give almost spherical impulse response
	if ($bl == uhe) then
		if ($argv[2] == beta) then	# For varying beta
			set	niter	= 100
			set	nsave	= 25
			set	saver	= save,$nsave

			set l2bx =	(-8	-7	-6	-5	-4	-3	-2	-1	0	1	2	3)
			set l2bz =	(-2.5	-2.5	-2.5	-2	-1.5	-1	-0.5	0.5	1.5	2.5	3.5	4.5)
			set nb = 12
		endif

	else if ($bl == tc) then
		if ($argv[2] == beta) then	# For varying beta
			set	niter	= 100
			set	nsave	= 25
			set	saver	= save,$nsave

			set l2bx =	(-8	-7	-6	-5	-4	-3	-2	-1	0	1	2	3)
			set l2bz =	(0	0	0	0	0	0	0.5	1	1.5	2.5	3.5	4.5)
			set nb = 9
		else if ($argv[2] == post) then	# For single beta, varying post-filtering
			set	niter	= 1500
			set	nsave	= $niter
			set	saver	= -

			set l2bx =	(-8)
			set l2bz =	(-3)
			set nb = 1
		else if ($argv[2] == iter) then	# For single beta, varying iteration
			set	niter	= 50
			set	nsave	= 5
			#set	niter	= 4
			#set	nsave	= 1
			set	saver	= save,$nsave

			set l2bx =	(-8)
			set l2bz =	(-3)
			set nb = 1
		endif
	endif
endif

if ($btype == conv1sym) then
	set suff = $niter
else
	set suff = $niter,pix
endif

if ($argv[1] == ospsn) then	# For OSPS w/ single beta
	if ($nb == 1)	set suff = $l2bx,$l2bz,$suff
endif

set seed1 = 0
set seed2 = 200

if (0 && ($#argv < 3)) then

	set seed = $seed1
	while ($seed < $seed2)

		set seed = `printf '%02d' $seed`

		# Generate noisy measurements
		set yi = $indir/yi,$seed.fld

		if (!(-e $yi) || $arg == yi) then
			op sim pet $yi - $ybi 1 $seed - $scatpc 1

			op range $yi
			js $yi
		endif

		# FBP reconstruction of noisy projections
		set	fbpwin	= hanning,0.6,0,1
		set	fbpn	= $fbpdir/fbp,$seed,$fbpwin.fld

		if ($arg == fbpn || !(-e $fbpn)) then
			op sub $tmp/t0 $yi - $scatmean	# ideal scatter correction
			op div $tmp/t0 $tmp/t0 $ci	# scan time correction
	#		op comp $tmp/t0 $proj
			op transpose $tmp/t0 $tmp/t0 1,2

			set cut = 0.5
			set ic = 0
			set nc = 11

			while ($ic < $nc)
				set fbpwin = hanning,$cut,0,1
				#set fbpwin = boxcar,1,1
				set fbpn = $fbpdir/fbp,$seed,$fbpwin.fld

				op fbp $tmp/t1 $tmp/t0 user,t $nx $ny $fbpwin 1 1 0  0 0 $orbit $ostart
				op mul $fbpn $tmp/t1 $chang

				op range $phantom $fbpn
				ji $fbpn

				set cut = `echo "$cut + 0.05" | bc`
				set cut = `printf '%g' $cut`
				@ ic += 1
			end
		endif

		if ($#argv > 1) then
			if ($argv[2] == iter) then	# For varying iteration
				set	init	= -$nx,$ny,$nz,$backg	# Uniform with true soft tissue value
			else						# For varying beta or post-filtering
				set	fbpwin	= hanning,0.6,0,1
				set	fbpn	= $fbpdir/fbp,$seed,$fbpwin.fld
				set	init	= $fbpn			# FBP
			endif
		endif

		# Regularized iterative reconstruction of noisy projections
		if ($arg == ospsn) then

			set ib = 1
			while ($ib <= $nb)

				set l2b = $l2bx[$ib],$l2bz[$ib]
				echo "seed = $seed, l2b = $l2b @ " `date`

				set	alg	= ospsc,pc,$nsubset,$na,1,0.1
				set	penal	= 3d,$l2b,quad,5,-
				set	method	= @$niter@$alg@$penal
				set	obj	= 0

				set	os0out	= $outdir/os0,$seed,$nsubset,$l2b.fld
				set	osout	= $outdir/os,$seed,$nsubset,$l2b,%02d.fld
				if ($saver == -) then
					set osout = `printf "$osout" $niter`
				endif

				#set init0 = $fbpout
				if (($seed == -1000) && (!(-e $os0out))) then
        				i -chat 20 empl3 $os0out $init0 $ybi - 1 -shift 1 $ftype $mask \
						$method $saver $obj 1 1e9 0 -
				endif

				if (!(-e `printf "$osout" $niter`)) then
					i -chat 20 empl3 $osout $init $yi $ci 1 - $scatmean $ftype $mask \
						$method $saver $obj 1 1e9 0 -
				endif

				ji `printf "$osout" $niter`

				if ($btype != conv1sym) then	# Transform output from blob domain to pixel domain
					set pixout = $outdir/os,$seed,$nsubset,$l2b,%02d,pix.fld

					echo "for it = [${nsave}:${nsave}:$niter], in = fld_read(sprintf('$osout', it)); out = blob_display(in, $sx, $R, $alpha, $kbm); mat_write(sprintf('$pixout', it), out); end" | matlab -nodesktop -nosplash

					ji `printf "$pixout" $niter`
				endif

				@ ib += 1
			end
		endif

		# Unregularized iterative reconstruction of noisy projections
		if ($arg == osemn) then

			echo "seed = $seed @ " `date`

			set	alg	= osemc,fast,$nsubset,$na,1.0
			set	penal	= -
			set	method	= @$niter@$alg@$penal
			set	obj	= 0

			set	osout	= $outdir/osem,$seed,$nsubset,%02d.fld
			if ($saver == -) then
				set osout = `printf "$osout" $niter`
			endif

			if (!(-e `printf "$osout" $niter`)) then
				i -chat 20 empl3 $osout $init $yi $ci 1 - $scatmean $ftype $mask \
					$method $saver $obj 1 1e9 0 -
			endif

			ji `printf "$osout" $niter`

			if ($btype != conv1sym) then	# Transform output from blob domain to pixel domain
				set pixout = $outdir/osem,$seed,$nsubset,%02d,pix.fld

				echo "for it = [${nsave}:${nsave}:$niter], in = fld_read(sprintf('$osout', it)); out = blob_display(in, $sx, $R, $alpha, $kbm); mat_write(sprintf('$pixout', it), out); end" | matlab -nodesktop -nosplash

				ji `printf "$pixout" $niter`
			endif
		endif

		@ seed += 1
	end

endif

# Post-filter noisy OSPS/OSEM reconstructions #########################
if ($#argv == 2) then
	if ($argv[2] == post) then

		set out = $pref,%02d,$nsubset,$suff.fld
	
		echo "dir='$outdir/'; infile='$out'; real=[${seed1}:$seed2-1]; fw=[0.5:0.5:7.0]; postfilt" \
			| matlab -nojvm -nosplash
	endif
endif

# Average realizations ################################################
if ($#argv > 2) then
	if ($argv[3] == mean) then

	if ($argv[2] == beta) then		# For varying beta

		while ($ib <= $nb)

			set l2b = $l2bx[$ib],$l2bz[$ib]

			set out = $outdir/os,*[0-9][0-9],$nsubset,$l2b,$suff.fld
			set meanout = $outdir/os,mean,$nsubset,$l2b,$suff.fld
			set stdout = $outdir/os,std,$nsubset,$l2b,$suff.fld

			op stat4 $meanout $stdout + $out

			@ ib += 1
		end

	else if ($argv[2] == post) then		# For varying post-filtering

		set fw = 0.5
		set nb = 14
		set ib = 0

		set out = $outdir/$pref,*[0-9][0-9],$nsubset,$suff.fld
		set meanout = $outdir/$pref,mean,$nsubset,$suff.fld
		set stdout = $outdir/$pref,std,$nsubset,$suff.fld

		op stat4 $meanout $stdout + $out

		while ($ib < $nb)

			set out = $outdir/$pref,*[0-9][0-9],$nsubset,$suff,g$fw.fld
			set meanout = $outdir/$pref,mean,$nsubset,$suff,g$fw.fld
			set stdout = $outdir/$pref,std,$nsubset,$suff,g$fw.fld

			op stat4 $meanout $stdout + $out

			set fw = `echo "$fw + 0.5" | bc`
			set fw = `printf '%g' $fw`
			@ ib += 1
		end

	else if ($argv[2] == iter) then		# For varying iteration

		#set iter = 1
		set iter = 5

		#set ni = 4
		set ni = 10
		set ii = 0

		while ($ii < $ni)

			set iter = `printf '%02d' $iter`
			if ($btype == conv1sym) then
				set suff = $iter
			else
				set suff = $iter,pix
			endif
			if ($arg == ospsn) then
				set suff = $l2bx,$l2bz,$suff
			endif

			set out = $outdir/$pref,*[0-9][0-9],$nsubset,$suff.fld
			set meanout = $outdir/$pref,mean,$nsubset,$suff.fld
			set stdout = $outdir/$pref,std,$nsubset,$suff.fld

			op stat4 $meanout $stdout + $out

			set iter = `echo "$iter + 5" | bc`
			@ ii += 1
		end
	endif
	endif
endif

# Average realizations for FBP ########################################
if $arg == meanfbp then

	set cut = 0.5

	set nb = 11
	set ib = 0

	while ($ib < $nb)

		set fbpwin = hanning,$cut,0,1
		set out = $fbpdir/fbp,*[0-9][0-9],$fbpwin.fld
		set meanout = $fbpdir/fbp,mean,$fbpwin.fld
		set stdout = $fbpdir/fbp,std,$fbpwin.fld

		op stat4 $meanout $stdout + $out

		set cut = `echo "$cut + 0.05" | bc`
		set cut = `printf '%g' $cut`
		@ ib += 1
	end
endif

# Calculate resolution for OSPS/OSEM ##################################
if ($#argv > 2) then
	if ($argv[3] == resol) then

	set resol = $outdir/resol,$pref,$argv[2].txt

	if ($argv[2] == beta) then		# For varying beta

		while ($ib <= $nb)

			set l2b = $l2bx[$ib],$l2bz[$ib]

			set out = $outdir/os,mean,$nsubset,$l2b,$suff.fld

			echo "out = fopen('$resol', 'a'); in0 = fld_read('$phantom'); in = fld_read('$out'); fwhm_gauss_x, fprintf(out, '$l2b\t\t%f\n', fw(1)); fclose(out);" | matlab -nojvm -nosplash

			@ ib += 1
		end

	else if ($argv[2] == post) then		# For varying post-filtering

		if 0 then
			set out = $outdir/$pref,mean,$nsubset,$suff.fld

			echo "out = fopen('$resol', 'a'); in0 = fld_read('$phantom'); in = fld_read('$out'); fwhm_gauss_x, fprintf(out, '0\t\t%f\n', fw(1)); fclose(out);" | matlab -nojvm -nosplash

		endif

		set fw = 0.5
		set ib = 0
		set nb = 10

		while ($ib < $nb)
			set out = $outdir/$pref,mean,$nsubset,$suff,g$fw.fld

			echo "out = fopen('$resol', 'a'); in0 = fld_read('$phantom'); in = fld_read('$out'); fwhm_gauss_x, fprintf(out, '$fw\t\t%f\n', fw(1)); fclose(out);" | matlab -nojvm -nosplash

			set fw = `echo "$fw + 0.5" | bc`
			set fw = `printf '%g' $fw`
			@ ib += 1
		end

	else if ($argv[2] == iter) then		# For varying iteration

		#set iter = 2
		set iter = 5

		#set ni = 3
		set ni = 10
		set ii = 0

		while ($ii < $ni)

			set iter = `printf '%02d' $iter`
			if ($btype == conv1sym) then
				set suff = $iter
			else
				set suff = $iter,pix
			endif
			if ($arg == ospsn) then
				set suff = $l2bx,$l2bz,$suff
			endif

			set out = $outdir/$pref,mean,$nsubset,$suff.fld

			echo "out = fopen('$resol', 'a'); in0 = fld_read('$phantom'); in = fld_read('$out'); fwhm_gauss_x, fprintf(out, '$iter\t\t%f\n', fw(1)); fclose(out);" | matlab -nodesktop -nosplash

			set iter = `echo "$iter + 5" | bc`
			@ ii += 1
		end
	endif
	endif
endif

# Calculate resolution for FBP ########################################
if $arg == resolfbp then
	set resol = $fbpdir/resol,fbp.txt

	set cut = 0.5
	set nb = 11

	set ib = 0

	while ($ib < $nb)

		set fbpwin = hanning,$cut,0,1
		set out = $fbpdir/fbp,mean,$fbpwin.fld

		echo "out = fopen('$resol', 'a'); in0 = fld_read('$phantom'); in = fld_read('$out'); fwhm_gauss_x, fprintf(out, '$cut\t\t%f\n', fw(1)); fclose(out);" | matlab -nodesktop -nosplash

		set cut = `echo "$cut + 0.05" | bc`
		set cut = `printf '%g' $cut`
		@ ib += 1
	end
endif

# Compute average noise in lung ROI for OSPS/OSEM #####################
if ($#argv > 2) then
	if ($argv[3] == noise) then

	set noise = $outdir/noise,$pref,$argv[2].txt

	if ($argv[2] == beta) then		# For varying beta

		while ($ib <= $nb)

			set l2b = $l2bx[$ib],$l2bz[$ib]

			set std = $outdir/os,std,$nsubset,$l2b,$suff.fld

			op slice t0 $std - 0 7 0 0 0 0			# 8 slices
			op roi t1 lung 20 31 8 8 1 t0			# 8 x 8 ROI
			printf $l2b >> $noise
			op range t1 | perl -ne 'if (/^x sum mean/) {chop; @f = split; print "\t\t$f[4]\n"}' >> $noise

			@ ib += 1
		end

	else if ($argv[2] == post) then		# For varying post-filtering

		set std = $outdir/$pref,std,$nsubset,$suff.fld

		op slice t0 $std - 0 7 0 0 0 0			# 8 slices
		op roi t1 lung 20 31 8 8 1 t0			# 8 x 8 ROI
		printf 0 >> $noise
		op range t1 | perl -ne 'if (/^x sum mean/) {chop; @f = split; print "\t\t$f[4]\n"}' >> $noise

		set fw = 0.5
		set ib = 0
		set nb = 10

		while ($ib < $nb)

			set std = $outdir/$pref,std,$nsubset,$suff,g$fw.fld

			op slice t0 $std - 0 7 0 0 0 0			# 8 slices
			op roi t1 lung 20 31 8 8 1 t0			# 8 x 8 ROI
			printf $fw >> $noise
			op range t1 | perl -ne 'if (/^x sum mean/) {chop; @f = split; print "\t\t$f[4]\n"}' >> $noise

			set fw = `echo "$fw + 0.5" | bc`
			set fw = `printf '%g' $fw`
			@ ib += 1
		end

	else if ($argv[2] == iter) then		# For varying iteration

		#set iter = 2
		set iter = 5

		#set ni = 3
		set ni = 10
		set ii = 0

		while ($ii < $ni)

			set iter = `printf '%02d' $iter`
			if ($btype == conv1sym) then
				set suff = $iter
			else
				set suff = $iter,pix
			endif
			if ($arg == ospsn) then
				set suff = $l2bx,$l2bz,$suff
			endif

			set std = $outdir/$pref,std,$nsubset,$suff.fld

			op slice t0 $std - 0 7 0 0 0 0			# 8 slices
			op roi t1 lung 20 31 8 8 1 t0			# 8 x 8 ROI
			printf $iter >> $noise
			op range t1 | perl -ne 'if (/^x sum mean/) {chop; @f = split; print "\t\t$f[4]\n"}' >> $noise

			set iter = `echo "$iter + 5" | bc`
			@ ii += 1
		end
	endif
	endif
endif

# Compute average noise in lung ROI for FBP ###########################
if $arg == noisefbp then
	set noise = $fbpdir/noise,fbp.txt

	set cut = 0.5
	set nb = 11

	set ib = 0

	while ($ib < $nb)

		set fbpwin = hanning,$cut,0,1
		set std = $fbpdir/fbp,std,$fbpwin.fld

		op slice t0 $std - 0 7 0 0 0 0			# 8 slices
		op roi t1 lung 20 31 8 8 1 t0			# 8 x 8 ROI
		printf $cut >> $noise
		op range t1 | perl -ne 'if (/^x sum mean/) {chop; @f = split; print "\t\t$f[4]\n"}' >> $noise

		set cut = `echo "$cut + 0.05" | bc`
		set cut = `printf '%g' $cut`
		@ ib += 1
	end
endif

# Compute STD of heart ROI activity for OSPS/OSEM #####################
set voi = $indir/voi,heart.fld
if (!(-e $voi)) then
	echo "x = fld_read('$phantom'); h = (x > 3.5); mat_write('$voi', h);" \
		| matlab -nodesktop -nosplash
endif
op mul t20 $phantom $voi
set true = `op range t20 | grep '^x sum mean' | awk '{print $4}'`
set true = `printf '%g' "$true"`

if ($#argv > 2) then
	if ($argv[3] == std) then

	set nseed = `echo "$seed2 - $seed1" | bc`
	op rep $tmp/t0 $voi $nseed

	set bias = $outdir/bias,$pref,$argv[2].txt
	set stdh = $outdir/stdheart,$pref,$argv[2].txt

	cp /dev/null $tmp/tmp.txt

	if ($argv[2] == beta) then		# For varying beta

		while ($ib <= $nb)

			set l2b = $l2bx[$ib],$l2bz[$ib]

			set out = $outdir/os,*[0-9][0-9],$nsubset,$l2b,$suff.fld

			op -chat 0 stack4 $tmp/t1 float $out
			op mul $tmp/t1 $tmp/t1 $tmp/t0
			op range-frame $tmp/t1 | grep '^x sum mean' | awk '{print $4}' >> $tmp/tmp.txt

			@ ib += 1
		end

		echo "bet = [$l2bx]; outm = '$bias'; outs = '$stdh'; nseed = $nseed; true = $true; heart;" \
			| matlab -nojvm -nosplash

	else if ($argv[2] == post) then		# For varying post-filtering

		set out = $outdir/$pref,*[0-9][0-9],$nsubset,$suff.fld

		op -chat 0 stack4 $tmp/t1 float $out
		op mul $tmp/t1 $tmp/t1 $tmp/t0
		op range-frame $tmp/t1 | grep '^x sum mean' | awk '{print $4}' >> $tmp/tmp.txt

		set fw = 0.5
		set ib = 0
		set nb = 14

		while ($ib < $nb)

			set out = $outdir/$pref,*[0-9][0-9],$nsubset,$suff,g$fw.fld

			op -chat 0 stack4 $tmp/t1 float $out
			op mul $tmp/t1 $tmp/t1 $tmp/t0
			op range-frame $tmp/t1 | grep '^x sum mean' | awk '{print $4}' >> $tmp/tmp.txt

			set fw = `echo "$fw + 0.5" | bc`
			set fw = `printf '%g' $fw`
			@ ib += 1
		end

		echo "bet = [0:0.5:7]; outm = '$bias'; outs = '$stdh'; nseed = $nseed; true = $true; heart;" \
			| matlab -nojvm -nosplash

	else if ($argv[2] == iter) then		# For varying iteration

		#set iter = 1
		set iter = 5

		#set ni = 4
		set ni = 10
		set ii = 0

		while ($ii < $ni)

			set iter = `printf '%02d' $iter`
			if ($btype == conv1sym) then
				set suff = $iter
			else
				set suff = $iter,pix
			endif
			if ($arg == ospsn) then
				set suff = $l2bx,$l2bz,$suff
			endif

			set out = $outdir/$pref,*[0-9][0-9],$nsubset,$suff.fld

			op -chat 0 stack4 $tmp/t1 float $out
			op mul $tmp/t1 $tmp/t1 $tmp/t0
			op range-frame $tmp/t1 | grep '^x sum mean' | awk '{print $4}' >> $tmp/tmp.txt

			set iter = `echo "$iter + 5" | bc`
			@ ii += 1
		end

		echo "bet = [5:5:50]; outm = '$bias'; outs = '$stdh'; nseed = $nseed; true = $true; heart;" \
		#echo "bet = [1:4]; outm = '$bias'; outs = '$stdh'; nseed = $nseed; true = $true; heart;" \
			| matlab -nojvm -nosplash
	endif
	endif
endif

# Compute STD and bias of heart ROI activity for FBP ##################
if $arg == stdfbp then

	set nseed = `echo "$seed2 - $seed1" | bc`
	op rep $tmp/t0 $voi $nseed

	set bias = $fbpdir/bias,fbp.txt
	set stdh = $fbpdir/stdheart,fbp.txt

	cp /dev/null $tmp/tmp.txt

	set cut = 0.5
	set nb = 11

	set ib = 0

	while ($ib < $nb)

		set fbpwin = hanning,$cut,0,1
		set out = $fbpdir/fbp,*[0-9][0-9],$fbpwin.fld

		op -chat 0 stack4 $tmp/t1 float $out
		op mul $tmp/t1 $tmp/t1 $tmp/t0
		op range-frame $tmp/t1 | grep '^x sum mean' | awk '{print $4}' >> $tmp/tmp.txt

		set cut = `echo "$cut + 0.05" | bc`
		set cut = `printf '%g' $cut`
		@ ib += 1
	end

	echo "bet = [.5:.05:1]; outm = '$bias'; outs = '$stdh'; nseed = $nseed; true = $true; heart;" \
		| matlab -nojvm -nosplash
endif

exit	##########################################################

