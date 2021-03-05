2001-06-12
This directory has MR data from Doug Noll
for testing PWLS type of phase fitting.

In `phfit.mat` there is a complex image (variable: `m1`)
where the phase represents the magnetic field inhomogeneity and the
magnitude is just the image.

The echo-time difference associated with this scan is 2ms.

I don't remember what the file `phfile.mat` is.

For the benefit of those who cannot read .mat files,
I also include the big endian 64x64 float files:
* `phfit_im.raw`
* `phfit_re.raw`

and the corresponding `.fld` header files (read with `fld_read`):
* `phfit_im.fld`
* `phfit_re.fld`

These phase data were used in my
2005 IEEE T-SP paper on Toeplitz based MR image reconstruction
[http://doi.org/10.1109/TSP.2005.853152]


2008-01-12
for Amanda Funai's field map estimation paper
[http://doi.org/10.1109/TMI.2008.923956]

```matlab
load biggerbrain
fld_write('mag128.fld', single(bigger_mag))
fld_write('fieldmap128.fld', single(bigger_w/2/pi))
```
