This directory has various data files used in the examples.

To keep the github repo small, the git version does not include the actual
data files.  Instead those files are stored in a separate repo:
https://github.com/JeffFessler/MIRTdata/

The files are fetched when needed from that data repo via `ir_get_data`
and stored in the `downloads/` folder.

Data files include:

* George Zubal's phantom from http://noodle.med.yale.edu/zubal/

* some slices of the NCAT phantom from UNC (or maybe now JHU)
`ncat,256,slice,140,ct,dens.fld`
see:
@a segars:02:sot ieee-t-ns 49 3 675-679 jun 2002
W P Segars  B M W Tsui
Study of the efficacy of respiratory gating in myocardial SPECT using the new 4-D NCAT phantom
http://doi.org/10.1109/TNS.2002.1039548

Please cite Paul Segar's paper(s) if you use that phantom image!

See subdirectories too.
