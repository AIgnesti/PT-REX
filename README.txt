###########################
"""
Author: Alessandro Ignesti
Ver 2.2
Contacts: https://aignesti.github.io/ | alessandro.ignesti@inaf.it

Point-to-point TRend EXtractor [PT-REX]
For further information check https://www.sciencedirect.com/science/article/pii/S1384107621001457?via%3Dihub
##########################
What's new:
-v2.2 : Bug fix. 
-v2.1 : Improved fitting algorithm to handle datasets with large dynamical range (> 2 orders of magnitude);
-v2.0 : First public release.
##########################
REQUIRMENTS:
-Python 3.*
-CASA 6.0 Python modular version (https://casa.nrao.edu/casadocs/casa-6.1.0/usingcasa/obtaining-and-installing)
[ CASA 6.0 or lower is strongly advised to use the imviewer task in parallel with PT-REX]
-matplotlib.pyplot as plt
-astropy
-numpy 
-scipy
-bces (https://github.com/rsnemmen/BCES)
##########################

DATA PREPARATION:
The input maps have to re-named as:
-radio map: [name]-radio.fits
-X-rays instensity (counts) map: [name]-img.fits
-X-rays background map: [name]-bmap.fits
-X-rays exposure map: [name]-emap.fits

The region and mask box files must be in CASA region format (.crtf).
##########################

USAGE:

1) Run PT-REX.py in the same folder of the images and the region files:

$ python3 PT-REX.py

2) Provide the requested inputs, that are fitting method and calibration error. They are going to be used during the analysis, and they can be updated by using the corresponding task [r]. 

$ Calibration error of the radio image [es. 0.05]: 0.1
$ Fitting method [LS | BCES_ort | BCES_bisec | LinMix]: LS

3) Select and run a task from the list:

$TASK LIST
$-Create a new mask.image from a mask.reg [1]
$-Create a new J2000 mesh.crtf for the radio map [2]
$-Single mesh analysis [3]
$-Monte-Carlo analysis [4]
$UTILITIES:
$--Reload images and configuration [r]
$--quit [q]
