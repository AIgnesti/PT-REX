###########################
"""
Author: Alessandro Ignesti
Ver 2.0
Contacts: https://aignesti.github.io/

Point-to-point TRend EXtractor [PT-REX]

##########################

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

Run PT-REX.py in the same folder of the images and the region files
$ python3 PT-REX.py
