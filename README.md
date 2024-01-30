# PT-REX

Point-to-point TRend EXtractor
v 3.0

Investigating the spatial correlation between different emissions in an extended astrophysical source can provide crucial insights into their physical connection, hence it can be the key to understand the nature of the system. The point-to-point analysis of surface brightness is a reliable method to do such an analysis. The PT-REX code is designed to carry out these studies between different emissions in extended sources. For further information check [the paper](https://www.sciencedirect.com/science/article/pii/S1384107621001457?via%3Dihub).
<p align="center">
<img src="https://github.com/AIgnesti/PT-REX/blob/master/images/%20logo2.png" width=30% heigth=30%>
</p>

## Requirements

- Python 3.7
- matplotlib 3.5.3
- numpy 1.21.5
- astropy 4.2.1
- regions 0.5
- bces 1.0.3
- scipy 1.7.3
- nmmn 1.2.1


## Running PT-REX

### Launching the code

PT-REX analyzes two separate images in FITS format. Run the code in the same folder as the two images ans use the following keyword to set the parameters:

Required:
- -im1: Set IMAGE1 
- -im2: Set IMAGE2 
- -cel_w: Set cell width [arcsec] 
- -cel_h: Set cell height [arcsec] 
- -thr: Set IMAGE1 threshold [IMAGE1 units]
Optional:
- -lm: Set grid overlap with mask [0.0-0.99, default 0.5]
- -sm: Set Gaussian smoothing sigma size for IMAGE2 [arcsec]
- -h Print help

EXAMPLE: 
```bash
python PT-REX -im1 image1.fits -im2 image2.fits -cel_w 10.0 -cel_h 10.0 -thr 42.0
```
### Performing the analysis

- Left panel: IMAGE1. The silver contour show the threshold level
- Right panel: IMAGE2 

Interactive commands:
click  -->  release on Left panel: Define region of interest (ROI)
Then, press:
- W/w: Create rectangular grid in ROI
- H/h: Create hexagonal grid in ROI
- D/d: Add mask 
- I/i: Run PtP analysis with active cells. Surface brigthenss is computed as sum of pixel values divided by cell area.  Output: Correlation plot in out.jpg, maps in out_plot.png, data series in out.dat, and grid in DS9 format in out_grid.reg
- X/x: Recenter images
+/-: Increase/decrease cell size by 0.5 arcsec
- C/c: Clear grid
<p align="center">
<img src="https://github.com/AIgnesti/PT-REX/blob/master/images/out_plot.png" width=80% heigth=80%>
</p>

## Author
Author: Alessandro Ignesti

Contacts: https://aignesti.github.io/ | alessandro.ignesti@inaf.it


## Cite this work

 If PT-REX was helpful to your research an acknowledgment would be appreciated:
 -Cite the paper: 
 ```bibtex
@ARTICLE{2022NewA...9201732I,
       author = {{Ignesti}, A.},
        title = "{Introducing PT-REX, the point-to-point TRend EXtractor}",
      journal = {\na},
     keywords = {Techniques, Image processing, Methods, Observational, Statistical, Radio continuum, General, Astrophysics - Instrumentation and Methods for Astrophysics},
         year = 2022,
        month = apr,
       volume = {92},
          eid = {101732},
        pages = {101732},
          doi = {10.1016/j.newast.2021.101732},
archivePrefix = {arXiv},
       eprint = {2110.12720},
 primaryClass = {astro-ph.IM},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2022NewA...9201732I},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```
 
 -Indicate the URL of the GitHub repository as a footnote: \footnote{https://github.com/AIgnesti/PT-REX}



