# Point-to-point TRend EXtractor (PT-REX)

*v 3.0*

Investigating the spatial correlation between different emissions in an extended astrophysical source can provide crucial insights into their physical connection, hence it can be the key to understanding the nature of the system. The point-to-point analysis of surface brightness is a reliable method to do such an analysis. The PT-REX code is designed[^1] to carry out these studies between different emissions in extended sources. For further information check [the paper](https://www.sciencedirect.com/science/article/pii/S1384107621001457?via%3Dihub).

<p align="center">
<img src="https://github.com/AIgnesti/PT-REX/blob/master/images/%20logo2.png" width=30% heigth=30%  />
</p>

## Requirements
PT-REX is a Python 3 code which requires the following packages:

- matplotlib 3.5.3
- numpy 1.21.5
- astropy 4.2.1
- regions 0.5
- bces 1.0.3
- scipy 1.7.3
- nmmn 1.2.1

They can be installed by running:

```bash
python -m pip install -r requirements.txt
```

## Running PT-REX

### Launching the code

PT-REX analyzes two separate images in FITS format. Run the code in the same folder as the two images and use the following keywords to set the parameters:

*Required*:
- **-im1**: Set IMAGE1 
- **-im2**: Set IMAGE2 
- **-cel_w**: Set cell width [arcsec] 
- **-cel_h**: Set cell height [arcsec] 
- **-thr**: Set IMAGE1 threshold [IMAGE1 units]
  
*Optional*:

- **-lm**: Set the minimum emission overlap required for each cell [0.0-0.9, default 0.5]
- **-sm**: Set a Gaussian smoothing sigma size for IMAGE2 [arcsec]
- **-h** Print help

EXAMPLE: 
```bash
python PT-REX.py -im1 image1.fits -im2 image2.fits -cel_w 10.0 -cel_h 10.0 -thr 42.0
```
### Performing the analysis

- Left panel: IMAGE1. The silver contour indicates the threshold level.
- Right panel: IMAGE2 

Interactive commands:
click  -->  release on Left panel: Define Region of Interest (ROI)

Then, with the cursor on the left panel, press:

- **W/w**: Create a rectangular grid in ROI
- **H/h**: Create a hexagonal grid in ROI
- **D/d**: Mask ROI
- **I/i**: Run PtP analysis with active cells. Surface brightness $I_1$ and $I_2$ are computed as the sum of pixel values divided by cell area. The corresponding error, $\sigma_1$ and $\sigma_2$ are derived from the rms in each cell. The best-fit correlation log $I_2=k\cdot$ log $I_1+A$ is derived with the orthogonal BCES algorithm. An estimate of the best-fitting $k$ and the Spearman and Pearson ranks are presented in the legend of the scatter plot.
  
  Output:
  
  - out.jpg: Scatter plot with best-fit slope
  - out_plot.png: IMAGE1 and IMAGE2 with cells used in PtP analysis
  - out.dat: Data readout in the format $I_1$  $\sigma_1$  $I_2$  $\sigma_2$  
  - out_grid.reg: grid in DS9 format for IMAGE1
- **X/x**: Recenter images in ROI
- **+/-**: Increase/decrease cell size by 0.5 arcsec
- **C/c**: Clear grid and masks
<p align="center">
<img src="https://github.com/AIgnesti/PT-REX/blob/master/images/out_plot.png" width=80% heigth=80%>

 <img src="https://github.com/AIgnesti/PT-REX/blob/master/images/out.jpg" width=50% heigth=50%>
</p>

## Author
Author: Alessandro Ignesti

Contacts: https://aignesti.github.io/ | alessandro.ignesti@inaf.it


## Cite this work

 If PT-REX was helpful to your research an acknowledgment would be appreciated:
 
 - Cite the paper: 
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
 
 - Indicate the URL of the GitHub repository as a footnote: 
 ```latex
 \footnote{\url{https://github.com/AIgnesti/PT-REX}}
```
- Refer to the [Astrophysics Source Code Library](https://ascl.net/2110.021)

## License
See LICENSE .

[^1]: I study stuff in galaxy clusters, I am not a professional programmer, and this code was developed for research purposes only. So be aware that the code is not completely stable and use it at your own risk. If you encounter any issues, please consider them as unplanned features, and if you would like to contribute, find things that are broken, or have any suggestions for this work, you can contact me at alessandro.ignesti@inaf.it .





