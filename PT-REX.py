#!usr/bin/python
###########################
#Author: Alessandro Ignesti
#Point-to-point TRend EXtractorn V. 3.0
#For reference https://www.sciencedirect.com/science/article/pii/S1384107621001457 
###########################

from matplotlib.widgets import RectangleSelector
from matplotlib import colors
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
from matplotlib.colors import LogNorm, SymLogNorm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.patches as patches
from regions.core import PixCoord
from regions import PixCoord, PolygonSkyRegion, PolygonPixelRegion
from regions import RectangleSkyRegion, RectanglePixelRegion
from astropy.coordinates import SkyCoord  # High-level coordinates
from astropy.wcs import WCS
from astropy.wcs import utils
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
import bces.bces as BCES
import scipy.ndimage as ndimage
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
from regions import Regions
import astropy.units as u
import nmmn.stats
import sys
import os
from matplotlib import cm
from scipy import stats
import warnings
warnings.filterwarnings("ignore")

cmpa = cm.get_cmap('Blues', 128)
newcolors = cmpa(np.linspace(0, 1, 10))
newcmp = ListedColormap(newcolors)



region_type='rectangle'
stepx = "none"
stepy ='none'
thresh='none'
sm_g='none'
limit=0.5
image_file1='none'#'JW100new_Ha_5x5.fits'
image_file2='none'#'16136-img.fits'
if("-cel_w" in  sys.argv):
	stepx = float(sys.argv[sys.argv.index("-cel_w") + 1]  )  

if("-cel_h" in  sys.argv):
	stepy = float(sys.argv[sys.argv.index("-cel_h") + 1] )
if("-sm" in  sys.argv):
	sm_g= float(sys.argv[sys.argv.index("-sm") + 1] )   

if("-thr" in  sys.argv):
	thresh = float(sys.argv[sys.argv.index("-thr") + 1] )   

if("-im1" in  sys.argv):
	image_file1 = sys.argv[sys.argv.index("-im1") + 1] 
if("-im2" in  sys.argv):
	image_file2 = sys.argv[sys.argv.index("-im2") + 1] 
if("-lm" in sys.argv):
	limit= float(sys.argv[sys.argv.index("-lm") + 1] )
if("-h" in  sys.argv or len(sys.argv)==1):
	print('---------------------')
	print('      PT-REX 3.0     ')
	print('---------------------')
	print('Command line:')
	print('-im1: Set IMAGE1 ')
	print('-im2: Set IMAGE2 ')
	print('-cel_w: Set cell width [arcsec] ')
	print('-cel_h: Set cell height [arcsec] ')
	print('-thr: Set IMAGE1 threshold [IMAGE1 units]')
	print('Optional:')
	print('-lm: Set grid overlap with mask [0.0-0.99, default 0.5]')
	print('-sm: Set Gaussian smoothing sigma size for IMAGE2 [arcsec]')
	print('-h Print help')
	print('EXAMPLE: python PTREX_3 -im1 image1.fits -im2 image2.fits -cel_w 10.0 -cel_h 10.0 -thr 42.0')
	print('---------------------')
	exit()
if('-pt-rex' in sys.argv):
	print(r'       __            ')
	print(r'      /oo\           ')
	print(r'     |    |          ')
	print(r' ^^  (vvvv)   ^^     ')
	print(r' \\  /\__/\  //      ')
	print(r'  \\/      \//       ')
	print(r'   /        \        ')
	print(r'  |          |    ^  ')
	print(r'  /          \___/ | ')
	print(r' (            )     |')
	print(r'  \----------/     / ')
	print(r'    //    \\_____/   ')
	print(r'   W       W         ')
	exit()
#from matplotlib.nxutils import points_inside_poly
celle=[]


def hexagon(xc,yc,r):
	p1=[xc+r*np.sqrt(3.)/2.,yc+r/2.]
	p2=[xc,yc+r]
	p3=[xc-r*np.sqrt(3.)/2.,yc+r/2.]
	p4=[xc-r*np.sqrt(3.)/2.,yc-r/2.]
	p5=[xc,yc-r]
	p6=[xc+r*np.sqrt(3.)/2.,yc-r/2.]
	vertices = PixCoord(x=[p1[0],p2[0],p3[0],p4[0],p5[0],p6[0]],y=[p1[1],p2[1],p3[1],p4[1],p5[1],p6[1]])

	return PolygonPixelRegion(vertices=vertices)




def thresh_mapper(thresh,img):
	global thresh_map
	thresh_map=np.zeros_like(img)
	thresh_map[img>thresh]=1.
	return thresh_map
def thr_map_update(x1,y1,x2,y2,thresh_map_l):
	reg_d = RectanglePixelRegion(PixCoord((x2+x1)/2.,(y2+y1)/2.), width=x2-x1,height=y2-y1)
	mask_d = reg_d.to_mask()
	thresh_map_l=np.subtract(thresh_map_l,mask_d.to_image(thresh_map.shape))
	

	return reg_d,thresh_map_l

def pow(x,a,b):
	return a*x+b
def func(x): return x[1]*x[0]+x[2]
def drop2axes(filename, outname):
	hdu = fits.open(filename)[0]
	for kw in "CTYPE", "CRVAL", "CRPIX", "CDELT", "CUNIT","CROTA":
		try:
			for n in 3, 4:
				hdu.header.remove(f"{kw}{n}")
		except:
			continue
	fits.writeto(outname, hdu.data[0,0], hdu.header, overwrite=True)


def line_select_callback(eclick, erelease):
	'eclick and erelease are the press and release events'
	global x1,x2,y1,y2
	x1, y1 = eclick.xdata, eclick.ydata
	x2, y2 = erelease.xdata, erelease.ydata
	
	print("(%3.2f, %3.2f) --> (%3.2f, %3.2f)" % (x1, y1, x2, y2))
	

def toggle_selector(event):
	global thresh_map
	global region_type
	global stepx
	global stepy
	if event.key in ['Q', 'q'] and toggle_selector.RS.active:
		print('Bye')
		toggle_selector.RS.set_active(False)
		for file in os.listdir("./"):
			if file.endswith("_2ax_temp.fits"):
				os.system('rm '+file)
	
	if event.key in ['+'] and toggle_selector.RS.active:
		stepx=stepx+0.5
		stepy=stepy+0.5
		print('----')
		print('Cell width: '+str(stepx)+' arcsec')
		print('Cell height: '+str(stepy)+' arcsec')
		print('----')
	if event.key in ['-'] and toggle_selector.RS.active:
		stepx=stepx-0.5
		stepy=stepy-0.5
		print('----')
		print('Cell width: '+str(stepx)+' arcsec')
		print('Cell height: '+str(stepy)+' arcsec')
		print('----')


	if event.key in ['H', 'h'] and toggle_selector.RS.active:
		region_type='hexagon'
		print('- Hexagon binning -')
		r=np.sqrt(stepx**2+stepy**2)/scale1/2.
		x_grid=np.arange(x1,x2,r*np.sqrt(3.)/2.)
		xi=x1
		yi=y1
		n=0.

		while xi<x2:
			xi=x1+n*r*np.sqrt(3.)/2.		
			yi=y1
			m=0.
			while yi<y2-3.*r:
				if (n % 2) == 0:
					yi=y1+m*3.*r
				else:
		
					yi=y1+3.*r*m+1.5*r
				center = PixCoord(xi, yi)
				reg1=hexagon(xi,yi,r)
				mask = reg1.to_mask(mode='center')
	
				weighted_data_1d = mask.get_values(thresh_map)
			if np.sum(weighted_data_1d)>limit*np.sqrt(3.)*3./2.*r*r:
					celle.append(center)
					
					wx, wy = w1.wcs_pix2world(center.x,center.y, 1)
#
					px, py = w2.wcs_world2pix(wx, wy, 1)
#
					reg1.plot(ax=ax1, facecolor='blue', edgecolor='dimgrey', lw=2)
					reg1.plot(ax=ax1,color='gold',fill=True,alpha=0.5)
					reg2 = hexagon(px,py,r*scale1/scale2)# RectanglePixelRegion(PixCoord(px,py), width=stepx/scale2,height=stepy/scale2)
					reg2.plot(ax=ax2, facecolor='blue', edgecolor='dimgrey', lw=2)
					reg2.plot(ax=ax2,color='gold',fill=True,alpha=0.5)
				
					
					
				m=m+1
				
			n=n+1
			
		plt.draw()
	if event.key in ['X', 'x'] and toggle_selector.RS.active:
		ax1.set_xlim(x1,x2)
		ax1.set_ylim(y1,y2)


		wx1, wy1 = w1.wcs_pix2world(x1,y1, 1)
		wx2, wy2 = w1.wcs_pix2world(x2,y2, 1)
		px1, py1 = w2.wcs_world2pix(wx1, wy1, 1)
		px2, py2 = w2.wcs_world2pix(wx2, wy2, 1)

		ax2.set_xlim(px1,px2)
		ax2.set_ylim(py1,py2)


		plt.draw()
	
	if event.key in ['D', 'd'] and toggle_selector.RS.active:
		reg_d,thresh_map=thr_map_update(x1,y1,x2,y2,thresh_map)
		reg_d.plot(ax=ax1, color='red',alpha=0.4, fill=True)
		thresh_map[thresh_map<0.]=0.
		plt.draw()

	if event.key in ['W','w'] and toggle_selector.RS.active:
		print('- Rectangle binning -')
		region_type='rectangle'
		x_grid=np.arange(x1,x2,stepx/scale1)
		y_grid=np.arange(y1,y2,stepy/scale1)
		x_cen=[(x_grid[k]+x_grid[k+1])/2. for k in range(0,len(x_grid)-1)]
		y_cen=[(y_grid[k]+y_grid[k+1])/2. for k in range(0,len(y_grid)-1)]
		
		for xi in x_cen:
			for yi in y_cen:
				center = PixCoord(xi, yi)

				reg1 = RectanglePixelRegion(center, width=stepx/scale1, height=stepy/scale1)

				mask = reg1.to_mask(mode='center')


				weighted_data_1d = mask.get_values(thresh_map)
				if np.sum(weighted_data_1d)>limit*stepx*stepy/scale1**2:
					celle.append(center)
					
					wx, wy = w1.wcs_pix2world(center.x,center.y, 1)

					px, py = w2.wcs_world2pix(wx, wy, 1)

					reg1.plot(ax=ax1, facecolor='blue', edgecolor='dimgrey', lw=2)
					reg1.plot(ax=ax1, color='gold',fill=True,alpha=0.5)
					reg2 = RectanglePixelRegion(PixCoord(px,py), width=stepx/scale2,height=stepy/scale2)
					reg2.plot(ax=ax2, facecolor='blue', edgecolor='dimgrey', lw=2)
					reg2.plot(ax=ax2, color='gold',fill=True,alpha=0.5)
		plt.draw()
	if event.key in ['I','i'] and toggle_selector.RS.active:
		print('- Grid reading -')
		ax1.patches.clear()
		ax2.patches.clear()
		
		grid=open('out_grid.reg', 'w')
		grid.write('#Region file format: DS9 version 4.1\n')
		grid.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
		
		if region_type=='rectangle':
			grid.write('fk5\n')
		if region_type=='hexagon':
			grid.write('image\n')
		S1=[]
		S2=[]
		e_s1=[]
		e_s2=[]
		
		for center in celle:

			wx, wy = w1.wcs_pix2world(center.x,center.y, 1)
			
			px, py = w2.wcs_world2pix(wx, wy, 1)
			if region_type=='rectangle':
				reg1 = RectanglePixelRegion(center, width=stepx/scale1,height=stepy/scale1)
				reg2 = RectanglePixelRegion(PixCoord(px,py), width=stepx/scale2,height=stepy/scale2)
				area=stepx*stepy
				grid.write('box('+str(wx)+','+str(wy)+','+str(stepx)+'",'+str(stepy)+'",0.0)\n')

			if region_type=='hexagon':
				r=np.sqrt(stepx**2+stepy**2)/scale1/2.
				reg1=hexagon(center.x,center.y,r)
				reg2=hexagon(px,py,r/scale2*scale1)

				area=3./2.*np.sqrt(3.)*r*r*scale1**2
				grid.write('polygon('+str(reg1.vertices[0].x)+','+str(reg1.vertices[0].y)+','+str(reg1.vertices[1].x)+','+str(reg1.vertices[1].y)+','+str(reg1.vertices[2].x)+','+str(reg1.vertices[2].y)+','+str(reg1.vertices[3].x)+','+str(reg1.vertices[3].y)+','+str(reg1.vertices[4].x)+','+str(reg1.vertices[4].y)+','+str(reg1.vertices[5].x)+','+str(reg1.vertices[5].y)+')\n')
			mask1 = reg1.to_mask(mode='center')

			mask2 = reg2.to_mask(mode='center')
			serie_data_1 = mask1.get_values(image_data1)
			serie_data_2 = mask2.get_values(image_data2)
			v1=np.nansum(serie_data_1)/area
			v2=np.nansum(serie_data_2)/area
			e_v1=np.sqrt(np.nanmean(np.square(serie_data_1)))
			e_v2=np.sqrt(np.nanmean(np.square(serie_data_2)))

			if v1>0. and v2>0. and np.isnan(v1)==False and np.isnan(v2)==False and e_v1>0. and e_v2>0. and np.isnan(e_v1)==False and np.isnan(e_v2)==False:

				S1.append(v1)
				S2.append(v2)
				e_s1.append(e_v1)
				e_s2.append(e_v2)
				reg1.plot(ax=ax1, facecolor='blue', edgecolor='dimgrey', lw=2)
				reg1.plot(ax=ax1, color='gold',fill=True,alpha=0.5)
				reg2.plot(ax=ax2, facecolor='blue', edgecolor='dimgrey', lw=2)
				reg2.plot(ax=ax2, color='gold',fill=True,alpha=0.5)
				
		ran=len(S1)
		grid.close()
		plt.draw()
		plt.savefig('out_plot.png',bbox_inches='tight',pad_inches=0.1,dpi=200)
		
		
		x_range=np.linspace(np.min(S1),np.max(S2),100)
		with open(r'out.dat', 'w') as fp:
			for i in range(0,ran):
				# write each item on a new line
				fp.write(str(S1[i])+' '+str(e_s1[i])+' '+str(S2[i])+' '+str(e_s2[i])+'\n')
		
		e_s1=np.array(e_s1)/np.array(S1)
		e_s2=np.array(e_s2)/np.array(S2)
		S1=np.log(S1)
		S2=np.log(S2)	
		spear = stats.spearmanr(S1, S2)
		per=stats.pearsonr(S1,S2)
		x_range=np.linspace(np.min(S1),np.max(S1),100)
		fig,ax3=plt.subplots()
		a,b,aerr,berr,covab=BCES.bcesp(S1,e_s1,S2,e_s2,np.zeros_like(S1),10000)
		fitm=np.array([a[3],b[3]])	
		covm=np.array([ (aerr[3]**2,covab[3]), (covab[3],berr[3]**2) ])	
		lcb,ucb,x_range=nmmn.stats.confband(S1,S2,a[3],b[3],conf=0.68)
		ax3.errorbar(S1,S2,xerr=e_s1,yerr=e_s2,linewidth=0,elinewidth=1,capsize=1,color='dimgrey',marker='h')
		
		ax3.plot(x_range,pow(x_range,a[3],b[3]),linewidth=2,color='dodgerblue',label=r' $k$='+str(round(a[3],1))+r'$\pm$'+str(round(aerr[3],1))+r' $A$='+str(round(b[3],1))+r'$\pm$'+str(round(berr[3],1))+'\n Spearman: '+str(round(spear[0],2))+'\n Pearson: '+str(round(per[0],2)))
		plt.fill_between(x_range, lcb, ucb, alpha=0.3, facecolor='blue')
		ax3.legend(frameon=False)
		ax3.tick_params(which='major', width=1.00, length=5,right=True,top=True)

		ax3.set_xlabel(r'Log (sum$_{IMG1}$/arcsec$^2$)')
		ax3.set_ylabel(r'Log (sum$_{IMG2}$/arcsec$^2$)')
		fig.savefig('out.jpg')

	if event.key in ['C', 'c'] and toggle_selector.RS.active:
		print('- Grid clearing -')
		ax1.patches.clear()
		ax2.patches.clear()
		plt.draw()
		celle.clear()
		thresh_map=thresh_mapper(thresh,image_data1)

		


plt.style.use(astropy_mpl_style)
plt.rcParams['axes.grid'] = False

f1 = fits.open(image_file1)
head1=f1[0].header

if head1['NAXIS']>2:
	drop2axes(image_file1,image_file1.replace('.fits','_2ax_temp.fits'))
	f1 = fits.open(image_file1.replace('.fits','_2ax_temp.fits'))
	w1 = WCS(f1[0].header)
	image_data1= fits.getdata(image_file1.replace('.fits','_2ax_temp.fits'), ext=0)
else:
	f1 = fits.open(image_file1)
	w1 = WCS(f1[0].header)
	image_data1= fits.getdata(image_file1, ext=0)

f2 = fits.open(image_file2)
head2=f2[0].header
if head2['NAXIS']>2:
	drop2axes(image_file2,image_file2.replace('.fits','_2ax_temp.fits'))
	f2 = fits.open(image_file2.replace('.fits','_2ax_temp.fits'))
	w2 = WCS(f2[0].header)
	image_data2= fits.getdata(image_file2.replace('.fits','_2ax_temp.fits'), ext=0)
else:
	f2 = fits.open(image_file2)
	w2 = WCS(f2[0].header)
	image_data2= fits.getdata(image_file2, ext=0)


pixel_scale_1 = utils.proj_plane_pixel_scales(w1)
pixel_scale_2 = utils.proj_plane_pixel_scales(w2)

scale1=pixel_scale_1[0]*3600.  
scale2=pixel_scale_2[0]*3600.  

if sm_g!='none':
	image_data2 = ndimage.gaussian_filter(image_data2, sigma=(sm_g*scale2, sm_g*scale2), order=0)

thresh_map=thresh_mapper(thresh,image_data1)
click = [None,None]
release = [None,None]

celle=[]
len_celle=0
fig, (ax1, ax2) = plt.subplots(1, 2,figsize=(15, 8))



px_fov, py_fov = w2.wcs_world2pix(w1.wcs.crval[0], w1.wcs.crval[1],1)
w_fov=len(image_data1[:,0])*scale1/scale2/2.
h_fov=len(image_data1[0,:])*scale1/scale2/2.


ax1.imshow(image_data1,origin='lower',aspect='equal',cmap=newcmp,norm=SymLogNorm(vmin=np.nanmean(image_data1)*0.005,vmax=np.nanmean(image_data1)*500.,linthresh=np.nanmean(image_data1)))#vmin=thresh,vmax=10.*thresh
ax1.contour(thresh_map,levels=[thresh],colors='silver',linewidths=0.7)
ax2.imshow(image_data2,origin='lower',cmap=newcmp,norm=SymLogNorm(vmin=np.nanmean(image_data2)*0.005,vmax=np.nanmean(image_data2)*500.,linthresh=np.nanmean(image_data2)))#,
ax2.set_xlim(px_fov-w_fov,px_fov+w_fov)
ax2.set_ylim(py_fov-h_fov,py_fov+h_fov)

ax1.set_title('IM1: '+image_file1.replace('.fits',''))
ax2.set_title('IM2: '+image_file2.replace('.fits',''))

ax1.tick_params(which='major', width=1.00, length=5,right=True,top=True)
ax2.tick_params(which='major', width=1.00, length=5,right=True,top=True)
ax1.set_xlabel('Pixel coord.')
ax1.set_ylabel('Pixel coord.')
ax2.set_xlabel('Pixel coord.')
ax2.set_ylabel('Pixel coord.')

print('---------------------')
print('Interactive commands:')
print("click  -->  release: Define region of interest (ROI)")
print('W: Create rectangular grid in ROI')
print('H: Create hexagonal grid in ROI')
print('D: Mask map in ROI')
print('I: Run PtP analysis with active grids. Output: Plot in out.jpg, out_plot.png, and data series in out.dat')
print('X: Recenter images')
print('+/-: Increase/decrease cell size by 0.5 arcsec')
print('C: Clear grid')

print('---------------------')
celle=[]
# drawtype is 'box' or 'line' or 'none'
toggle_selector.RS = RectangleSelector(ax1, line_select_callback,
									    useblit=True,
									   button=[1,3],  # don't use middle button
									   minspanx=5, minspany=5,
									   spancoords='pixels',
									   interactive=True,
									   props=dict(linestyle='-', color='dodgerblue', 
                                       fill=True, alpha=.4,linewidth=2))
plt.connect('key_press_event', toggle_selector)
plt.tight_layout()

plt.show()
