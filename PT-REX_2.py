#!usr/bin/python
# Filename:PT-REX.py

###########################
"""
Author: Alessandro Ignesti
Ver 2.0

Point-to-point TRend EXtractor [PT-REX]

DATA PREPARATION:
The input maps have to be in CASA image format and they have to finish with the proper code:
-radio map: [name]-radio.fits
-X-rays total instensity map: [name]-img.fits
-X-rays background map: [name]-bmap.fits
-X-rays total exposure map: [name]-emap.fits

The region and mask box files must be in CASA region format (.crtf).

NOTES:
-The X-ray maps regridding procedure (1) is slow and quite useless, it may
 be used only to produce images for publications.
----NEW UPDATE PYTHON 3 CASA 6
      __
      /oo\
     |    |
 ^^  (vvvv)   ^^
 \\  /\__/\  //
  \\/      \//
   /        \
  |          |    ^
  /          \___/ |
 (            )     |
  \----------/     /
    //    \\_____/
   W       W


"""
###########################
import sys
sys.path.append('casa6/lib64/python3.6/site-packages')
sys.path.append('casa6/lib/python3.6/site-packages')

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from astropy.io import fits
import numpy as np
import os

import scipy
import scipy.stats
from scipy.optimize import curve_fit
import bces.bces
import linmix
import nmmn.stats
import gc
from casatasks import *
from casatools import *
import readline
readline.parse_and_bind("tab: complete")

casalog.filter('SEVERE')
os.system('mv *.log casa6/logs')
		

def fit_func(x,a,b):
	return b*x**a

def fit_alg(key,x,ex,y,ey):#k,kerr,A,Aerr,
	
	
	if key=='LS':
		i_g=[np.min(x),1.0]
		z,pcov=curve_fit(fit_func,x,y,sigma=ey,p0=i_g,maxfev=1000)
		k,A=z[0],z[1]
		perr=np.sqrt(np.diag(pcov))
		kerr,Aerr=perr[0],perr[1]
		
	if key=='BCES':
		x_norm,y_norm=x/np.mean(x),y/np.mean(y)
		a,b,aerr,berr,covab=bces.bces.bces(np.log10(x_norm),ex/x_norm,np.log10(y_norm),ey/y_norm,np.zeros(len(y)))
		b[2]=np.mean(y)/(np.mean(x)**a[2])
		k,kerr,A,Aerr=a[2],aerr[2],b[2],berr[2]

	if key=='Ort':
		x_norm,y_norm=x/np.mean(x),y/np.mean(y)
		a,b,aerr,berr,covab=bces.bces.bces(np.log10(x_norm),ex/x_norm,np.log10(y_norm),ey/y_norm,np.zeros(len(y)))
		b[3]=np.mean(y)/(np.mean(x)**a[3])
		k,kerr,A,Aerr=a[3],aerr[3],b[3],berr[3]

	if key=='LinMix':
		x_norm,y_norm=x/np.mean(x),y/np.mean(y)
		lm = linmix.LinMix(np.log10(x), np.log10(y), xsig=np.array(ex)/np.array(x), ysig=np.array(ey)/np.array(y),K=2 )#, delta=delta, K=K, nchains=nchains)xycov=xycov
		lm.run_mcmc(silent=True)#, silent=silent)miniter=10, maxiter=100
		k,kerr=lm.chain['beta'].mean(), lm.chain['beta'].std()
		Am,Aerrm=lm.chain['alpha'].mean(), lm.chain['alpha'].std()
		A=np.mean(y)/(np.mean(x)**k)
		Aerr=A*Aerrm/Am
	return k,kerr,A,Aerr



def start():
	img=[]
	bmap=[]
	emap=[]

	for file in os.listdir("./"):
		if file.endswith("-radio.fits"):
			map_radio=file
		if file.endswith("-img.fits"):
			img.append(file)
		if file.endswith("-bmap.fits"):
			bmap.append(file)
		if file.endswith("-emap.fits"):
			emap.append(file)
		
	if len(bmap)==0:
		
		radio= fits.open(map_radio)
		data=radio[0].data
		header=radio[0].header
		bmap=data*0
		fits.writeto('dummy_b.fits',bmap,header,overwrite=True)
		bmap=['dummy_b.fits']
		print ('No background map was found. X-ray background is set to 0.')
	if len(emap)==0:
		
		radio= fits.open(map_radio)
		header=radio[0].header
		data=radio[0].data
		emap=data/data
		fits.writeto('dummy_e.fits',emap,header,overwrite=True)
		emap=['dummy_e.fits']
		print ('No exposure map was found. X-ray exposure is set to 1.')
	img.sort()
	bmap.sort()
	emap.sort()

	head=imhead(imagename=map_radio)

	if head['restoringbeam']['major']['unit']=='arcsec':
		bmaj,bmin,scala=head['restoringbeam']['major']['value'],head['restoringbeam']['minor']['value'],head['incr'][1]*206264.8
		beam_area=3.14*bmaj*bmin/4.0/0.693147181/scala**2
		beam_area_arcsec=3.14*bmaj*bmin/4.0/0.693147181
	if head['restoringbeam']['major']['unit']=='deg':
		bmaj,bmin,scala=head['restoringbeam']['major']['value']*3600,head['restoringbeam']['minor']['value']*3600,head['incr'][1]*206264.8
		beam_area=3.14*bmaj*bmin/4.0/0.693147181/scala**2
		beam_area_arcsec=3.14*bmaj*bmin/4.0/0.693147181
	camp=np.sqrt(bmaj*bmin)/scala
	print ('-----------------------------------------')
	print ('Radio map: ',map_radio)
	print ('X-rays map(s):',img)
	print ('X-rays background map(s):',bmap)
	print ('X-rays exposure map(s):',emap)
	print ('Beam radio: ',bmaj,'X',bmin,' arcsec')
	print ('Beam area [px]: ',beam_area)
	print ('Beam area [arcsec]: ',beam_area_arcsec)
	print ('1 px= ',scala,' arcsec')
	print ('Beam sampling factor: ',camp)
	print ('------------------------------------------')
	cal_e=np.float32(input('Calibration error of the radio image [es. 0.05]: '))
	stat_fit=input('Fitting method [LS | BCES | Ort | LinMix]: ')

	return(img,bmap,emap,map_radio,bmaj,bmin,scala,beam_area,beam_area_arcsec,camp,cal_e,stat_fit)


logo=mpimg.imread('casa6/logo2.png')
imgplot = plt.imshow(logo)
plt.axis('off')
plt.show(block=False)

print ('\n---Point-to-point TRend EXtractor---')
img=[]
bmap=[]
emap=[]
img,bmap,emap,map_radio,bmaj,bmin,scala,beam_area,beam_area_arcsec,camp,cal_e,stat_fit=start()



extra=', linewidth=1, linestyle=-, color=magenta\n'

while True:
	print ('-------------------------')
	print ('TASK LIST')
	print ('Utilities: ')
	print ('-Regrid X-ray maps [1]')
	print ('-Create a new mask.image from a mask.reg [2]')
	print ('-Create a new J2000 mesh over the radio map [3]')
	print ('Point-to-point analysis:')
	print ('-Single mesh analysis [4]')
	print ('-Monte-Carlo analysis [5]')
	print ('-Spectral index analysis (cooming soon) [6]')
	print ('--Viewer [v]')
	print ('--Reload images and configuration [r]')
	print ('--quit [q]')
	print ('   \n')
	task=input('Task to execute: ')

	if task=='1':
		for k in range(0,len(img)):
			imsmooth(imagename=img[k],beam=head['restoringbeam'],outfile=img[k]+'_smoothed',overwrite=True)
			imsmooth(imagename=bmap[k],beam=head['restoringbeam'],outfile=bmap[k]+'_smoothed',overwrite=True)
			imsmooth(imagename=emap[k],beam=head['restoringbeam'],outfile=emap[k]+'_smoothed',overwrite=True)
			imregrid(imagename=img[k]+'_smoothed',output=img[k]+'_regrid',template=map_radio,overwrite=True)
			imregrid(imagename=bmap[k]+'_smoothed',output=bmap[k]+'_regrid',template=map_radio,overwrite=True)
			imregrid(imagename=emap[k]+'_smoothed',output=emap[k]+'_regrid',template=map_radio,overwrite=True)

		print ('New _regrid maps created')

	if task=='2':
		prjct=input('New mask name: ')
		msk_r=input('Mask regions (.crtf): ')
		os.system('rm -r *_temp.image')
		importfits(fitsimage=map_radio,imagename='radio_temp.image',overwrite=True)
		immath(imagename='radio_temp.image',region=msk_r,outfile='mask_temp.image',expr='IM0/IM0')
		imregrid(imagename='mask_temp.image',template='radio_temp.image',decimate=10,interpolation='nearest',output=prjct+'.image',overwrite=True)
		
		exportfits(imagename=prjct+'.image',fitsimage=prjct+'.fits',overwrite=True)
		os.system('rm -r *_temp.image')
		print ('New mask created: ',prjct,'.fits')


	if task=='3':
		flusso=0.0
		prjct=input('New mesh name: ')
		region=input('Working region [.crtf]: ')
		reg=imstat(imagename=map_radio,region=region)
		xtr,ytr,xbl,ybl=reg['trc'][0],reg['trc'][1],reg['blc'][0],reg['blc'][1]
		msk=input('Mask: ')
		deltax=np.float32(input('Box x-size [pixel]: '))
		deltay=np.float32(input('Box y-size [pixel]: '))
		thresh=np.float32(input('Radio flux density threshold [Jy/beam]: '))
		if beam_area>deltax*deltay:
			choice=(input('WARNING: the box si smaller than the beam. Do you want to choose a new box size? (y/n) '))
			if choice=='y':
				print ('Smallest allowed size: ',np.sqrt(beam_area))
				deltax=np.float32(input('Box x-size [pixel]: '))
				deltay=np.float32(input('box y-size [pixel]: '))
			if choice=='n':
				print ('Proceeding anyway')
		area_arcsec=deltax*deltay*scala**2

		g = open(prjct+'.crtf',"w")

		x=np.arange(xbl,xtr,deltax)
		y=np.arange(ybl,ytr,deltay)
		g.write('#CRTFv0 CASA Region Text Format version 0\n')
		print ('Raw mesh: ',len(x),'x',len(y))
		print ('Total box: ',len(x)*len(y))
		print ('Box size: ',deltax,'x',deltay,'pix= ',deltax*deltay,' pix^2= ',area_arcsec,'arcsec^2')
		print ('-Adapting mesh-')
		i=0

		for k in range(0,len(x)-1):

			for j in range(0,len(y)-1):
				
				box='box [['+str(x[k+1])+'pix,'+str(y[j+1])+'pix], ['+str(x[k])+'pix, '+str(y[j])+'pix]]'
				try:
					mystat=imstat(imagename=map_radio,region=box,listit=False,verbose=False)
				#flusso=np.float(mystat['flux'])

					if np.float32(mystat['flux'])/area_arcsec>thresh/beam_area_arcsec:
							mystat_m=imstat(imagename=msk,region=box,listit=False,verbose=False)['sum']
							if not mystat_m>0.0:

								g.write('box [['+mystat['trcf'][0:12]+', '+mystat['trcf'][14:27]+'], ['+mystat['blcf'][0:12]+', '+mystat['blcf'][14:27]+']] coord=J2000'+extra)
				except:
					continue

				
				s =  ((np.int(i*20/((len(x)-1)*(len(y)-1)))*'#')+(np.int(20-i*20/((len(x)-1)*(len(y)-1)))*' ')+(' ')+str(np.int(i*100/((len(x)-1)*(len(y)-1)-1)))+(' %'))
				sys.stdout.write("\rAdapting - Status: ["+s+"]" )
				sys.stdout.flush()
				i=i+1

		g.close()
		gc.collect()
		print ('\nNew mesh created: ',prjct+'.crtf')

	if task=='4':
		prjct=input('Project name: ')
		grid=(input('Mesh name: '))
		griglia = [line.rstrip('\n') for line in open(grid)]
		
		rms=np.float32(input('RMS of the radio map [Jy/beam]: '))
		p1=open(prjct+'.dat',"w")
		f_x=[]
		e_x=[]
		f_r=[]
		e_r=[]
		
		top=len(griglia)
		for k in range (1,top):
			
			mystat_R=imstat(imagename=map_radio,region=str(griglia[k]),listit=False,verbose=False)
			area_arcsec=(mystat_R['trc'][0]-mystat_R['blc'][0])*(mystat_R['trc'][1]-mystat_R['blc'][1])*scala**2
			flusso_R=np.float32(mystat_R['flux'])/area_arcsec

			error_r=np.sqrt((cal_e*flusso_R)**2+(rms*np.sqrt(np.float32(mystat_R['npts'])/beam_area)/beam_area_arcsec)**2)


		##################################
			l_img=[]
			l_bmap=[]
			l_emap=[]
			for t in range (0,len(img)):
				mystat_x=imstat(imagename=img[t],region=str(griglia[k]),listit=False,verbose=False)['sum']
				l_img.append(np.float32(mystat_x))
				mystat_x_bkg=imstat(imagename=bmap[t],region=str(griglia[k]),listit=False,verbose=False)['sum']
				l_bmap.append(np.float32(mystat_x_bkg))
				mystat_x_exp=imstat(imagename=emap[t],region=str(griglia[k]),listit=False,verbose=False)['mean']
				l_emap.append(np.float32(mystat_x_exp))

			flusso_x=(sum(l_img)-sum(l_bmap))/sum(l_emap)/area_arcsec
			error_x=(np.sqrt(sum(l_img)+sum(l_bmap)))/sum(l_emap)/area_arcsec
			if flusso_x[0]>0 and flusso_R[0]>0 and error_x[0]<flusso_x[0] and error_r[0]<flusso_R[0]:
				f_x.append(flusso_x[0])
				f_r.append(flusso_R[0])
				e_x.append(error_x[0])
				e_r.append(error_r[0])
			if flusso_x<0 or flusso_R<0 or error_x>flusso_x:
				continue


			################
			p1.write(('{} {} {} {}\n'.format(float(flusso_x),float(error_x),float(flusso_R),float(error_r))))

			
			s = ((np.int(k*20/top)*'#')+((np.int(20-k*20/top))*' ')+(' ')+str(np.int(k*100/(top-1)))+(' %'))
			sys.stdout.write("\rAnalysys - Status: ["+s+"]" )
			sys.stdout.flush()
		#print '\nVIB trovati: ',VIB
		print ('\nData file created: ',prjct,'.dat')
		p1.close()
		pers=scipy.stats.pearsonr(np.log10(f_x),np.log10(f_r))
		sper=scipy.stats.spearmanr(np.log10(f_x),np.log10(f_r))
		fit=fit_alg(stat_fit,f_x,e_x,f_r,e_r)
		print ('Best-fit slope k: ',fit[0],'(',fit[1],')')
		print ('Person coeff: ',str(round(pers[0],2)),' P-val: ',str('{:0.3e}'.format(pers[1])))
		print ('Spearman coeff: ',str(round(sper[0],2)),' P-val: ',str('{:0.3e}'.format(sper[1])))
		x2=np.linspace(0.5*np.min(f_x), 1.5*np.max(f_x), 10)
		fit3=fit_func(x2,fit[0],fit[2])
		
		########PLOT###################
		plt.clf()
		plt.rc('font', size=12)
		plt.rc('axes', labelsize=15)
		plt.rc('legend', fontsize=15)
		k_mock=[fit[0]-fit[1],fit[0]+fit[1]]
		
		A_mock=[np.mean(f_r)/(np.mean(f_x)**k_mock[0]),np.mean(f_r)/(np.mean(f_x)**k_mock[1])]
		fit_low=fit_func(x2,k_mock[0],A_mock[0])
		fit_up=fit_func(x2,k_mock[1],A_mock[1])
		plt.fill_between(x2,fit_low,fit_up,facecolor='blue',alpha=0.3)
		
		plt.errorbar(f_x,f_r,xerr=e_x,yerr=e_r,linewidth=0,elinewidth=1,color='black',capsize=2)
		plt.plot(x2,fit3,color='blue',linewidth=1.0,label='k$_{SM}$='+str(round(fit[0],2))+'$\pm$'+str(round(fit[1],2)))
		plt.xscale("log")
		plt.yscale("log")
		plt.xlim(0.7*np.min(f_x), 1.3*np.max(f_x))
		plt.ylim(0.5*np.min(f_r), 1.5*np.max(f_r))
		plt.ylabel('$I_{R}$ [Jy arcsec$^{-2}$]')
		plt.xlabel('$I_{X}$ [ph cm$^{-2}$ s$^{-1}$ arcsec$^{-2}$]')
		#plt.title('k= '+str(a[2])+' st.dev= '+str(aerr[2])+' corr_p= '+str(pers[0]))
		plt.legend()
		plt.show(block=False)
		plt.savefig('plot_'+prjct+'.pdf')
		gc.collect()
		#gc.collect()
	if task=='5':
		prjct=input('Project name: ')
		niter=input('Number of MC iterations: ')
		niter=int(niter)
		region=input('Working region [.crtf]: ')
		reg=imstat(imagename=map_radio,region=region)
		xtr,ytr,xbl,ybl=reg['trc'][0],reg['trc'][1],reg['blc'][0],reg['blc'][1]
		xc,yc,width,height=(xtr+xbl)/2,(ytr+ybl)/2,xtr-xbl,ytr-ybl
		msk=input('Mask: ')
		deltax=np.float32(input('Box x-size [pixel]: '))
		deltay=np.float32(input('Box y-size [pixel]: '))
		area_arcsec=deltax*deltay*scala**2
		rms=np.float32(input('RMS of the radio map [Jy/beam]: '))
		thresh=np.float32(input('Radio flux density threshold [Jy/beam]: '))
		if beam_area>deltax*deltay:
			choice=(input('WARNING: the box size is smaller than the beam. Do you want to choose a new box size? (y/n) '))
			if choice=='y':
				print ('Smallest allowed size: ',np.sqrt(beam_area))
				deltax=np.float32(input('Box x-size [pixel]: '))
				deltay=np.float32(input('box y-size [pixel]: '))
			if choice=='n':
				print ('Proceeding anyway')
		area_arcsec=deltax*deltay*scala**2
		rand_x=np.random.randint(xc-width/4.0,xc+width/4.0,niter)
		rand_y=np.random.randint(yc-height/4.0,yc+height/4.0,niter)
		#slope=[]
		def fit_func(x,a,b):
			return a*x**b
		i=0
		nbox=[]
		s =  ((np.int(i*20/(niter))*'#')+(np.int(20-i*20/(niter))*' ')+(' ')+str(np.int(i*100/(niter-1)))+(' %'))
		sys.stdout.write("\rMC analysis - Status: ["+s+"]" )
		sys.stdout.flush()
		for i in range(0,niter):
			
			f_x=[]
			e_x=[]
			f_r=[]
			e_r=[]
			box_tot=[]
			cont=0.0
			x=np.arange(rand_x[i]-3*width/4,rand_x[i]+3*width/4,deltax)
			y=np.arange(rand_y[i]-3*height/4,rand_y[i]+3*height/4,deltay)
			
			for k in range(0,len(x)-1):
				for j in range(0,len(y)-1):
					box='box [['+str(x[k+1])+'pix,'+str(y[j+1])+'pix], ['+str(x[k])+'pix, '+str(y[j])+'pix]]'
					box_tot.append(box)
					try:
						mystat_R=imstat(imagename=map_radio,region=box,listit=False,verbose=False)
					#flusso=np.np.float32(mystat_R['flux'])
						if mystat_R['flux']/area_arcsec>thresh/beam_area_arcsec:
							mystat_m=np.float16(imstat(imagename=msk,region=box,listit=False,verbose=False)['sum'])
							if not mystat_m>0:
								cont=cont+1
								box_j2000='box [['+mystat_R['trcf'][0:12]+', '+mystat_R['trcf'][14:27]+'], ['+mystat_R['blcf'][0:12]+', '+mystat_R['blcf'][14:27]+']] coord=J2000'
							#g.write(box_j2000+extra2)
								flusso_R=np.float32(mystat_R['flux'])/area_arcsec

								error_r=np.sqrt((cal_e*flusso_R)**2+(rms*np.sqrt(np.float32(mystat_R['npts'])/beam_area)/beam_area_arcsec)**2)


								l_img=[]
								l_bmap=[]
								l_emap=[]
								for t in range (0,len(img)):
									mystat_x=np.float32(imstat(imagename=img[t],region=box_j2000,listit=False,verbose=False)['sum'])
									l_img.append(mystat_x)
									mystat_x_bkg=np.float32(imstat(imagename=bmap[t],region=box_j2000,listit=False,verbose=False)['sum'])
									l_bmap.append(mystat_x_bkg)
									mystat_x_exp=np.float32(imstat(imagename=emap[t],region=box_j2000,listit=False,verbose=False)['mean'])
									l_emap.append(mystat_x_exp)

								flusso_x=(sum(l_img)-sum(l_bmap))/sum(l_emap)/area_arcsec
								error_x=(np.sqrt(sum(l_img)+sum(l_bmap)))/sum(l_emap)/area_arcsec
								if flusso_x>0 and flusso_R>0 and error_x<flusso_x:
									f_x.append(flusso_x[0])
									f_r.append(flusso_R[0])
									e_x.append(error_x[0])
									e_r.append(error_r[0])
								if flusso_x<0 or flusso_R<0 or error_x>flusso_x:
									continue
							#del mystat_R,mystat_m,mystat_x,mystat_x_bkg,mystat_x_exp
					except:
						continue
			#g.close()
			nbox.append(cont)
			fit=fit_alg(stat_fit,f_x,e_x,f_r,e_r)
			k_boot=np.random.normal(fit[0],fit[1],1) #BOOTSTRAPING del valore da usare nella distribuzione!!!!
			gc.collect()
			with open (prjct+'.dat',"a") as myfile:
				myfile.write(str(float(k_boot))+'\n')
			#del f_x,f_r,e_x,e_r,l_img,l_emap,l_bmap,x,y,box,box_j2000
			
			s =  ((np.int(i*20/(niter))*'#')+(np.int(20-i*20/(niter))*' ')+(' ')+str(np.int(i*100/(niter-1)))+(' %'))
			sys.stdout.write("\rMC analysis - Status: ["+s+"]")
			sys.stdout.flush()
			#print i,z[1]
			#gc.collect()



					#g.close()
		slope=np.genfromtxt(prjct+'.dat')
		print ('\nMean slope: ',np.mean(slope))
		print ('st. dev: ',np.std(slope))
		print ('Boxes: ',np.mean(nbox),' (',np.std(nbox),')')

		########PLOT###################
		#sc3=(input('Pulire plot? (s/n) '))
		#if sc3=='s':
		#	plt.cla()
		plt.clf()
		plt.hist(slope)
		plt.title('$k_{MC}=$'+str(round(np.mean(slope),2))+'$\pm$'+str(round(np.std(slope),2)))
		#plt.text(0.1,0.9,'media: '+str(np.mean(slope))+' st.dev: '+str(np.std(slope)))
		plt.show(block=False)
		plt.savefig('plot_'+prjct+'.pdf')
		#gc.collect()
	
	if task=='6':
		prjct=input('Project name: ')
		mappe=[]
		up,mid,low=100,10,120
		error=[]
		freq=[]
		rms_m=[]
		head1=imhead(imagename=map_radio)
		freq1=head1['refval'][2]
		mappe.append(map_radio)
		freq.append(np.float(freq1)*1e-9)
		error.append(cal_e)
		
		rms=np.float32(input('RMS of the radio map [Jy/beam]: '))
		rms_m.append(rms)
		n=input('How many ohter radio map do you want to combine [1: alpha-X-ray ptp 2: color-color plot]: ')
		n=int(n)
		for i in range(int(n)):
			map=input('Fits file: ')
			rms_t=np.float32(input('RMS of the radio map [Jy/beam]: '))
			er=input('Calibration error: ')
			head_t=imhead(imagename=map)
			freq.append(np.float(head_t['refval'][2])*1e-9)
			mappe.append(map)
			error.append(np.float(er))
			rms_m.append(np.float(rms_t))
		if n==1:
			if freq[0]>freq[1]:
				up,low=0,1
				
				
			else:
				up,low=1,0
				
		if n==2:
			if freq[0]>freq[1] and freq[0]>freq[2]:
				if freq[1]>freq[2]:
					#print(a[0],b[0],c[0])
					up,mid,low=0,1,2
				else:
					#print(a[0],c[0],b[0])
					up,mid,low=0,2,1
			if a[1]<b[1] and a[1]<c[1]:
				if b[1]>c[1]:
					#print(b[0],c[0],a[0])
					up,mid,low=1,2,0
				else:
					#print(c[0],b[0],a[0])
					up,mid,low=2,1,0
			if b[1]>a[1]>c[1]:
				#print(b[0],a[0],c[0])
				up,mid,low=1,0,2
			if c[1]>a[1]>b[1]:
				#print(c[0],a[0],b[0])
				up,mid,low=2,0,1
		
		grid=(input('Mesh name: '))
		griglia = [line.rstrip('\n') for line in open(grid)]

		p1=open(prjct+'.dat',"w")
		alpha1=[]
		e_alpha1=[]
		alpha2=[]
		e_alpha2=[]
		
		f_x=[]
		e_x=[]
		top=len(griglia)
		
		for k in range (1,top):
			fl_R=[]
			e_R=[]
			mystat_R1=imstat(imagename=mappe[0],region=str(griglia[k]),listit=False,verbose=False)
			area_arcsec1=(mystat_R1['trc'][0]-mystat_R1['blc'][0])*(mystat_R1['trc'][1]-mystat_R1['blc'][1])*scala**2
			#flusso_R1=np.float32(mystat_R1['flux'])/area_arcsec1
			#flusso_R1=np.float32(mystat_R1['flux'])
			#error_r1=np.sqrt((cal_e*flusso_R1)**2+(rms*np.sqrt(np.float32(mystat_R1['npts'])/beam_area))**2)
			for i in range(0,len(mappe)):
				mystat_Rt=imstat(imagename=mappe[i],region=str(griglia[k]),listit=False,verbose=False)
				flusso_Rt=np.float(mystat_Rt['flux'])
				error_rt=np.sqrt((error[i]*flusso_Rt)**2+(rms_m[i]*np.sqrt(np.float(mystat_Rt['npts'])/beam_area))**2)
				fl_R.append(flusso_Rt)
				e_R.append(error_rt)

			
			
			
		
			l_img=[]
			l_bmap=[]
			l_emap=[]
			for t in range (0,len(img)):
				mystat_x=imstat(imagename=img[t],region=str(griglia[k]),listit=False,verbose=False)['sum']
				l_img.append(np.float32(mystat_x))
				mystat_x_bkg=imstat(imagename=bmap[t],region=str(griglia[k]),listit=False,verbose=False)['sum']
				l_bmap.append(np.float32(mystat_x_bkg))
				mystat_x_exp=imstat(imagename=emap[t],region=str(griglia[k]),listit=False,verbose=False)['mean']
				l_emap.append(np.float32(mystat_x_exp))

			flusso_x=(sum(l_img)-sum(l_bmap))/sum(l_emap)/area_arcsec1
			error_x=(np.sqrt(sum(l_img)+sum(l_bmap)))/sum(l_emap)/area_arcsec1
			if flusso_x[0]>0 and error_x[0]<flusso_x[0]:# and e_alfio1<alfio1 and e_alfio2<alfio2:
				f_x.append(flusso_x[0])
				
				e_x.append(error_x[0])
				#alpha.append(alfio)
				#e_alpha.append(np.float(e_alfio))
				#p1.write(('{} {} {} {} {} {}\n'.format(float(flusso_x),float(error_x),float(alfio1),float(e_alfio1),float(alfio2),float(e_alfio2))))
			if flusso_x<0 or error_x>flusso_x:# or e_alfio1>alfio1 or e_alfio2>alfio2:
				continue

			if n==1:
				alfio1=np.log10(fl_R[up]/fl_R[low])/np.log10(freq[up]/freq[low])
				e_alfio1=np.sqrt((e_R[up]/fl_R[up])**2+(e_R[low]/fl_R[low])**2)/np.log10(freq[up]/freq[low])
				alpha1.append(alfio1)
				e_alpha1.append(e_alfio1)

				p1.write('{} {} {} {}\n'.format(float(flusso_x),float(error_x),float(alfio1),float(e_alfio1)))
			if n==2:

				alfio1=np.log10(fl_R[up]/fl_R[mid])/np.log10(freq[up]/freq[mid])
				e_alfio1=np.sqrt((e_R[up]/fl_R[up])**2+(e_R[mid]/fl_R[mid])**2)/np.log10(freq[up]/freq[mid])
				alpha1.append(alfio1)
				e_alpha1.append(e_alfio1)
				alfio2=np.log10(fl_R[mid]/fl_R[low])/np.log10(freq[mid]/freq[low])
				e_alfio2=np.sqrt((e_R[low]/fl_R[low])**2+(e_R[mid]/fl_R[mid])**2)/np.log10(freq[mid]/freq[low])
				alpha2.append(alfio2)
				e_alpha2.append(e_alfio2)
				p1.write('{} {} {} {} {} {}\n'.format(float(flusso_x),float(error_x),float(alfio1),float(e_alfio1),float(alfio2),float(e_alfio2)))
				
				
			################
#			p1.write(('{} {} {} {}\n'.format(float(flusso_x),float(error_x),float(alfio))))


			s = ((np.int(k*20/top)*'#')+((np.int(20-k*20/top))*' ')+(' ')+str(np.int(k*100/(top-1)))+(' %'))
			sys.stdout.write("\rAnalysys - Status: ["+s+"]" )
			sys.stdout.flush()
		#print '\nVIB trovati: ',VIB
		print ('\nData file created: '+prjct+'.dat')

		p1.close()
		if n==1:
			pers=scipy.stats.pearsonr((f_x),(alpha1))
			sper=scipy.stats.spearmanr((f_x),(alpha1))
			print ('Person coeff: ',str(round(pers[0],2)),' P-val: ',str('{:0.3e}'.format(pers[1])))
			print ('Spearman coeff: ',str(round(sper[0],2)),' P-val: ',str('{:0.3e}'.format(sper[1])))
		
			plt.clf()
			plt.rc('font', size=12)
			plt.rc('axes', labelsize=15)
			plt.rc('legend', fontsize=15)
			plt.clf()

			plt.errorbar(f_x,alpha1,xerr=e_x,yerr=e_alpha1,linewidth=0,capsize=3,elinewidth=1,color='black')
			
			
			
			plt.ylabel('Spectral index '+str(round(freq[low],2))+'-'+str(round(freq[up],2))+' GHz')
			plt.xlabel('$I_{X}$ [ph cm$^{-2}$ s$^{-1}$ arcsec$^{-2}$]')
		
			plt.legend()
			plt.show(block=False)
			plt.savefig('plot_'+prjct+'.pdf',ddp=200)
		if n==2:
			plt.cla()
			colors = cm.rainbow(f_x)
			x_plot=np.linspace(0.1*min(np.min(alpha2),np.min(alpha1)),10.*max(np.max(alpha2),np.max(alpha1)),10)
			norm = matplotlib.colors.Normalize(vmin=min(f_x), vmax=max(f_x), clip=True)
			mapper = cm.ScalarMappable(norm=norm, cmap='plasma')
			time_color = np.array([(mapper.to_rgba(v)) for v in f_x])

			#loop over each data point to plot
			for x, y, e1,e2, color in zip(alpha2, alpha1, e_alpha2,e_alpha1, time_color):
				if e2>0.45*np.abs(y):

					plt.errorbar(x, y, xerr=e1,uplims=True,yerr=0.3,lw=0,elinewidth=.85, capsize=3, color=color)
				#if e2>0.45*y and e1>0.45*x:
				#	plt.errorbar(x, y, xlolims=True, xerr=0.3,lolims=True,yerr=0.3,lw=0,elinewidth=.85, capsize=3, color=color)

				if e1>0.45*np.abs(x):
					plt.errorbar(x, y, xuplims=True,yerr=e2 ,lw=0,elinewidth=.85, capsize=3, color=color)
			    
				else:
					plt.errorbar(x, y, xerr=e1,yerr=e2 ,lw=0,elinewidth=.85, capsize=3, color=color)


			plt.scatter(alpha2,alpha1,c=f_x,cmap='plasma',alpha=1,linewidth=0.6,marker='h',s=100,zorder=100)
			#plt.scatter(alpha1,alpha2,c=f_x,cmap='inferno',alpha=1,linewidth=0.6,edgecolor='silver',marker='h',s=200,zorder=100)
			#plt.errorbar(alpha1,alpha2,xerr=e_alpha1,yerr=e_alpha2,linewidth=0,elinewidth=1.,capsize=1,color='silver',marker='',zorder=0)
			plt.plot(x_plot,x_plot,"k--",linewidth=.5)
			#f_x.sort()
			plt.colorbar(label='$I_{X}$ [ph cm$^{-2}$ s$^{-1}$ arcsec$^{-2}$]')
			
			plt.xlim(0.3*min(np.min(alpha2),np.min(alpha1)),7.*max(np.max(alpha2),np.max(alpha1)))
			plt.ylim(0.3*min(np.min(alpha2),np.min(alpha1)),7.*max(np.max(alpha2),np.max(alpha1)))
			plt.ylabel('Spectral index '+str(round(freq[mid],2))+'-'+str(round(freq[up],2))+' GHz')
			plt.xlabel('Spectral index '+str(round(freq[low],2))+'-'+str(round(freq[mid],2))+' GHz')#$I_{X}$ [ph cm$^{-2}$ s$^{-1}$ arcsec$^{-2}$]
			plt.show(block=False)
			plt.savefig('plot_'+prjct+'.pdf',ddp=200,bbox_size='tight')
					



	if task=='v':
		#viewer(infile=map_radio)
		os.system('./casa6/casa-release-5.4.1-32.el7/bin/casaviewer &')
	
	if task=='r':
		
		img=[]
		bmap=[]
		emap=[]
		img,bmap,emap,map_radio,bmaj,bmin,scala,beam_area,beam_area_arcsec,camp,cal_e,stat_fit=start()
		
		
	if task=='q':

		
		print ('Good bye!')
		break
