
'''
	This module contains all the functions required by selfcal_LTS.py to perform selfcal.

'''
import sys
import os
import numpy as np
import time
import glob
import pickle,pyfits
#from numpy import nanstd, nanmedian #For CASA 5.1
from scipy.stats import nanstd,nanmedian
import datetime as dt
import numpy.ma as ma
from numpy import matrix
import matplotlib.pyplot as plt

import selfcal_input
reload(selfcal_input)
del selfcal_input
from selfcal_input import *

def apply_acorr(dat, ants1, ants2, idx_ant, elen, tile, uvw, BW, fout):
    """
    Amplitude correction for the beamformer-to-receiver cable length
    differences (""decorrelation"")
    """
    c_speed_of_light = 299792458. # m/s
    tauw_sign = +1.
    iau = np.where(ants1 == ants2)[0] # Pointers into ants where a1 == a2
    uant = ants1[iau]  # Unique antennas in the same order as they appear
    nant = len(uant)                  # Number of unique antennas
    antmax = ants1.max()              # Max antenna numeric *name*
    #
    # Fill in the array of antenna indices idx_ant
    #
    iuant = 0
    for ia in uant:
        idx_ant[ia] = iuant
        iuant = iuant + 1

    #
    # Baseline u,v,w coordinates and lengths
    #
    u = uvw[0,:]
    v = uvw[1,:]
    w = uvw[2,:]
    blen = np.sqrt(u**2 + v**2)

    #
    # Amplitude-corrected data block
    #
    cdat = np.zeros_like(dat)

    ndat = len(dat[0,0,:])  # Number of correlations per a pol and a chan

    #
    # Loop over all the baselines
    #

    fout.write('# a1  tl1     elen1     tau1        a2  tl2     ' \
               'elen2     tau2         blen       w      tauw       acor\n') 
    for idat in xrange(ndat):
        a1 = ants1[idat]                  # First antenna 'name'
        a2 = ants2[idat]                  # Second antenna 'name'
        ia1 = idx_ant[a1]                 # First antenna index
        ia2 = idx_ant[a2]                 # Second antenna index

        #
        # Here we calculate the cable delays for each pair of tiles,
        # tau1 and tau2, from their electrical cable lengths, elen. 
        # Also, from the wave number (in meters), w, we get the
        # geometric delay between the tiles, tauw.
        # The delays are combined in tau:
        #
        tau1 = elen[ia1]/c_speed_of_light
        tau2 = elen[ia2]/c_speed_of_light

        # Let us not take this into account for a time being
        tauw = w[idat]/c_speed_of_light
        
        tau = tau2 - tau1 + tauw_sign*tauw

        #
        # The correction factor, acor, is the result of integration:
        #
        # acor = {1/BW}\int_{-BW/2}^{+BW/2} cos(2\pi f tau) df = sinc(BW tau).
        #
        acor = np.sinc(BW*tau)

        line = '%4d %4d %9.2f  %12.5e %4d %4d %9.2f  %12.5e ' \
               '%8.2f %8.2f %12.5e %9.6f\n' % \
               (a1, tile[ia1], elen[ia1], tau1,  \
                a2, tile[ia2], elen[ia2], tau2,  \
                blen[idat], w[idat], tauw, acor)

        fout.write(line)
        
        #print 'tau2=%g, tau1=%g, tauw=%g, w=%g, tau=%g, acor=%g' % \
        #      (tau2, tau1, tauw, w[idat], tau, acor)
        
        cdat[:,:,idat] = dat[:,:,idat] / acor
     
    return cdat

def calc_flagged_ant(msname): 
	ms.open(msname)
	ms.selectinit(datadescid=0)
	flag	=	ms.getdata("flag",ifraxis=True)["flag"]
	ms.close()
	flagged_ant	=	0
	for ant1 in range(128):
		temp=0
		for ant2 in range(128):
			if (ant1==ant2):
				continue
			min_ant=min(ant1,ant2)
			max_ant=max(ant1,ant2)
		
			index=128*min_ant-(min_ant-1)*min_ant/2+(max_ant-min_ant)
			if flag[0,0,index,0]==True or flag[1,0,index,0]==True:
				temp+=1
		if (temp>100):
			flagged_ant+=1
	print 'Flagged antenna: ',flagged_ant
	del flag
	del max_ant
	del min_ant
	del index
	del ant1
	del ant2
	return flagged_ant
	
def calc_psf(freq):
	lambda1		=	299792458.0/freq
	psf		=	1.42*lambda1/max_baseline*180/np.pi*3600   #### this is in arcsecs
	del lambda1	
	return psf

def check_inellipse(x,y,x0,y0,a,b,R=np.eye(2)):
	x_v=np.matrix([[x-x0],[y-y0]])
	R=np.matrix(R)
	xp_v=R*x_v
	xp_v=np.array(xp_v.transpose()).reshape(2)
	if xp_v[0]**2/a**2 + xp_v[1]**2/b**2<1:
		out=True
	else:	
		out=False
	del x_v
	del R
	del xp_v
	return out

def check_source_true_loc(true_loc_folder,msname,ref_time_str,timerange,cellsize,imsize,STOKES,fitstim,chan,taper):
	'''
	This code takes in the MS name and time chosen as reference (ref_pos_time) to solar position on the sky. Code assumes the MS to have calibrator cal tables applied. This code will apply the phase solutions coresponding to the ref_pos_time MS and cleans. It notes the location of the brightest source. 
	'''
	strt=time.time()
	cal_tabs=glob.glob(ref_time_str)
	
	if len(cal_tabs)==0:
		raise RuntimeError('Calibration tables for reference position time slice don\'t exist')
	applycal(vis=msname,gaintable=cal_tabs,timerange=timerange)
	print 'Applycal done..',cal_tabs[0:2],'..etc'
	junkname='test_loc_img'
	clean(vis=msname,threshold='10Jy',niter=1000,imagename=junkname,timerange=timerange,cell=cellsize,imsize=imsize,weighting='natural',stokes=STOKES,psfmode='hogbom',uvtaper=True,outertaper=[taper+'lambda', taper+'lambda', '0deg'],innertaper=['0deg','0deg','0deg'],robust=0.0,gain=0.05)
	exportfits(imagename=junkname+'.image',fitsimage=true_loc_folder+'/Sun_trueloc_img_'+fitstim+'_'+chan+'_full_'+STOKES+'.fits',overwrite=True)
	print 'True location image written out: ','Sun_trueloc_img_',fitstim+'_',chan,'_full_',STOKES,'.fits'
	st=imstat(junkname+'.image')
	maxpix=st['max'][0]
	maxpos=st['maxpos'][0:2]
	clearcal(msname)
	del st
	### Test for dominant source...Reliability of the max pos and derived centroid shift depends on this.###
	st=imhead(junkname+'.image')
	circ1=st['restoringbeam']['major']['value']/cellsize*1.05
	circ2=st['restoringbeam']['minor']['value']/cellsize*1.05
	posang=st['restoringbeam']['positionangle']['value']*np.pi/180.	
	rad=np.max([circ1,circ2])/2*1.6
	del st
	ia.open(junkname+'.image')
	mat=ia.getchunk([maxpos[0]-rad,maxpos[1]-rad],[maxpos[0]+rad,maxpos[1]+rad])
	mat=mat.reshape(mat.shape[0],mat.shape[1])
	ia.close()
	m=mat.shape[0]
	n=mat.shape[1]
	xb,yb=np.where(mat==np.max(mat))
	xb=xb[0]
	yb=yb[0]
	R=np.array([[np.cos(posang),-np.sin(posang)],[np.sin(posang),np.cos(posang)]])
	for i in range(m):
		for j in range(n):
			if(check_inellipse(i,j,xb,yb,circ1/2,circ2/2,R)):
				mat[i,j]=np.nan
			if(check_inellipse(i,j,xb,yb,circ1*3/4,circ2*3/4,R)==False):
				mat[i,j]=np.nan
	medpix=nanmedian(mat.flatten())
	if medpix<maxpix*0.6:
		reliability=True
	else:
		reliability=False
	os.system('rm -rf '+junkname+'*')
	print 'Source True loc found...'
	print 'Time taken: ',(time.time()-strt)/60.,' min.'
	del mat
	del R
	del circ1
	del circ2
	del posang
	del rad
	del strt
	del cal_tabs
	del medpix
	del maxpix
	del xb
	del yb
	del junkname
	return reliability,maxpos


def choose_scales(cellsize,psf):
	sun_dia	=	16*60 ### in arcsecs
	psf_pix	=	int(psf/cellsize)
	scale=[0,psf_pix,3*psf_pix,int(sun_dia/cellsize)]  ### choosing scale to be [0,psf,3*psf, sun_dia/2]
	del sun_dia
	del psf_pix
	return scale	

def clean_as_necessary(want_deep,msname,final_image_poln,cellsize,imsize,weighting,sigma,maskfile,fits_imagename,fits_modelname,BOX,header_dict,scales,threshold='',selfcal_success=True):
	print 'Final deep clean being done..'
	junkimg='junk_'+msname[:-3]
	if (want_deep==True):
		saturate_clean=False
		dyn_range=0
		clean_iter=0
		while (saturate_clean==False):
			clean(vis=msname,imagename=junkimg,outlierfile="",field="",spw="",selectdata=True,timerange='',uvrange="",antenna="",scan="", observation="",intent="",mode="mfs",resmooth=False,gridmode="widefield",wprojplanes=1,facets=1, cfcache="cfcache.dir",rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="", interpolation="linear", niter=2000,gain=0.05,threshold="0Jy",psfmode="hogbom",imagermode="",ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=scales,negcomponent=-1, smallscalebias=0.6,interactive=False,mask=[maskfile],nchan=-1,start=0,width=1,outframe="",veltype="radio", imsize=imsize,cell=cellsize, phasecenter="",restfreq="",stokes=final_image_poln,weighting=weighting,robust=0.0,uvtaper=False,outertaper=[],innertaper=[], modelimage='',restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False, noise="1.0Jy",npixels=0,npercycle=100, cyclefactor=5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,flatnoise=True,allowchunk=False)

			rms	=	imstat(imagename=junkimg+".image",axes=-1,region=region_file,box=BOX,chans="",stokes="",listit=True,verbose=True,mask=[],stretch=False,logfile="",append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto")["rms"][0]
	
			max_pix	=	imstat(imagename=junkimg+".image",axes=-1,region="",box="",chans="",stokes="",listit=True,verbose=True,mask=[],stretch=False,logfile="",append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto")["max"][0]


			if (max_pix/rms<dyn_range or clean_iter>150000):
				saturate_clean=True
			dyn_range=max_pix/rms


			clean_iter+=2000


		exportfits(imagename=junkimg+".image",fitsimage=fits_imagename,history=False,dropdeg=False,overwrite=True)	
		exportfits(imagename=junkimg+".model",fitsimage=fits_modelname,history=False,dropdeg=False,overwrite=True)	
		print 'Fits image and model written out..'
		os.system("rm -rf "+junkimg+".*")
		del rms
		del max_pix
	else:
		if len(threshold)==0:
			imgtag='_'.join(msname.split('_')[:-1])+'*.image'
			lastimg=sorted(glob.glob(imgtag))[-1]
			rms=imstat(imagename=lastimg,box=BOX,listit=True)['rms'][0]
			threshold=str(sigma*rms)+"Jy"	
		clean(vis=msname,imagename=junkimg,outlierfile="",field="",spw="",selectdata=True,timerange='',uvrange="",antenna="",scan="",observation="", intent="",mode="mfs",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0,aterm=True,psterm=False, mterm=True,wbawp=False,conjbeams=True,epjtable="",interpolation="linear",niter=1000000000,gain=0.05,threshold=threshold,psfmode="hogbom",imagermode="",ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=scales,negcomponent=-1,smallscalebias=0.6, interactive=False,mask=[maskfile],nchan=-1,start=0,width=1,outframe="",veltype="radio",imsize=imsize,cell=cellsize,phasecenter="",restfreq="", stokes=final_image_poln,weighting=weighting,robust=0.0,uvtaper=False,outertaper=[],innertaper=[],modelimage='',restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",npixels=0,npercycle=100,cyclefactor=5, cyclespeedup=-1, nterms=1,reffreq="",chaniter=False, flatnoise=True,allowchunk=False)
		
		exportfits(imagename=junkimg+".image",fitsimage=fits_imagename,history=False,dropdeg=False,overwrite=True)
		if selfcal_success==True:	
			exportfits(imagename=junkimg+".model",fitsimage=fits_modelname,history=False,dropdeg=False,overwrite=True)
		os.system("rm -rf "+junkimg+".*")
		rms_val=imstat(imagename=fits_imagename,box=BOX)['rms'][0]
		max_val=imstat(imagename=fits_imagename)['max'][0]
		dyn_range=max_val/rms_val
		del max_val
		del rms_val
	hdu=pyfits.open(fits_imagename,mode='update')
	head=hdu[0].header
	head.update(key='Creators',value='Surajit Mondal, Atul Mohan and Divya Oberoi')
	params,vals=np.genfromtxt('clean.last',delimiter='=',dtype=str,autostrip=True,skip_header=13,skip_footer=True,usecols=(0,1),unpack=True)
	for key,value in zip(params,vals):
		if value=='""':
			value='NIL'
		head.update(key=key[:8],value=value)
	#head.update(key='Input',value=header_dict)
	hdu.flush()
	hdu.close()
	del head
	del hdu
	del params
	del vals
	del junkimg
	return round(dyn_range,2)

def count_tot_antenna_for_clean(antenna_to_use):
	if '~' in antenna_to_use:
		ends=np.array(antenna_to_use.split('~')).astype(int)
		count=ends[1]-ends[0]+1
	if ',' in antenna_to_use:
		ends=antenna_to_use.split(',')
		count=len(ends)
	del ends
	return count

def do_uvsub(msname):
	print 'Doing uvsub..',msname
	uvsub(vis=msname)
	ms.open(msname)
	ms.selectinit(datadescid=0)
	data=ms.getdata('corrected_data')['corrected_data']
	flag=ms.getdata('flag')['flag']
	ms.close()
	MAD,threshold,median=flagger(msname,data,flag)
	del data
	del flag
	return MAD,threshold,median

def decor(msname,metafits,n_tblk,single_time):	# n_tblk: Size of a read/write block: how many individual times in it
	c_speed_of_light = 299792458. # m/s

	#
	# The sign of geometric delay term
	#
	tauw_sign = +1.
	#tauw_sign = -1.
		
	#
	# Read the cable delays (electrical cable lengths, elen)
	# from the metafits file.
	#
	ms_dir	=	msname
	hl = pyfits.open(metafits)
	#
	# The second HDU is a fits table with 256 (nPol*nTile) rows.
	# The order of the rows is a little strange, but they can be
	# ordered in the same order as casa as follows:
	# Find the CASA order of the tiles and elen:
	#
	h1 = hl[1]
	tile_order = np.argsort(h1.data['Tile'])
	#
	# The lengths are prefixed by "EL_" which indicates that the
	# lengths given (in metres) are electrical lengths.
	# Cut off the 'EL_' prefixes and  reorder the e-lengths:
	#
	tile = [til for til in h1.data['Tile'][tile_order][::2]]
	tile = np.array(tile, dtype=int)
	elen = [float(Len[3:]) for Len in h1.data['Length'][tile_order][::2]]
	elen = np.array(elen, dtype=np.float64)
	hl.close()

	    

	#raise  SystemExit

	sqrt_2 = np.sqrt(2.)
	r2d = 180./np.pi
	before = time.time()
	fout = open('decor_'+msname.split('/')[-1]+'.txt', 'w')

	ms.open(ms_dir, nomodify=False)

	#
	# Find exact number of individual times in MS
		#
	md = ms.metadata()
	
	if single_time:
	    n_times = 1
	else:
	    tims = md.timesforfield(0)
	    n_times = len(tims)
	
	# Convert MS times from metadata to yyyymmdd hh:mm:ss 
	ymds = []
	for tim in tims:
	    q = qa.quantity(tim, 's')
	    ymd = qa.time(q, form='ymd', prec=8) # 8 means 12:34:56.75 format
	    ymds.append(ymd)
	
	
	cwids = md.chanwidths(0)       # Hz, 64 channel widths (40000.Hz = 4 kHz)
	cwid = cwids[0]
	#cwid0 = 40000.                 # Hz, channel width (40000.Hz = 4 kHz)
	IntgTime = md.exposuretime(scan=1,  spwid=0,  polid=0)
	tintg = IntgTime['value']     # Integration time. IntgTime['unit'] = 's' 
	an_samp = cwid*tintg       # = 20000, cplx samples in 0.5s of 40kHz channel
	#an_samp = cwid0*tintg       # = 20000, cplx samples in 0.5s of 40kHz channel
	n_samp = int(an_samp)
	ms.selectinit(datadescid=0)  # Untranslatable CASA dirty swearword

	if single_time:
	    # Select the whole MS dataset
	    rec =  ms.getdata(['data'])
	    dat =  np.copy(rec['data'])
	    ants1 = ms.getdata(['antenna1'])['antenna1']
	    ants2 = ms.getdata(['antenna2'])['antenna2']
	    tims =  ms.getdata(['time'])['time']
	    uvw =   ms.getdata(['uvw'])['uvw']
	    idx_ant = 999999*np.ones(512, dtype=int) # Fill with 999999 where empty
	    
	    cdat = apply_acorr(dat, ants1, ants2, idx_ant, elen, tile, \
	                       uvw, cwid, fout)
	
	    rec['data'][:,:,:] = cdat
	    
	    print '%5d, time = %.2f' % (0, tims[0])  # Only one single time
	    
	    #
	    # Put the corrected data back to the MS database
	    #
	    ms.putdata(rec)
	
	    
	else:
	    n_io = n_times//n_tblk
	    if n_times % n_tblk <> 0: n_io = n_io + 1
	    tblk = np.arange(n_io + 1)*n_tblk
	    tblk[-1] = n_times
	
	    for iblk in xrange(n_io):
	        it0 = tblk[iblk]      # Block start (points at the first)
	        it1 = tblk[iblk+1]    # Block end (points at the next after the last)
	        tim0 = tims[it0]
	        tim1 = tims[it1-1]
	        ms.selectinit(reset=True)          # Reset the cumulative selection!!!
	        ms.select({'time':[tim0, tim1]})
	        rec =  ms.getdata(['data'])
	        blkdat =  np.copy(rec['data'])
	        blkants1 = ms.getdata(['antenna1'])['antenna1']
	        blkants2 = ms.getdata(['antenna2'])['antenna2']
	        blktims =  ms.getdata(['time'])['time']
	        blkuvw =   ms.getdata(['uvw'])['uvw']
	        blkcdat = np.zeros_like(blkdat)   # Corrected data to be written to MS
	
	        #
	        # Start/end indices of same times into blktims[]
	        #
	        ibt = np.where(np.diff(blktims) <> 0)[0] + 1
	        ibt = np.hstack((0, ibt, len(blkants1)))
	
	        idx_ant = 999999*np.ones(512, dtype=int) # Fill with 999999 where empty
	
	        #raise  SystemExit
	        #sys.exit(0)
	
	        #
	        # Apply amplitude corrections to the (large) piece of correlation data 
	        # loaded to RAM. For each single time 'tim' and full set of 
	        # baselines (ants1, ants2), a corresponding data block 'dat' is
	        # extracted 
	        # from 'blkdata', corrected, and stored in 'vdat'.
	        #
	
	        for itim in xrange(it0,it1):
	            tim = tims[itim]
	            ibt0 = ibt[itim-it0]
	            ibt1 = ibt[itim-it0+1]
	            dat = blkdat[:,:,ibt0:ibt1]
	            ants1 = blkants1[ibt0:ibt1]
	            ants2 = blkants2[ibt0:ibt1]
	            uvw = blkuvw[:,ibt0:ibt1]
	
	            cdat = apply_acorr(dat, ants1, ants2, idx_ant, elen, tile, \
	                               uvw, cwid, fout)
	
	            blkcdat[:,:,ibt0:ibt1] = cdat
	            
	            tymd = qa.time(qa.quantity(tim, 's'), form='ymd', prec=8)[0]
	
	            print '%5d, time = %.2f or %s' % (itim, tim, tymd)
	
	        print 'Block %d-%d done.' % (it0,it1-1)
	
	
	        #
	        # Store the corrected correlations in the measurement set
	        #
	        rec['data'][:,:,:] = blkcdat
	
	        #
	        # Put the corrected data back to the MS database
	        #
	        ms.putdata(rec)	
	ms.close()
	fout.close()
	return
		
	
def error_msgs(err_code):
	if err_code==1:
		return "split problem"
	elif err_code==2:
		return "dirty image not produced"
	elif err_code==3:
		return "large number of flagged antennas"
	elif err_code==4:
		return "failed to make image in a selfcal iteration"
	elif err_code==5:
		return "decreasing dynamic range"
	elif err_code==6:
		return "no good solution found"
	elif err_code==7:
		return "maximum iteration exceeded."
	elif err_code==8:
		return "Bad choice of reference time."
	else:
		return "succeeded"

def extend_calibrator_flags(calibrator_cal,num_ant=128,need_both_poln=True):
	bpass,flag, success	= read_bandpass(calibrator_cal,num_ant)
	num_chan		=	np.size(flag[0,:,0,0])

	for ant in range(num_ant):
		pos1	=	np.where(flag[0,:,0,ant]==True)[0]
		pos2	=	np.where(flag[1,:,0,ant]==True)[0]
		if ((len(pos1)>0.6*num_chan or len(pos2)>0.6*num_chan) and need_both_poln==True):
			flag[:,:,0,ant]=True
		elif (len(pos1)>0.6*num_chan):
			flag[0,:,0,ant]=True
		elif (len(pos2)>0.6*num_chan):
			flag[0,:,0,ant]=True

	BCAL		=	tb.open(calibrator_cal,nomodify=False)
	flag1		=	tb.getcol('FLAG') # Extract FLAG column from bandpass table
	dims		= 	np.shape(flag1)			
	flag 		=	np.reshape(flag,dims)
	flag1		=	flag
	tb.putcol('FLAG',flag1)
	tb.flush()
	tb.close()

	del flag1
	del flag
	del BCAL
	del bpass
	del success
	return	

def field_of_view(freq):
	FOV=np.sqrt(610)*150*10**6/freq  #610 deg^2 is the image FoV at 150MHz. So extrapolating this to cntrfq
	return FOV*3600 ### in arcsecs

def find_bpass(msname,spwr,UVRANGE,TIMERANGE,snr,ref_ant,APMODE):
	bandpass(vis=msname,caltable=msname[:-3]+'.cal',field="",spw=spwr,intent="",selectdata=True,timerange=TIMERANGE,uvrange=UVRANGE,antenna="",\
				scan="",observation="",msselect="",solint="inf",combine="scan",refant=str(ref_ant),minblperant=4,minsnr=snr,solnorm=False,bandtype="B",smodel=[],\
				append=False,fillgaps=0,degamp=3,degphase=3,visnorm=False,maskcenter=0,maskedge=5,docallib=False,callib="",gaintable=[],gainfield=[],\
				interp=[],spwmap=[],parang=False)

def find_phase_sol(msname,spwr,UVRANGE,TIMERANGE,snr,ref_ant,APMODE):
	cb.open(msname)
	cb.selectvis(chanmode='none',spw=spwr,uvrange=UVRANGE,time=TIMERANGE)
	PHASE_CAL_TABLE_OUT=msname[:-3]+'.cal'
	cb.setsolve(type='B',t='inf',refant=str(ref_ant),apmode=APMODE,table=PHASE_CAL_TABLE_OUT,append=False,minsnr=snr)
	cb.solve()
	cb.close()

def flag_antennas(msname,antenna_to_use):
	if '~' in antenna_to_use:
		return
	antenna_used=np.array(antenna_to_use.split(',')).astype(int)
	print 'Antenna used in clean: ',antenna_used
	flagged_antennas=','.join(np.setdiff1d(np.arange(1,129,1),antenna_used,assume_unique=True).astype(str))
	print 'Antenna to be Flagged: ',flagged_antennas
	print "MS name: ",msname
	if flagged_antennas!='':
		flagdata(vis=msname,mode='manual',antenna=flagged_antennas)
	del flagged_antennas
	del antenna_used
	return

def flagger(msname,data,flag):
	mask_data=ma.array(np.abs(data),mask=flag)	
	median=ma.median(mask_data)
	MAD=ma.median(np.abs(mask_data-median))
	sigma=MAD*1.4826
	flagdata(vis=msname,mode='clip',clipminmax=[0,median+6.5*sigma],clipoutside=True,datacolumn='corrected',flagbackup=False)
	del mask_data
	return MAD,median+6.5*sigma,median

def gen_time_array(starttime,ref_timeslice_done,times,ref_time,tdt,endtime):   #### This function generates the indices of the times
									#### which are to be imaged and also index of the 
									#### reference timeslice.
	yrstr=str(starttime.date())
	wholerng=[]
	stt=starttime
	DSx=[]
	i=0

	while stt<=endtime:
		tempt=stt.time().strftime('%H:%M:%S.%f')[:-5]
		stt+=tdt
		wholerng+=[tempt]
		DSx+=[i]
		i+=1

	time_DSx=dict(zip(wholerng,DSx))
	if ref_timeslice_done==False:
		times=[ref_time]+times
	ref_time_index=time_DSx[ref_time]
	timeranges=[]
	for i in times:	
		if '~' not in i:
			timeranges+=[i]
			continue
		ts=i.split('~')
		strt=ts[0]
		endt=ts[1]
		stt=dt.datetime.strptime(yrstr+' '+strt,'%Y-%m-%d %H:%M:%S.%f')
		endtime=dt.datetime.strptime(yrstr+' '+endt,'%Y-%m-%d %H:%M:%S.%f')
		while stt<=endtime:
			tmp=stt.time().strftime('%H:%M:%S.%f')[:-5]
			ntmp=stt+tdt
			timeranges+=[tmp]  # Change the image averaging time delta here
			stt+=tdt
		
	time_beg=[time_DSx[i] for i in timeranges]
	return time_beg,ref_time_index


def generate_new_bookkeeping_array():
	max_pix_array=np.zeros(20)*np.nan
	dyn_range_array=np.zeros(20)*np.nan
	return max_pix_array,dyn_range_array

def getnearpos(array,value):
	### returns index of two elements nearest to a given number #####

	a = abs(array-value)
    	b=np.argsort(a)
	del a
    	return b[0],b[1]

def model_phase_rms_curve(phase_rms,iteration,time):
	### here I assume that the rms curve is similar for both polarisations.
	### will model rms for XX only ####

	num_cal_done=len(phase_rms)
	start_new_cal=num_cal_done+1
	time_temp=scan_timeranges(time)	
	filename="time_"+"_".join(time_temp.split("~")[0].split(":"))+"_self_"+str(start_new_cal)+".cal"	
	while (os.path.isdir(filename)):
		bpass,flag,success=read_bandpass(filename,128)
		pos		=	np.where(flag==True)[0]
		bpass[pos]	=	np.nan
		polx_data	=	np.angle(bpass[0,0,0,:])
		phase_rms.append(nanstd(polx_data))
		iteration.append(start_new_cal)
		start_new_cal+=1
		filename="time_"+"_".join(time_temp.split("~")[0].split(":"))+"_self_"+str(start_new_cal)+".cal"
	temp=np.log10(np.array(phase_rms))	
	iter_temp=iteration[-np.size(iteration)/3:-5]
	temp=temp[-np.size(temp)/3:-5]
	if (np.size(temp)>10):
		line=np.polyfit(iter_temp,temp,deg=1)	
		return line[0]*180/np.pi
	else:
		return 0	

def num_pixels(FOV,cellsize):
	num	=	FOV/cellsize
	pow2	=	int(np.log2(num))
	possibility=	np.array([2**(pow2-1)*3,2**(pow2-2)*5,2**(pow2-2)*7,2**(pow2+1)])
	return possibility[getnearpos(possibility,num)[0]]

def pixelsize(freq):
	psf	=	calc_psf(freq)	
	pixel	=	int(psf/3)
	return pixel

def read_bandpass(CALTABLE,NANT):	
	BCAL=tb.open(CALTABLE,nomodify=True) # Open the bandpass table 
	bpass=tb.getcol('CPARAM') # Extract the CPARAM=GAIN column from the bandpass table
	flag=tb.getcol('FLAG') # Extract FLAG column from bandpass table
	snr=tb.getcol('SNR');
	dims = bpass.shape	
	NPOL = dims[0]
	NCHAN = dims[1]
	NTIME = dims[2]/NANT			
	flag = np.reshape(flag,[NPOL,NCHAN,NTIME,NANT])
	bpass=	np.reshape(bpass,[NPOL,NCHAN,NTIME,NANT])
	success=True
	tb.close()
	del BCAL
	del snr
	del NPOL
	del NCHAN
	del NTIME
	return bpass,flag, success

def scan_timeranges(casa_times):
	frac_sec	=	casa_times-int(casa_times)
	time_string	=	time.gmtime(casa_times)   #### Remember gmtime ignores fractional seconds.
	time1		=	str(time_string[3])+":"+str(time_string[4])+":"
	
	if (frac_sec<0.5): #### MWA times are in .25 or in .75 seconds. 
		time1	=	time1+str(time_string[5])+"~"+time1+str(time_string[5])+".5"
	else:
		time1	=	time1+str(time_string[5])+".5"+"~"+time1+str(time_string[5]+1)
	del frac_sec
	del time_string
	return time1

def sigma_in_unknown(slope):
	#### calculates by how much sigma should go down in unknown ######		
	if slope==0:
		return 0
	elif slope>-0.17:
		return 0
	else:
		return 0.2

def stop_sigma(max_pix,intelligence_file):
	#### intelligence array has first element as flux and second element as sigma

	intelligence_array=pickle.load(open(intelligence_file,"rb"))
	intelligence_array=np.array(intelligence_array)
	idx1,idx2=getnearpos(intelligence_array[:,0],max_pix)
	intell_pix1=intelligence_array[idx1,0]
	intell_pix2=intelligence_array[idx2,0]
	intell_sig1=np.log10(intelligence_array[idx1,1])
	intell_sig2=np.log10(intelligence_array[idx2,1])
	sigma=intell_sig2+(intell_sig2-intell_sig1)/(intell_pix2-intell_pix1)*(max_pix-intell_pix2)
	del intelligence_array
	del idx1
	del idx2
	del intell_pix1
	del intell_pix2
	del intell_sig1
	del intell_sig2
	return max(10**sigma,1.5)

def update_bookkeeping_arrays(max_pix,dyn_range,max_pix_array,dyn_range_array):
	pos=np.where(np.isnan(dyn_range_array)==True)[0]
	if np.size(pos)!=0:  #### that means that the book-keeping arrays are not filled
		dyn_range_array[pos[0]]=dyn_range  ### numpy where stores in a sequential order.
						   ### hence using first element for updating
		max_pix_array[pos[0]]=max_pix
	else:
		dyn_range_array[:-1]=dyn_range_array[1:]
		max_pix_array[:-1]=max_pix_array[1:]
		dyn_range_array[-1]=dyn_range
		max_pix_array[-1]=max_pix
	print '$$$$$$$$$$$$$$$$$$$$$$'
	print 'Updated Dyn Range array: ',dyn_range_array[0:5]
	print '$$$$$$$$$$$$$$$$$$$$$$$'
	return max_pix_array,dyn_range_array

####################################################### Main Selfcal Functions #########################################################

def skip_selfcal(channel_msname,part_time,header_dict,chan,final_image_poln,modeldir='',basedir='',imagedir='',corr="XX,YY",region_file='',BOX='',cell=0,imsize=[0],scales=[-1],maskfile='',phase_snr=3,refant='1',selfcal_failed=False,continuous=True,):
	print "entered skip selfcal\n"
	log_file=open('CASA_imaging.txt','a')
	print '#### scan_timeranges ' 
	gen_str=scan_timeranges(part_time).split('~')[0]
	#channel=channel_msname.split('~')[-1].split('.ms')[0]
	channel=channel_msname.split('.ms')[0].split('_')[-1]
	chan=channel
	ftim=scan_timeranges(part_time).split('~')[0].split(':')
	if len(ftim[-1])==1:
		ftim[-1]='0'+ftim[-1]
		ftim[-1]=ftim[-1]+'.0'
	elif len(ftim[-1])==2 and '.' not in ftim[-1]:
		ftim[-1]+='.0'
	elif len(ftim[-1])==3 and '.' in ftim[-1]:
		ftim[-1]='0'+ftim[-1]
	if len(ftim[0])!=2:
		ftim[0]='0'+ftim[0]
	if len(ftim[1])!=2:
		ftim[1]='0'+ftim[1]
	fitstim=''.join(ftim)
	time_to_image=fitstim.split('~')[0]

	
	time_to_image_object=dt.datetime(1,1,1,int(time_to_image[0:2]),int(time_to_image[2:4]),int(time_to_image[4:6]),int(time_to_image[7])*100000)  #### assume that the day is same
	
	final_model_name=basedir+'/'+modeldir+'/'+'Sun_Clean_Matrix_'+fitstim+'_'+chan+'_full_'+final_image_poln+'_model.fits'
	final_image_name=basedir+'/'+imagedir+'/'+'Sun_Clean_Matrix_'+fitstim+'_'+chan+'_full_'+final_image_poln+'.fits'
	final_uvfits='Sun_Clean_Matrix_'+fitstim+'_'+chan+'_full_'+final_image_poln+'.uvfits'

	relevant_files=sorted(glob.glob(basedir+'/'+modeldir+'/'+'Sun_Clean_Matrix_*'+'_'+chan+'_full_'+final_image_poln+'_model.fits'))
	
	temp=[i.split('Sun_Clean_Matrix_')[-1].split('_')[0] for i in relevant_files]

	times_done_object=[dt.datetime(1,1,1,int(i[0:2]),int(i[2:4]),int(i[4:6]),int(i[7])*100000) for i in temp]
	time_diff=[i-time_to_image_object for i in times_done_object]
	del times_done_object
	del temp
	absolute_time_diff=abs(np.array([i.days*86400+i.seconds+i.microseconds*1e-6 for i in time_diff]))
	del time_diff
	rel_pos1=np.argsort(absolute_time_diff)[0]
	time_diff=absolute_time_diff[rel_pos1]
	if time_diff>max_time_delta and selfcal_failed==False:
		log_file.write(gen_str+" time difference large\n")
		log_file.write(gen_str+" going for a selfcal\n")
		log_file.close()
		return 9

	#rel_pos1,rel_pos2=getnearpos(times_done,float(fitstim))
	#print rel_pos1,rel_pos2
	
	nearest_selfcal_model=relevant_files[rel_pos1]
	nearest_time=nearest_selfcal_model.split('/')[-1].split('_')[3]
	nearest_caltable='caltables/caltable_'+nearest_time+'.cal'
	nearest_selfcal_image=basedir+'/'+imagedir+'/'+nearest_selfcal_model.split('/')[-1].split('_model.fits')[0]+".fits"
	
	#nearest_time_imaged=nearest_selfcal_image.split('_')[3]
	

	filename="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))
	os.system("rm -rf "+filename+"_self_*")  #### this is to delete any folder which has a name same as this
	
	ms_out	=	filename+"_self_1.ms"
	
	num_split_trials=0
	ms_size=0
	ms_file_made=False
	while (ms_size==0 and ms_file_made==False and num_split_trials<5):
		if (os.path.isdir(ms_out)):
			os.system("rm -rf "+ms_out)	
		split(vis=channel_msname,outputvis=ms_out,keepmms=True,field="",spw="",scan="",antenna="",correlation=corr,timerange=scan_timeranges(part_time),intent="",array="",uvrange="",observation="",feed="",datacolumn="data",keepflags=True,width=1,timebin="0s",combine="")
		ms_file_made=os.path.isdir(ms_out)
		if ms_file_made==False:
			num_split_trials+=1
			continue		
		ms_size=os.stat(ms_out+"/table.f21_TSM0").st_size
		num_split_trials+=1
	if (num_split_trials==5):
		#### if initial split did not work not point in continuing #####
		log_file.write(gen_str+" split problem\n")
		log_file.close()
		return 1
	
	log_file.write(gen_str+" Applying solutions from "+nearest_selfcal_model+"\n")
	msname=ms_out
	

	try:
		os.mkdir('caltables')
	except OSError:
		pass
	
	if nearest_time[-1]=='0':
		timerange=nearest_time[0:2]+":"+nearest_time[2:4]+":"+nearest_time[4:]+"~"+nearest_time[0:2]+":"+nearest_time[2:4]+":"+nearest_time[4:6]+".5"
	else:
		timerange=nearest_time[0:2]+":"+nearest_time[2:4]+":"+nearest_time[4:]+"~"+nearest_time[0:2]+":"+nearest_time[2:4]+":"+nearest_time[4:6]+".9999"

	if os.path.isdir(nearest_caltable)==False:
		importfits(fitsimage=nearest_selfcal_model,imagename=nearest_selfcal_model[:-5]+".model")
		ft(vis=channel_msname,model=nearest_selfcal_model[:-5]+".model",usescratch=True)
		bandpass(vis=channel_msname, caltable=nearest_caltable, refant=str(refant), minsnr=phase_snr, timerange=timerange)
		delmod(vis=channel_msname,scr=True)
		os.system("rm -rf "+nearest_selfcal_model[:-5]+".model")

	applycal(vis=msname,gaintable=nearest_caltable,interp=["nearest","nearest"],timerange=scan_timeranges(part_time))
	
	rms=imstat(imagename=nearest_selfcal_image,region=region_file,box=BOX)["rms"][0]
	log_file.write(gen_str+" rms of nearby selfcal image: "+str(rms)+"\n")
	dyn_range=clean_as_necessary(want_deep,msname,final_image_poln,cell,imsize,weighting,deep_sigma,maskfile,final_image_name,final_model_name,BOX,header_dict,scales,threshold=str((deep_sigma+increase_threshold_no_selfcal)*rms)+"Jy",selfcal_success=False)	

	if (dyn_range<worst_DR_allowed and selfcal_failed==False):
		os.system("rm -rf "+final_image_name)
		log_file.write(gen_str+" DR very bad: "+str(dyn_range)+"\n")
		log_file.write(gen_str+" going for a selfcal\n")
		log_file.close()
		return 9
		
	exportuvfits(vis=msname,fitsfile=final_uvfits,datacolumn='corrected',overwrite=True)
	os.system("rm -rf "+msname[:-2]+"ms*")
	DRfil=open('Dynamic_Ranges.txt','a')
	DRfil.write(fitstim+'\t'+str(dyn_range)+'\n')
	DRfil.close()
	log_file.write(gen_str+" DR good: "+str(dyn_range)+"\n")
	log_file.close()
	return 0	
	

def do_selfcal(channel_msname,part_time,scratch=True,corr="XX,YY",ref_ant=1,max_iter=1000,num_ant=128,region_file='',BOX='',cell=0,imsize=[0],scales=[-1],maskfile='', start_sigma=10,phase_snr=3,beyond_intelligence=False,DR_delta1=20,DR_delta2=20,calc_final_sigma=True,intelligence_file='', limit_dyn=1e6, dyn_range_array=np.array([np.nan]), max_pix_array=np.array([np.nan]),bookkeeping_file='',worst_DR_allowed=100,safety_standard=0,final_sigma=10,need_uvsub=True,want_check_dyn=False,wavelength=[1.0]):

### channel_msname is the channel ms
#### returns True when selfcal ends normally
#### code donesn't go below 6 sigma by default
####dyn_range_array and max_pix_array should be supplied in second run.
####Otherwise it will assume that it is the first run. 
#### scratch=True means that you want to do full selfcal. During this set ref_time_index to -1.
### anyway it will not be used if scratch=True. Otherwise it should be time_index for which you want to apply the 
#### solutions.
#### bookkeeping file is the file in which we have stored the dyn_range_array and max_pix_array. if this is a blank string
####then I will assume that you have not done any imaging before for this dataset. Then I will generate a nan filled array for both
#### of these parameters and then populate them iteratively. But if you have done at least one imaging and want to use them for other
#### timeslices if possible, go ahead and give the dyn_range array. The order will be first the actual numbers and then nan. Don't put
#### nan in between numbers.

#### I am just checking here if the array is completely filled with nan or are there some actual numbers
	log_file=open('CASA_imaging.txt','a')
	gen_str=scan_timeranges(part_time).split('~')[0]
	global in_unknown
	sigma=start_sigma
	if np.isnan(dyn_range_array[0]):
		if bookkeeping_file=='':
			max_pix_array,dyn_range_array=generate_new_bookkeeping_array()
		else:
			bookarray=np.loadtxt(bookkeeping_file)
			max_pix_array=bookarray[:,0]
			dyn_range_array=bookarray[:,1]
			pos=np.where(max_pix_array<0)[0]
			max_pix_array[pos]=np.nan
			dyn_range_array[pos]=np.nan
			del pos
			del bookarray
	print '####\nScratch to begin: ',scratch,'\n#######'
	if (scratch==True):
		end_ant		=	50
	else:
		end_ant		=	128

	data_to_split='data'
	
#### take care that ref_time_index is non-negative when scratch=False
	if (scratch==False):
		available_caltables	=	glob.glob("time*_self_1.cal")
		junk_times		=	np.array([float(i.split('_')[1])*60+float(i.split('_')[2])*60+float(i.split('_')[3]) for i in available_caltables])	
		current_time_array	=	np.array(scan_timeranges(part_time).split('~')[0].split(':')).astype(float)
		time_difference		=	np.abs(junk_times-(current_time_array[0]*60+current_time_array[1]*60+current_time_array[2]))
		nearest_ref_time	=	np.where(time_difference==np.min(time_difference))[0]	
		file_time		=	'_'.join(available_caltables[nearest_ref_time[0]].split('_')[:-1])	
		data_to_split		=	"corrected"
		num_caltables=len(glob.glob(file_time+"*.cal"))
		num_iteration=num_caltables
		phasetable		=	[]
		timeslice=scan_timeranges(part_time)
		for  table_num in range(1,num_iteration):
			phasetable.append(file_time+"_"+str(table_num)+".cal")   #### the last caltable made is time*_self_(num_iteration).cal"
		#### All caltables except the last one will be applied in calonly mode #####
		#### The last one will be calflag. All we can apply these at one go on the whole dataset ####
		print 'File time: ',file_time	
		applycal(vis=channel_msname,field="",spw="",intent="",selectdata=True,timerange=timeslice,uvrange="",antenna="",scan="",observation="",msselect="",docallib=False, callib="",gaintable=phasetable,gainfield=[],interp='nearest',spwmap=[],calwt=[True],parang=False,applymode="calonly",flagbackup=False)
		
		phasetable=[file_time+"_"+str(num_iteration)+".cal"]
		print 'Phasetable: ',phasetable
		log_file.write(gen_str+' Phasetable: '+str(phasetable)+'\n')
		applycal(vis=channel_msname,field="",spw="",intent="",selectdata=True,timerange=timeslice,uvrange="",antenna="",scan="",observation="",msselect="",docallib=False, callib="",gaintable=phasetable,gainfield=[],interp='nearest',spwmap=[],calwt=[True],parang=False,applymode="calflag",flagbackup=False)

	filename="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))
	os.system("rm -rf "+filename+"_self_*")  #### this is to delete any folder which has a name same as this
	
	ms_out	=	filename+"_self_1.ms"
	
	num_split_trials=0
	ms_size=0
	ms_file_made=False
	while (ms_size==0 and ms_file_made==False and num_split_trials<5):
		if (os.path.isdir(ms_out)):
			os.system("rm -rf "+ms_out)	
		split(vis=channel_msname,outputvis=ms_out,keepmms=True,field="",spw="",scan="",antenna="",correlation=corr,timerange=scan_timeranges(part_time),intent="",array="",uvrange="",observation="",feed="",datacolumn=data_to_split,keepflags=True,width=1,timebin="0s",combine="")
		ms_file_made=os.path.isdir(ms_out)
		if ms_file_made==False:
			num_split_trials+=1
			continue		
		ms_size=os.stat(ms_out+"/table.f21_TSM0").st_size
		num_split_trials+=1
	if (num_split_trials==5):
		#### if initial split did not work not point in continuing #####
		log_file.close()
		del gen_str
		return 1,False,max_pix_array,dyn_range_array
	
	
	msname	=	ms_out
	del ms_out
	
	antenna_to_use	=	'1~'+str(end_ant)	
	do_bandpass	=	False	
	taper		=	str(0.95*max_baseline/wavelength[0])	

	clean(vis=msname,imagename=msname[:-3],outlierfile="",field="",spw="",selectdata=True,timerange=scan_timeranges(part_time),uvrange="",antenna=antenna_to_use,scan="",observation="",intent="",mode="mfs",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="",interpolation="linear",niter=0,gain=0.05,threshold="2mJy",psfmode="hogbom",imagermode="",ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1,smallscalebias=0.6,interactive=False,mask=[maskfile],nchan=-1,start=0,width=1,outframe="",veltype="radio",imsize=imsize,cell=cell,phasecenter="",restfreq="",stokes="I",weighting="natural",robust=0.0,uvtaper=True,outertaper=[taper+'lambda', taper+'lambda', '0deg'],innertaper=['0deg','0deg','0deg'],modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",npixels=0,npercycle=100,cyclefactor=5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,flatnoise=True,allowchunk=False)

	
	if os.path.isdir(msname[:-3]+".image"):
		rms	=	imstat(imagename=msname[:-3]+".image",axes=-1,region=region_file,box=BOX,chans="",stokes="",listit=True,verbose=True,mask=[],stretch=False,logfile="",append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto")["rms"][0]
	else:
	#### if initial clean did not work not point in continuing #####
		log_file.close()
		del gen_str
		return 2,False,max_pix_array,dyn_range_array
		

	max_pix	=	imstat(imagename=msname[:-3]+".image",axes=-1,region="",box='',chans="",stokes="",listit=True,verbose=True,mask=[],stretch=False,logfile="",append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto")["max"][0]

##### Here I will check for a match with the maximum pixel in this image and those stored in the max_pix_array.
#### Unless it is empty, we can always use it. Of course it is better once we populate it fully. So first I see 
#### if the first element is nan. If yes than all are nan. This also means nearest_dyn_range is defined only if
#### there is at least one element in the max_pix_array
	do_full_selfcal=False
	if np.isnan(max_pix_array[0])==False:
		ratio=max_pix_array/max_pix
		pos1=np.where(ratio>0.1)[0]
		pos2=np.where(ratio<10)[0]
		pos=np.intersect1d(pos1,pos2)
		print 'Ratio array: ',ratio
		log_file.write(gen_str+' Ratio array: '+str(ratio)+'\n')
		del ratio
		del pos1
		del pos2	
		if np.size(pos)!=0:
			pos1=0
			if (np.size(pos)>=2):
				pos1,pos2=getnearpos(max_pix_array[pos],max_pix)
				max_dyn_range=max(dyn_range_array[pos1],dyn_range_array[pos2])
				pos3=np.where(dyn_range_array==max_dyn_range)[0]
				nearest_dyn_range=dyn_range_array[pos3]*max_pix/max_pix_array[pos3] #### here I store the nearest dyn_range which we will try to match
				del pos1,pos2,pos3	
			else:
				nearest_dyn_range=dyn_range_array[pos1]*max_pix/max_pix_array[pos1]
				del pos1
		else:
		##### If this is the case I will do a full selfcal run #######
			do_full_selfcal=True
			nearest_dyn_range=100000000.0
		del pos
	else:
		do_full_selfcal=True
		nearest_dyn_range=10000000000.0
	print 'Do full selfcal: ',do_full_selfcal
	if do_full_selfcal==True and scratch==False:
		scratch=True
		end_ant=50
		antenna_to_use='1~'+str(end_ant)
		os.system(" rm -rf "+msname[:-3]+"*")
		data_to_split="data"
		split(vis=channel_msname,outputvis=msname,keepmms=True,field="",spw="",scan="",antenna="",correlation=corr,timerange=scan_timeranges(part_time),intent="",array="",uvrange="",observation="",feed="",datacolumn=data_to_split,keepflags=True,width=1,timebin="0s",combine="")	
		# Making dirty map...
		print 'Making dirty map for finding Max Pix and RMS.'
		log_file.write(gen_str+' Making dirty map for finding Max Pix and RMS.\n')
		clean(vis=msname,imagename=msname[:-3],outlierfile="",field="",spw="",selectdata=True,timerange=scan_timeranges(part_time),uvrange="",antenna=antenna_to_use,scan="",observation="",intent="",mode="mfs",resmooth=False,gridmode="", wprojplanes=-1, facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True, wbawp=False,conjbeams=True,epjtable="",interpolation="linear",niter=0,gain=0.05,threshold="2mJy", psfmode="hogbom",imagermode="",ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=[],negcomponent=-1, smallscalebias=0.6,interactive=False,mask=[maskfile],nchan=-1,start=0,width=1, outframe="",veltype="radio",imsize=imsize,cell=cell,phasecenter="",restfreq="",stokes="I",weighting="natural",robust=0.0,uvtaper=True,outertaper=[taper+'lambda', taper+'lambda', '0deg'],innertaper=['0deg','0deg','0deg'],modelimage="",restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",npixels=0,npercycle=100,cyclefactor=5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,flatnoise=True,allowchunk=False)	
		
		if os.path.isdir(msname[:-3]+".image"):
			rms	=	imstat(imagename=msname[:-3]+".image",axes=-1,region=region_file,box=BOX,chans="",stokes="",listit=True,verbose=True,mask=[],stretch=False,logfile="",append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto")["rms"][0]
		else:
		#### if initial clean did not work not point in continuing #####
			log_file.close()
			del gen_str
			return 2,False,max_pix_array,dyn_range_array

		max_pix	=	imstat(imagename=msname[:-3]+".image",axes=-1,region="",box='',chans="",stokes="",listit=True,verbose=True,mask=[],stretch=False,logfile="",append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto")["max"][0]
			
	
	del do_full_selfcal	
	print '############\nFinal scratch: ',scratch,'\n##########\n'
	log_file.write(gen_str+' Final scratch: '+str(scratch)+"\n")
	if scratch==True:
		want_uvsub=need_uvsub
	else:
		want_uvsub=False
	log_file.write(gen_str+' UVSUB needed: '+str(want_uvsub)+"\n")
	

	dyn_range_1	=	max_pix/rms

	dyn_range_2	=	dyn_range_1


	end_selfcal		=	False
	num_iteration	=	0
	num_iteration_fixed_ant=0
	flag_num		=	0
	dummy_1			=	'junk1'
	dummy_2			=	'junk2'

	increase_ant_flag	=	False #### This flag is used to check if I have addded more antennas in cleaning.
					#### This is because, it is expected to reduce the dynamic range. So we need
					##### to put different premature stopping criteria on dynamic range for two cases:-
					#### when we add antenna and when we do not add.
	if (calc_final_sigma==True):
		if scratch==False:
		#### Since you have already applied the initial calibrator solutions
		#### you will not be off by a very large amount. Hence your rms in 
		#### dirty map will not be as high. Hence use a larger stop threshold
			final_sigma=min(stop_sigma(max_pix,intelligence_file)+0.5,10)
			print 'Calc_final_sigma was set to true.. final sigma: ',final_sigma
			print 'Scratch: ',scratch
			log_file.write(gen_str+' Calc_final_sigma was set to true.. final sigma: '+str(final_sigma)+'\n'+'Scratch: '+str(scratch)+'\n')
		else:
			final_sigma=stop_sigma(max_pix,intelligence_file)
			print 'Scratch: ',scratch
			print 'Final Sigma calculated: ',final_sigma
			log_file.write(gen_str+' Final Sigma calculated: '+str(final_sigma)+'\nScratch'+str(scratch)+'\n')		
		if (scratch==False):
			#### Here I specify the start sigma if you have already done one selfcal round.
			#### Since you have already some good solutions, you can do a deep selfcal in
			#### the first round only.
			sigma=max(final_sigma,6)
			print 'Scratch :',scratch
			print 'Sigma set: ',sigma
			log_file.write(gen_str+' Scratch : '+str(scratch)+'\nSigma set: '+str(sigma)+'\n')
	elif threshold_check :
		final_sigma_default=min(stop_sigma(max_pix,intelligence_file)+0.5,10)
		final_sigma=max(final_sigma_default,6,final_sigma)
	if beyond_intelligence==True:
		niter_in_unknown=0

	modelimage=''
	num_iter_fixed_sigma=0  #### will check how many selfcal rounds it went through in the same sigma
	uvsub_counter=0
	while (end_selfcal==False):
		print 'Inside while loop of selfcal..'
		log_file.write(gen_str+' Inside while loop of selfcal..\n')
		num_ant_clean=count_tot_antenna_for_clean(antenna_to_use)
		if num_ant_clean<12:
			if scratch==True:
				file_tag="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*"
				os.system('rm -rf '+dummy_1)
				os.system('rm -rf '+dummy_2)
				os.system('rm -rf '+msname)
				log_file.close()
				del gen_str
				return 3,False,max_pix_array,dyn_range_array
			else:
				log_file.write(gen_str+" going for a selfcal with scratch=True because "+error_msgs(3)+"\n")
				print "\n$$$$$$$$$$$$$$$$$$\ngoing for a selfcal with scratch=True because "+error_msgs(3)+"\n$$$$$$$$$$$$$$$$$$$$$$$$\n"
				log_file.close()
				del gen_str
				msg,success,max_pix_array,dyn_range_array=do_selfcal(channel_msname,part_time,scratch=True,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=start_sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn, dyn_range_array=dyn_range_array,max_pix_array=max_pix_array,bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn)
				return msg,success,max_pix_array,dyn_range_array

		clean(vis=msname,imagename=msname[:-3],outlierfile="",field="",spw="",selectdata=True,timerange='',uvrange="",antenna=antenna_to_use,scan="",observation="",intent="",mode="mfs",resmooth=False,gridmode="",wprojplanes=-1,facets=1,cfcache="cfcache.dir",rotpainc=5.0,painc=360.0,aterm=True,psterm=False,mterm=True,wbawp=False,conjbeams=True,epjtable="",interpolation="linear",niter=1000000,gain=0.05
,threshold=str(sigma*rms)+"Jy",psfmode="hogbom",imagermode="",ftmachine="mosaic",mosweight=False,scaletype="SAULT",multiscale=scales,negcomponent=-1,smallscalebias=0.6,interactive=False,mask=[maskfile],nchan=-1,start=0,width=1,outframe="",veltype="radio",imsize=imsize,cell=cell,phasecenter="",restfreq="",stokes="I",weighting="natural",robust=0.0,uvtaper=True,outertaper=[taper+'lambda', taper+'lambda', '0deg'],innertaper=[],modelimage=modelimage,restoringbeam=[''],pbcor=False,minpb=0.2,usescratch=False,noise="1.0Jy",npixels=0,npercycle=100,cyclefactor=5,cyclespeedup=-1,nterms=1,reffreq="",chaniter=False,flatnoise=True,allowchunk=False)
		clearstat()
	
		if os.path.isdir(msname[:-3]+".image")==False:
			if scratch==True:
				file_tag="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*"
				os.system('rm -rf '+dummy_1)
				os.system('rm -rf '+dummy_2)
				flag_antennas(msname,antenna_to_use)
				log_file.close()
				del gen_str
				return 4,False,max_pix_array,dyn_range_array
			else:
				log_file.write(gen_str+" going for a selfcal with scratch=True because "+error_msgs(4)+"\n")
				print "\n$$$$$$$$$$$$$$$$$$\ngoing for a selfcal with scratch=True because "+error_msgs(4)+"\n$$$$$$$$$$$$$$$$$$$$$$$$\n"
				log_file.close()
				del gen_str
				msg,success,max_pix_array,dyn_range_array=do_selfcal(channel_msname,part_time,scratch=True,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=start_sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn, dyn_range_array=dyn_range_array,max_pix_array=max_pix_array,bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn,wavelength=wavelength)
				return msg,success,max_pix_array,dyn_range_array
		try:	
			temp_rms	=	imstat(imagename=msname[:-3]+".image",axes=-1,region=region_file,box=BOX,chans="",stokes="",listit=True,verbose=True,mask=[],stretch=False,logfile="",append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto")["rms"][0]
		except TypeError:
			if scratch==True:
				file_tag="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*"
				os.system('rm -rf '+dummy_1)
				os.system('rm -rf '+dummy_2)
				flag_antennas(msname,antenna_to_use)
				log_file.close()
				del gen_str
				return 4,False,max_pix_array,dyn_range_array
			else:
				log_file.write(gen_str+" going for a selfcal with scratch=True because "+error_msgs(4)+"\n")
				print "\n$$$$$$$$$$$$$$$$$$\ngoing for a selfcal with scratch=True because "+error_msgs(4)+"\n$$$$$$$$$$$$$$$$$$$$$$$$\n"
				log_file.close()
				del gen_str
				msg,success,max_pix_array,dyn_range_array=do_selfcal(channel_msname,part_time,scratch=True,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=start_sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn, dyn_range_array=dyn_range_array,max_pix_array=max_pix_array,bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn,wavelength=wavelength)
				return msg,success,max_pix_array,dyn_range_array	
		
			
		max_pix	=	imstat(imagename=msname[:-3]+".image",axes=-1,region="",box="",chans="",stokes="",listit=True,verbose=True,mask=[],stretch=False,logfile="",append=True,algorithm="classic",fence=-1,center="mean",lside=True,zscore=-1,maxiter=-1,clmethod="auto")["max"][0]	
		
		dyn_range	=	max_pix/temp_rms
			
	#### This is a check if the rms is diverging. then no point in continuing self-cal. I have allowed the selfcal iterationto continue for few iterations before doing this check
	#### Actually Huib Intenma suggested this. He said that often there is some jerkiness when starting selfcal. hence this was done.
		print '###########\n'
		print dyn_range_1, dyn_range_2, dyn_range,increase_ant_flag,num_iteration
		print '###########\n'
		log_file.write(gen_str+'  Dynamic range'+str(dyn_range_1)+'   '+str(dyn_range_2)+'   '+str(dyn_range)+'   '+str(increase_ant_flag)+'   '+str(num_iteration)+'\n')
		if ((dyn_range<0.85*dyn_range_1 and dyn_range<0.9*dyn_range_2 and dyn_range_1>dyn_range_2) and increase_ant_flag==False and end_ant>120):
			model_to_delete="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self_"+str(num_iteration+1)  #### for deleting the previous model
			os.system("rm -rf "+model_to_delete+".model")
			os.system("rm -rf "+model_to_delete+".image")
			if uvsub_counter<1 and want_uvsub==True:
				MAD,THRESHOLD,MEDIAN=do_uvsub(dummy_1)
				log_file.write(gen_str+ ' MAD : '+str(round(MAD,4))+'  MEDIAN: '+str(round(MEDIAN,4))+'  THRESHOLD: '+str(round(THRESHOLD,1))) 
				os.system('rm -rf '+dummy_2[:-3]+'*')
				os.system('rm -rf '+msname)
				num_iteration-=2
				split(vis=dummy_1,outputvis=dummy_2,datacolumn='data')
				os.system('rm -rf '+dummy_1[:-3]+'*')
				os.system('mv '+dummy_2+' '+dummy_1)
				msname=dummy_1
				uvsub_counter+=1
				dyn_range_2=dyn_range_1
				dummy_2='junk2'
				dummy_1='junk1'
				continue
			
			if scratch==True:
				print '##########\n DR 2nd last: ',dyn_range_1,'\nDR 1st last: ',dyn_range_2,'\n##############'
				log_file.write(gen_str+' DR 2nd last: '+str(dyn_range_1)+'\nDR 1st last: '+str(dyn_range_2)+'\n')
				file_tag="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*"
				if "time" in dummy_1 and "time" in dummy_2 and "time" in msname:	
					flag_antennas(dummy_1,antenna_to_use)	
					os.system('rm -rf '+dummy_2)
					os.system('rm -rf '+msname)
				success=False
				if dyn_range_1>worst_DR_allowed:
					success=True
					os.system('rm -rf '+file_tag+str(num_iteration)+'.cal '+file_tag+str(num_iteration-1)+".cal")	
				log_file.close()
				del gen_str	
				return 5,success,max_pix_array,dyn_range_array
			else:
				log_file.write(gen_str+" going for a selfcal with scratch=True because "+error_msgs(5)+"\n")
				print "\n$$$$$$$$$$$$$$$$$$\ngoing for a selfcal with scratch=True because "+error_msgs(5)+"\n$$$$$$$$$$$$$$$$$$$$$$$$\n"
				log_file.close()
				del gen_str
				msg,success,max_pix_array,dyn_range_array=do_selfcal(channel_msname,part_time,scratch=True,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=start_sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn, dyn_range_array=dyn_range_array,max_pix_array=max_pix_array,bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn,wavelength=wavelength)
				return msg,success,max_pix_array,dyn_range_array	
		elif ((dyn_range<0.8*dyn_range_1 and dyn_range<0.85*dyn_range_2 and dyn_range_1>dyn_range_2) and increase_ant_flag==True and end_ant>120):
			model_to_delete="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self_"+str(num_iteration+1)  #### for deleting the previous model
			os.system("rm -rf "+model_to_delete+".model")
			os.system("rm -rf "+model_to_delete+".image")
			if uvsub_counter<1 and want_uvsub==True:
				MAD,THRESHOLD,MEDIAN=do_uvsub(dummy_1)
				log_file.write(gen_str+ ' MAD : '+str(round(MAD,4))+'  MEDIAN: '+str(round(MEDIAN,4))+'  THRESHOLD: '+str(round(THRESHOLD,1))) 
				os.system('rm -rf '+dummy_2[:-3]+'*')
				os.system('rm -rf '+msname)
				num_iteration-=2
				split(vis=dummy_1,outputvis=dummy_2,datacolumn='data')
				os.system('rm -rf '+dummy_1[:-3]+'*')
				os.system('mv '+dummy_2+' '+dummy_1)
				msname=dummy_1
				uvsub_counter+=1
				dyn_range_2=dyn_range_1
				dummy_2='junk2'
				dummy_1='junk1'
				continue	
			if scratch==True:
				print '##########\n DR 2nd last: ',dyn_range_1,'\nDR 1st last: ',dyn_range_2,'\n##############'
				log_file.write(gen_str+' DR 2nd last: '+str(dyn_range_1)+'\nDR 1st last: '+str(dyn_range_2)+'\n')
				file_tag="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*"
				if "time" in dummy_1 and "time" in dummy_2 and "time" in msname:	
					flag_antennas(dummy_1,antenna_to_use)	
					os.system('rm -rf '+dummy_2)
					os.system('rm -rf '+msname)
				success=False
				if dyn_range_1>worst_DR_allowed:
					success=True
					os.system('rm -rf '+file_tag+str(num_iteration)+'.cal '+file_tag+str(num_iteration-1)+".cal")
				log_file.close()
				del gen_str		
				return 5,success,max_pix_array,dyn_range_array
			else:
				log_file.write(gen_str+" going for a selfcal with scratch=True because "+error_msgs(5)+"\n")
				print "\n$$$$$$$$$$$$$$$$$$\ngoing for a selfcal with scratch=True because "+error_msgs(5)+"\n$$$$$$$$$$$$$$$$$$$$$$$$\n"
				log_file.close()
				del gen_str
				msg,success,max_pix_array,dyn_range_array=do_selfcal(channel_msname,part_time,scratch=True,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=start_sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn, dyn_range_array=dyn_range_array,max_pix_array=max_pix_array,bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn,wavelength=wavelength)	
				return msg,success,max_pix_array,dyn_range_array

		elif ((dyn_range<0.9*dyn_range_2 and dyn_range_2>1.5*dyn_range_1) and increase_ant_flag==False and end_ant>120):
			model_to_delete="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self_"+str(num_iteration)  #### for deleting the previous model
			os.system("rm -rf "+model_to_delete+".model")
			os.system("rm -rf "+model_to_delete+".image")	
			if uvsub_counter<1 and want_uvsub==True:
				MAD,THRESHOLD,MEDIAN=do_uvsub(dummy_1)
				log_file.write(gen_str+ ' MAD : '+str(round(MAD,4))+'  MEDIAN: '+str(round(MEDIAN,4))+'  THRESHOLD: '+str(round(THRESHOLD,1))) 
				os.system('rm -rf '+dummy_2[:-3]+'*')
				os.system('rm -rf '+msname)
				num_iteration-=2
				split(vis=dummy_1,outputvis=dummy_2,datacolumn='data')
				os.system('rm -rf '+dummy_1[:-3]+'*')
				os.system('mv '+dummy_2+' '+dummy_1)
				msname=dummy_1
				uvsub_counter+=1
				dyn_range_2=dyn_range_1
				dummy_2='junk2'
				dummy_1='junk1'
				continue
			if scratch==True:
				print '##########\n DR 2nd last: ',dyn_range_1,'\nDR 1st last: ',dyn_range_2,'\n##############'
				log_file.write(gen_str+' DR 2nd last: '+str(dyn_range_1)+'\nDR 1st last: '+str(dyn_range_2)+'\n')
				file_tag="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*"
				if "time" in dummy_1 and "time" in dummy_2 and "time" in msname:	
					flag_antennas(dummy_1,antenna_to_use)
					os.system('rm -rf '+dummy_2)
					os.system('rm -rf '+msname)
				success=False	
				if dyn_range_1>worst_DR_allowed:
					success=True
					os.system('rm -rf '+file_tag+str(num_iteration)+'.cal '+file_tag+str(num_iteration-1)+".cal")	
				log_file.close()
				del gen_str	
				return 5,success,max_pix_array,dyn_range_array
			else:
				log_file.write(gen_str+" going for a selfcal with scratch=True because "+error_msgs(5)+"\n")
				print "\n$$$$$$$$$$$$$$$$$$\ngoing for a selfcal with scratch=True because "+error_msgs(5)+"\n$$$$$$$$$$$$$$$$$$$$$$$$\n"
				log_file.close()
				del gen_str
				msg,success,max_pix_array,dyn_range_array=do_selfcal(channel_msname,part_time,scratch=True,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=start_sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn, dyn_range_array=dyn_range_array,max_pix_array=max_pix_array,bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn,wavelength=wavelength)
				return msg,success,max_pix_array,dyn_range_array	

		### Here I give the criteria to stop self-cal or increase number of antennas.
		### Assume that if dynmaic range is does not change by 20%, change antenans or stop.

		increase_ant_flag=False
			
		if (dyn_range>=limit_dyn and end_ant>=120):	
			model_to_delete="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self_"+str(num_iteration)  #### for deleting the previous model
			os.system("rm -rf "+model_to_delete+".model")	
			os.system("rm -rf "+model_to_delete+".image")	
			if scratch==True:
				max_pix_array,dyn_range_array=update_bookkeeping_arrays(max_pix,dyn_range,max_pix_array,dyn_range_array)
				os.system('rm -rf '+dummy_1)
				os.system('rm -rf '+dummy_2)
				flag_antennas(msname,antenna_to_use)
				log_file.close()
				del gen_str
				return 0,True,max_pix_array,dyn_range_array
			else:
				if (want_check_dyn==True):
					if (dyn_range<0.85*nearest_dyn_range):
						log_file.write(gen_str+" going for a selfcal with scratch=True because dynamic range was below expectation. Expected="+str(0.85*nearest_dyn_range)+" Obtained="+str(dyn_range)+"\n")
						print "\n$$$$$$$$$$$$$$$$$$\ngoing for a selfcal with scratch=True because dynamic range was below expectation. Expected="+str(0.85*nearest_dyn_range)+" Obtained="+str(dyn_range)+" \n$$$$$$$$$$$$$$$$$$$$$$$$\n"
						log_file.close()
						del gen_str
						msg,success,max_pix_array,dyn_range_array=do_selfcal(channel_msname,part_time,scratch=True,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=start_sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn, dyn_range_array=dyn_range_array,max_pix_array=max_pix_array,bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn,wavelength=wavelength)	
						return msg,success,max_pix_array,dyn_range_array
					else:
						max_pix_array,dyn_range_array=update_bookkeeping_arrays(max_pix,dyn_range,max_pix_array,dyn_range_array)
						os.system("rm -rf "+"time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*.cal")
						os.system('rm -rf '+dummy_1)
						os.system('rm -rf '+dummy_2)
						flag_antennas(msname,antenna_to_use)
						log_file.close()
						del gen_str
						return 0,True,max_pix_array,dyn_range_array
				else:
					max_pix_array,dyn_range_array=update_bookkeeping_arrays(max_pix,dyn_range,max_pix_array,dyn_range_array)
					os.system("rm -rf "+"time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*.cal")
					os.system('rm -rf '+dummy_1)
					os.system('rm -rf '+dummy_2)
					flag_antennas(msname,antenna_to_use)
					log_file.close()
					del gen_str
					return 0,True,max_pix_array,dyn_range_array	
					
						
		elif (abs(dyn_range-dyn_range_2)<DR_delta1 and abs(dyn_range-dyn_range_1)<DR_delta2 and do_bandpass==True and abs(dyn_range/dyn_range_2-1)<0.08):
			if sigma==final_sigma and uvsub_counter<1 and want_uvsub==True:
				print 'doing uvsub here'
				log_file.write(gen_str+' doing uvsub here\n')
				MAD,THRESHOLD,MEDIAN=do_uvsub(msname)
				log_file.write(gen_str+ ' MAD : '+str(round(MAD,4))+'  MEDIAN: '+str(round(MEDIAN,4))+'  THRESHOLD: '+str(round(THRESHOLD,1))) 
				split(vis=msname,outputvis='junk.ms',datacolumn='data')
				os.system('rm -rf '+msname)
				os.system('mv junk.ms '+msname)
				uvsub_counter+=1
				continue
				
			if (sigma>final_sigma):
				if (safety_standard==0):
					min_num_iter_fixed_sigma=0
					if (scratch==True):
						min_iteration=10
					else:
						min_iteration=5
				elif (safety_standard==1):
					min_num_iter_fixed_sigma=5
					if (scratch==True):
						min_iteration=20
					else:
						min_iteration=10
				else:
					min_num_iter_fixed_sigma=10
					if (scratch==True):
						min_iteration=30
					else:
						min_iteration=20

				if (num_iter_fixed_sigma>min_num_iter_fixed_sigma and num_iteration >min_iteration):
					sigma=max(sigma-0.5,final_sigma)
					print 'Changed sigma to ',sigma
					log_file.write(gen_str+' Changed sigma to '+str(sigma)+'\n')
					num_iter_fixed_sigma=-1   ### will add 1 later in this for loop itself	
			else:
			#### Here the supplied intelligence has been used. Now if the user has given permission then 
			#### the code will try to go cautiously in unknown territory.########
			
				if (beyond_intelligence==False):
					model_to_delete="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self_"+str(num_iteration)  #### for deleting the previous model
					os.system("rm -rf "+model_to_delete+".model")
					os.system("rm -rf "+model_to_delete+".image")
					if scratch==True:
						max_pix_array,dyn_range_array=update_bookkeeping_arrays(max_pix,dyn_range,max_pix_array,dyn_range_array)
						os.system('rm -rf '+dummy_1)
						os.system('rm -rf '+dummy_2)
						flag_antennas(msname,antenna_to_use)
						log_file.close()
						del gen_str
						return 0,True,max_pix_array,dyn_range_array
					else:
						if (want_check_dyn==True):
							if (dyn_range<0.85*nearest_dyn_range):
								log_file.write(gen_str+" going for a selfcal with scratch=True because dynamic range was below expectation. Expected="+str(0.85*nearest_dyn_range)+" Obtained="+str(dyn_range)+"\n")
								print "\n$$$$$$$$$$$$$$$$$$\ngoing for a selfcal with scratch=True because dynamic range was below expectation. Expected="+str(0.85*nearest_dyn_range)+" Obtained="+str(dyn_range)+" \n$$$$$$$$$$$$$$$$$$$$$$$$\n"
								log_file.close()
								del gen_str
								msg,success,max_pix_array,dyn_range_array=do_selfcal(channel_msname,part_time,scratch=True,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=start_sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn, dyn_range_array=dyn_range_array,max_pix_array=max_pix_array,bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn,wavelength=wavelength)	
								return msg,success,max_pix_array,dyn_range_array
							else:
								max_pix_array,dyn_range_array=update_bookkeeping_arrays(max_pix,dyn_range,max_pix_array,dyn_range_array)	
								os.system("rm -rf "+"time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*.cal")
								os.system('rm -rf '+dummy_1)
								os.system('rm -rf '+dummy_2)
								flag_antennas(msname,antenna_to_use)
								log_file.close()
								del gen_str
								return 0,True,max_pix_array,dyn_range_array
						else:
							max_pix_array,dyn_range_array=update_bookkeeping_arrays(max_pix,dyn_range,max_pix_array,dyn_range_array)	
							os.system("rm -rf "+"time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*.cal")
							os.system('rm -rf '+dummy_1)
							os.system('rm -rf '+dummy_2)
							flag_antennas(msname,antenna_to_use)
							log_file.close()
							del gen_str
							return 0,True,max_pix_array,dyn_range_array
				else:
			#### the model_phase_rms_curve also takes phase_rms as a parameter so that we can give the existing list already
			#### that way we save time in it... 
					slope=model_phase_rms_curve([],[],time_num)
					sigma_diff=sigma_in_unknown(slope)
					sigma=sigma-sigma_diff
					in_unknown=True
					niter_in_unknown=0
					
	
		elif (abs(dyn_range/dyn_range_2-1)<0.08 and abs(dyn_range/dyn_range_1-1)<0.08 and num_iteration_fixed_ant>=5):
			if (end_ant<109):	
				end_ant	=	end_ant+20
				num_iteration_fixed_ant=-1
				flag_num	+=	3   #### I am reducing the number of flaggged antennas as I increase the antennas
				increase_ant_flag=True
			elif (end_ant<128):
				end_ant=128
				num_iteration_fixed_ant=-1
				increase_ant_flag=True
			else:
				if (do_bandpass==False):
					do_bandpass=True
	
		'''
		This is the part where I do things which will keep me safe in unknown territory.
		I will not go more than 10 iterations without having a substantial amount of increase in dynamic range	
		'''

		if (in_unknown==True):
			niter_in_unknown+=1
			if niter_in_unknown==1:
				max_dyn_in_unknown=dyn_range	
				max_ms_name=msname
				last_dyn_increase_iter=1    #### This variable is to check how many iterations
							    #### before did I have the last increase in dynamic 
							    #### range
			else:
			##### Here I am saving the maximum dynamic range obtained in unknown domain ######
			##### I am doing this so that I have something to fall back on if all hell breaks loose. ####
				if (max_dyn_in_unknown<dyn_range):	
					max_dyn_in_unknown=dyn_range	
					max_ms_name=msname
					last_dyn_increase_iter=niter_in_unknown
				else:
					if (niter_in_unknown-last_dyn_increase_iter==10):
			##### I will wait for 10 iterations to increase the dynamic range. After that I will just assume that
			#### I am going downhill.
						end_selfcal=True
						a	=	max_ms_name.split('_')
						b	=	str(-1+int(a[-1].split('.')[0]))+".ms"
			#### Here I am using the ms which was produced just before max_ms_name. This is because the clean task
			#### always cleans the corrected data column if given choice. The corrected data column of its previous 
			#### was actually used to produce the image which had the highest dynamic range.

						last_proper_self_iteration=int(a[-1].split('.')[0])
						del a[-1]
						a.append(b)
						msname	=	'_'.join(a)

						#### Here I delete the bad caltables and bad ms files
						filename="_".join(max_ms_name[:-3].split('_')[:-1])+"_"+str(last_proper_self_iteration)
						if (os.path.isdir(filename)):
							os.system("rm -rf "+filename+".cal")
							os.system("rm -rf "+filename+".ms")
							last_proper_self_iteration+=1
							filename="_".join(max_ms_name[:-3].split('_')[:-1])+"_"+str(last_proper_self_iteration)
					
						continue
				
			
		
		
		### put model in ms
		ft(vis=msname,field="",spw="",model=msname[:-3]+".model",nterms=1,reffreq="",complist="",incremental=False,usescratch=True)
	

		
		if dyn_range>dyn_range_1 and dyn_range>dyn_range_2 and end_ant==128 and num_iteration_fixed_ant>=5:
			modelimage=msname[:-3]+".model"
		else:
			modelimage=""


		snr	=	phase_snr
		
		if (do_bandpass==False):
			find_phase_sol(msname,'','','',snr,ref_ant,'p')	
		else:
			find_bpass(msname,'','','',snr,ref_ant,'p')
			
				
	
		### Here I read the caltable
		CALTABLE = msname[:-3]+'.cal'
		if (os.path.isdir(CALTABLE)==False):
			model_to_delete="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self_"+str(num_iteration)  #### for deleting the previous model
			os.system("rm -rf "+model_to_delete+".model")
			os.system("rm -rf "+model_to_delete+".image")
			if scratch==True:
				os.system("rm -rf "+"time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*.cal")
				os.system('rm -rf '+dummy_1)
				os.system('rm -rf '+dummy_2)
				os.system('rm -rf '+msname)
				log_file.close()
				del gen_str
				return 6,False,max_pix_array,dyn_range_array
			else:

				log_file.write(gen_str+" going for a selfcal with scratch=True because "+error_msgs(6)+"\n")
				print "\n$$$$$$$$$$$$$$$$$$\ngoing for a selfcal with scratch=True because "+error_msgs(6)+"\n$$$$$$$$$$$$$$$$$$$$$$$$\n"
				log_file.close()
				del gen_str
				msg,success,max_pix_array,dyn_range_array=do_selfcal(channel_msname,part_time,scratch=True,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=start_sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn, dyn_range_array=dyn_range_array,max_pix_array=max_pix_array,bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn,wavelength=wavelength)
				return msg,success,max_pix_array,dyn_range_array
		bpass,flag, success	= read_bandpass(CALTABLE,num_ant)
		#### if failed to find phase solutions then reduce snr and try again.
		size	=	0		
		pos1=np.where(flag[0,0,0,:]==True)[0]	
		pos2=np.where(flag[1,0,0,:]==True)[0]
		pos=np.union1d(pos1,pos2)
		size	=	np.size(pos)
		
					
				
		del bpass
		del flag
		del success	
		del CALTABLE

			### apply the found phase solutions
		applycal(vis=msname,field="",spw="",intent="",selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",msselect="",docallib=False,callib="",gaintable=msname[:-3]+".cal",gainfield=[],interp='nearest',spwmap=[],calwt=[True],parang=False,applymode="calonly",flagbackup=False)

		##### decided not to do split. So doing a split by hand by puting corrected data
		#### in place of data
		

		a	=	msname.split('_')
		b	=	str(1+int(a[-1].split('.')[0]))+".ms"
		del a[-1]
		a.append(b)
		ms_out	=	'_'.join(a)
		del a
		del b	
		
		os.system('cp -r '+msname+' '+ms_out)
		if (num_iteration >1):
			os.system('rm -rf '+dummy_1)
			dummy_1=dummy_2
			dummy_2=msname
		elif num_iteration==1:
			dummy_2=msname
		else:
			dummy_1=msname 
		tb.open(ms_out)
		corrected_data=tb.getcol("CORRECTED_DATA")
		tb.close()
		tb.open(ms_out,nomodify=False)
		data=tb.getcol("DATA")
		tb.putcol("DATA",corrected_data)
		tb.flush()
		tb.close()
		del data
		del corrected_data
		

	
		image_to_delete="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self_"+str(num_iteration+1) ### for deleting the current image
		model_to_delete="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self_"+str(num_iteration)  #### for deleting the previous model
		
		os.system("rm -rf "+model_to_delete+".image")
		os.system("rm -rf "+image_to_delete+".residual")
		os.system("rm -rf "+image_to_delete+".psf")
		os.system("rm -rf "+model_to_delete+".model")   
		os.system("rm -rf "+image_to_delete+".flux")
		os.system("rm -rf "+image_to_delete+".mask")
		

	

		min_flagged_ant	=	0  #### useful for checking if antenans flagged are included
							   ### in clean list

		if (np.size(pos)!=0):
			pos				=	pos+np.ones(size)
			min_flagged_ant	=	int(min(pos))
		else:
			min_flagged_ant	=	0
		
		del antenna_to_use
		antenna_to_use	=	''

		if (min_flagged_ant< end_ant):
			for i in range(1,end_ant):
				if i not in pos:
					antenna_to_use=antenna_to_use+str(i)+","
			if end_ant not in pos:
				antenna_to_use=antenna_to_use+str(end_ant)
			else:
				antenna_to_use=antenna_to_use[:-1]
		
		else:
			antenna_to_use	=	'1~'+str(end_ant)
		log_file.write(gen_str+' Antennas used for next clean:\n '+antenna_to_use+'\n')

		num_iteration	+=	1
		num_iter_fixed_sigma+=	1
		num_iteration_fixed_ant+=1
		if (num_iteration>max_iter):
			model_to_delete="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self_"+str(num_iteration)  #### for deleting the previous model
			os.system("rm -rf "+model_to_delete+".model")
			os.system("rm -rf "+model_to_delete+".image")
			if (dyn_range>worst_DR_allowed and scratch==True):
				os.system("rm -rf "+"time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*.cal")
				os.system('rm -rf '+dummy_1)
				os.system('rm -rf '+dummy_2)
				flag_antennas(msname,antenna_to_use)
				log_file.close()
				del gen_str
				return 7,True,max_pix_array,dyn_range_array
			else:
				log_file.write(gen_str+" going for a selfcal with scratch=True because "+error_msgs(7)+"\n")
				print "\n$$$$$$$$$$$$$$$$$$\ngoing for a selfcal with scratch=True because "+error_msgs(7)+"\n$$$$$$$$$$$$$$$$$$$$$$$$\n"
				log_file.close()
				del gen_str
				msg,success,max_pix_array,dyn_range_array=do_selfcal(channel_msname,part_time,scratch=True,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=start_sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn, dyn_range_array=dyn_range_array,max_pix_array=max_pix_array,bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn,wavelength=wavelength)
				return msg,success,max_pix_array,dyn_range_array
		
		print '\n######################\n Dyn Range : ',dyn_range,'\n'+'Time: '+scan_timeranges(part_time).split('~')[0]+'\n######################'
		log_file.write(gen_str+' Dyn Range : '+str(dyn_range)+'\n')
		rms		=	temp_rms
		dyn_range_1	=	dyn_range_2	
		dyn_range_2	=	dyn_range	


		
		del msname
		msname	=	ms_out
		del ms_out
	

def control_selfcal(basedir,header_dict,msname,intelligence_file,times,starttime,endtime,tdt,imagedir='',modeldir='',true_loc_imgs_folder='',out_dir='chan_',corr="XX,YY",need_both_poln=True, ref_ant=1,max_iter=1000,num_ant=128,region_file='',cell=-1.0,calc_imsize=True, imsize=[-1],calc_scale_flag=True,scales=[-1],maskfile='',sigma=10,phase_snr=3.0,beyond_intelligence=False,DR_delta1=20, DR_delta2=20,calc_final_sigma=True,limit_dyn=1e6,final_image_poln="I",worst_DR_allowed=100,weighting='natural', ref_timeslice_done=False,ref_time='',img_number=400,safety_standard=0,want_deep=False,deep_sigma=6,final_sigma=10,need_uvsub=True,want_check_dyn=False):

#### this msname is the name of the channel ms.
#### set the imagedir and model dir to the chanel folder by default
#### before calling this function create the channel folder and go inside it.
### now yref_time will be the time for which we want to do a selfcal and then
#### get the image locations
	## Loading the centroid shift file for further update with new time slice info.
	centroid_shift_dict_file='chan_'+chan+'_centroid_shift_test.p'
	log_file=open('CASA_imaging.txt','a')
	###########################################################################
	try:
		os.mkdir(basedir+'/'+imagedir)
	except OSError:
		pass

	try:
		os.mkdir(basedir+'/'+modeldir)
	except OSError:
		pass
	
	try:
		os.mkdir(basedir+'/'+true_loc_imgs_folder)
	except OSError:
		pass	
	
	if generate_time_indices==True:
		time_beg,ref_time_index=gen_time_array(starttime,ref_timeslice_done,times,ref_time,tdt,endtime)
	else:
		time_beg=time_indices
		
	
	#print 'Imaging: ',timeranges
	msmd.open(msname)
	all_times	=	msmd.timesforscans([1])		#### assume only 1 scan. which is true in our case. 1 ms file has one scan.
	freq	 	= 	msmd.chanfreqs(0) 
	wavelength	=	299792458.0/freq 
	taper		=	str(0.95*max_baseline/wavelength[0])
	msmd.close()

	if (calc_imsize):
		cell=pixelsize(freq)
		FOV	=field_of_view(freq)
		imsize	=[num_pixels(FOV,cell)]

	if (calc_scale_flag==True):
		psf=	calc_psf(freq)
		scales=choose_scales(cell,psf)
	cellsize=str(cell)


	 
	if region_file=='':
		BOX='50,50,'+str(imsize[0]-50)+','+str(int(imsize[0]/4))
		print 'Bad region in the image for rms noise estimation: ',BOX
	else:
		BOX=''
########################################################

	bad_clean_counter=0
	print 'Bad clean error counter set: ',bad_clean_counter,'\n Image counter set to 0.'
	num_imgs=0
	for time_num,num_time_beg in enumerate(time_beg):	### n has number and i has the time
		if num_imgs==img_number:	# Limit imaging to just 100 images at  time.
			break

		part_time=all_times[num_time_beg]
		############# Checking for the fits file's existence ##################
		
		ftim=scan_timeranges(part_time).split('~')[0].split(':')
		if len(ftim[-1])==1:
			ftim[-1]='0'+ftim[-1]
			ftim[-1]=ftim[-1]+'.0'
		elif len(ftim[-1])==2 and '.' not in ftim[-1]:
			ftim[-1]+='.0'
		elif len(ftim[-1])==3 and '.' in ftim[-1]:
			ftim[-1]='0'+ftim[-1]
		if len(ftim[0])!=2:
			ftim[0]='0'+ftim[0]
		if len(ftim[1])!=2:
			ftim[1]='0'+ftim[1]
		fitstim=''.join(ftim)
		final_uvfits='Sun_Clean_Matrix_'+fitstim+'_'+chan+'_full_'+final_image_poln+'.uvfits'
		final_img_name=basedir+'/'+imagedir+'/'+'Sun_Clean_Matrix_'+fitstim+'_'+chan+'_full_'+final_image_poln+'.fits'	
		fl=glob.glob(final_img_name[:-4]+"*")
		if len(fl)==1:
			if ref_timeslice_done==False:
				ref_timeslice_done=True
			#print fl[0]+' already exists..!!'
			continue	
		badtime_fil=open('bad_times.txt','a') # Selfcal related errors. Erroneus time steps where selfcal failed with the corresponding error is saved in this file.
		timfil=open('Time_taken_till_final_imgs.txt','a')
		bad_finfil=open('Unreliable_times.txt','a') # This stores the time and DR obtained for the image at that time if the DR<worst case DR
		DRfil=open('Dynamic_Ranges.txt','a')
		bad_loc_fil=open('Bad_loc_shift_times.txt','a')
		testfil=open('Maxpix_location_info.txt','a')

		if os.path.isfile(centroid_shift_dict_file)==False:					# Create new dict file.
			print 'Initialising file, ',centroid_shift_dict_file
			location_shift=[]
			time_keys=[]
		else:										# Load earlier shift values and times
			print 'Loading existing file, ',centroid_shift_dict_file
			fulldict=pickle.load(open(centroid_shift_dict_file,'rb'))
			location_shift=fulldict.values()
			time_keys=fulldict.keys()
			del fulldict
		
		print 'Starting Imaging ..',final_img_name
		final_model_name=basedir+'/'+modeldir+'/'+'Sun_Clean_Matrix_'+fitstim+'_'+chan+'_full_'+final_image_poln+'_model.fits'
		################# Removing .ms and other files relating to this time_num #######

		msfilt='*_'+'_'.join(scan_timeranges(part_time).split('~')[0].split(':'))+'_*'
		print 'Deleting MS pertaining to ..',msfilt
		allthere=glob.glob(msfilt+'.ms')
		print ' Files before deletion: ',allthere
		os.system('rm -rf '+msfilt)
		allleft=glob.glob(msfilt+'.ms')
		print ' Files after deletion: ',allleft
		del allthere
		del allleft
		del msfilt
		del ftim
		del fl	#######################################################################################				
		begin_clock=time.time()
		

		if os.path.isfile("previous_images_records.txt")==False:
			print 'Making book keeping file..'
			bookkeeping_file=''
		else:
			print 'Reloading bookkeeping file..'
			bookkeeping_file="previous_images_records.txt"	
		beg_time=time.time()
		scratch=True
		if (ref_timeslice_done==False): # Do selfcal if this is ref time slice. ##?? Max_pix_array n dyn_range_array not passed.
			print 'Ref time slice imaging to be done.. calling do_selfcal..'
			
			msg,success,max_pix_array,dyn_range_array=do_selfcal(msname,part_time,scratch=True,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn, bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn,wavelength=wavelength)	
			pos=np.where(np.isnan(max_pix_array)==True)[0]
			max_pix_array[pos]=-1
			dyn_range_array[pos]=-1	
			np.savetxt("previous_images_records.txt",np.array([max_pix_array,dyn_range_array]).T)
			del dyn_range_array
			del max_pix_array
			
		else:	# If it is not the ref time slice do the max pixel location noting first. Then do selfcal.
			print 'Ref time slice is done so checking the max pos (Xc,Yc)_i'
			file_str="time_"+'_'.join((scan_timeranges(all_times[ref_time_index]).split('~')[0]).split(':'))+"_self*.cal"
			reliability,max_pix_coords_before_selfcal=check_source_true_loc(basedir+'/'+true_loc_imgs_folder,msname,file_str,scan_timeranges(part_time),cell,imsize,final_image_poln,fitstim,chan,taper)	
			print 'Initial max loc: (',max_pix_coords_before_selfcal[0],',',max_pix_coords_before_selfcal[1],').'
			print 'Reliability: ',reliability
			scratch= not reliability
			scratch=False
			print 'Going ahead with selfcal..'
			if reliability==False:
				bad_loc_fil.write(fitstim+'\n')
			if max_time_delta!=0:	
				print "Calling skip selfcal\n"
				msg=skip_selfcal(msname,part_time,header_dict,chan, final_image_poln, modeldir=modeldir, basedir=basedir, imagedir=imagedir, corr="XX,YY", region_file=region_file, BOX=BOX, cell=cell, imsize=imsize, scales=scales,maskfile=maskfile,phase_snr=phase_snr,refant=ref_ant,selfcal_failed=False)
			
				if msg==0:
					max_pix_coords_final=imstat(imagename=final_img_name)['maxpos'][0:2] 
					shift_factor=max_pix_coords_final-max_pix_coords_before_selfcal
					testfil.write(fitstim+'\t ('+str(max_pix_coords_before_selfcal[0])+','+str(max_pix_coords_before_selfcal[1])+').'+'\t ('+str(max_pix_coords_final[0])+','+str(max_pix_coords_final[1])+').\n')
					print 'Time: ',fitstim
					print 'Final max pixel loc found: (',max_pix_coords_final[0],',',max_pix_coords_final[1],').'
					print 'Shift : ',shift_factor,' written out to dictionary..'			
					location_shift+=[shift_factor]
					time_keys+=[fitstim]
					bad_finfil.close()
					badtime_fil.close()
					timfil.close()
					bad_loc_fil.close()
					DRfil.close()
					testfil.close()
					shift_dict=dict(zip(time_keys,location_shift))					
					pickle.dump(shift_dict,open(centroid_shift_dict_file,'wb')) # Write out new centroid shift array
					del shift_dict
					del final_uvfits
					del fitstim
					del final_img_name
					continue		
			msg,success,max_pix_array,dyn_range_array=do_selfcal(msname,part_time,scratch=scratch,corr=corr,ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,BOX=BOX, cell=cell,imsize=imsize,scales=scales,maskfile=maskfile,start_sigma=sigma,phase_snr=phase_snr,beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1,DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,intelligence_file=intelligence_file,limit_dyn=limit_dyn,bookkeeping_file=bookkeeping_file,worst_DR_allowed=worst_DR_allowed, safety_standard=safety_standard,final_sigma=final_sigma,need_uvsub=need_uvsub,want_check_dyn=want_check_dyn,wavelength=wavelength)
			pos=np.where(np.isnan(max_pix_array)==True)[0]
			max_pix_array[pos]=-1
			dyn_range_array[pos]=-1	
			np.savetxt("previous_images_records.txt",np.array([max_pix_array,dyn_range_array]).T)
			del dyn_range_array
			del max_pix_array
		
		file_str="time_"+'_'.join((scan_timeranges(part_time).split('~')[0]).split(':'))+"_self*"
		gen_str=scan_timeranges(part_time).split('~')[0]
		if msg!=2:	# Reset bad clean error counter if error message isn't 2.
			bad_clean_counter=0
		if msg==1: #split problem
			num_imgs+=1
			os.system("touch "+final_img_name[:-5]+".junk")
			log_file.write(gen_str+" "+error_msgs(msg)+"\n")
			log_file.write("exiting CASA\n")
			os.system('rm -rf '+file_str)
			shift_dict=dict(zip(time_keys,location_shift))
			pickle.dump(shift_dict,open(centroid_shift_dict_file,'wb')) # Write out the updated pickle file before ending selfcal
			bad_finfil.close()
			badtime_fil.close()
			timfil.close()
			DRfil.close()
			testfil.close()
			bad_loc_fil.close()
			return msg,num_imgs
		elif msg==2:# clean problem
			log_file.write(gen_str+" "+error_msgs(msg)+"\n")
			bad_clean_counter+=1
			os.system('rm -rf '+file_str)
			badtime_fil.write(scan_timeranges(part_time).split('~')[0]+'\t'+error_msgs(msg)+'\n')
			num_imgs+=1
			os.system("touch "+final_img_name[:-5]+".junk")
			if bad_clean_counter>10: # IF continuously bad clean error is raised for 10 times. 
				shift_dict=dict(zip(time_keys,location_shift))
				pickle.dump(shift_dict,open(centroid_shift_dict_file,'wb')) # Write out the updated pickle file before ending selfcal	
				bad_finfil.close()
				badtime_fil.close()
				timfil.close()
				DRfil.close()
				bad_loc_fil.close()
				testfil.close()
				log_file.write("exiting CASA\n")
				return msg,num_imgs	
		elif msg==3:# Many antennas flagged
			num_imgs+=1	
			log_file.write(gen_str+" "+error_msgs(msg)+"\n")
			os.system("touch "+final_img_name[:-5]+".junk")	
			os.system('rm -rf '+file_str)
			badtime_fil.write(scan_timeranges(part_time).split('~')[0]+'\t'+error_msgs(msg)+'\n')
		elif msg==4:# selfcal failed at an iteration
			log_file.write(gen_str+" "+error_msgs(msg)+"\n")
			last_msname=glob.glob(file_str+'.ms')[0]
			last_img_name=glob.glob(file_str+'.image')[-1]
			max_t=imstat(imagename=last_img_name)['max'][0]		 
			rms_t=imstat(imagename=last_img_name,box=BOX)['rms'][0]
			DR_temp=max_t/rms_t
			print 'Selfcal failed at a stage. The final image DR: ',DR_temp
			log_file.write(gen_str+' Selfcal failed at a stage. The final image DR: '+str(DR_temp)+"\n")
			badtime_fil.write(scan_timeranges(part_time).split('~')[0]+'\t'+error_msgs(msg))
			if DR_temp>worst_DR_allowed:
				num_imgs+=1
				badtime_fil.write(' Last image has reasonable DR. Hence made.'+'\n')
				log_file.write(gen_str+' Last image has reasonable DR. Hence made.'+'\n')
				DR_final=clean_as_necessary(want_deep,last_msname,final_image_poln,cell,imsize,weighting,deep_sigma,maskfile,final_img_name,final_model_name,BOX,header_dict,scales)
				exportuvfits(vis=last_msname,fitsfile=final_uvfits,datacolumn='data',overwrite=True)
				DRfil.write(fitstim+'\t'+str(DR_final)+'\n')
				os.system('rm -rf '+file_str)
			else:
				os.system("touch "+final_img_name[:-5]+".junk")
				bad_finfil.write(scan_timeranges(part_time).split('~')[0]+'\t'+error_msgs(msg)+'\n')
				badtime_fil.write('\n')
				os.system('rm -rf '+file_str)
				num_imgs+=1	
		elif msg==5:# DR decreasing
			log_file.write(gen_str+" "+error_msgs(msg)+"\n")
			last_msname=glob.glob(file_str+'.ms')[0]
			last_img_name=glob.glob(file_str+'.image')[0]
			max_t=imstat(imagename=last_img_name)['max'][0]		 
			rms_t=imstat(imagename=last_img_name,box=BOX)['rms'][0]
			DR_temp=max_t/rms_t
			print 'DR of the image was steadily decreasing.. The final DR obtained: ',DR_temp
			log_file.write(gen_str+' DR of the image was steadily decreasing.. The final DR obtained: '+str(DR_temp)+"\n")
			badtime_fil.write(scan_timeranges(part_time).split('~')[0]+'\t'+error_msgs(msg))
			if DR_temp>worst_DR_allowed:
				num_imgs+=1
				print 'Usable DR found, so saving the image..'
				log_file.write(' Last image has reasonable DR. Hence made.'+'\n')
				badtime_fil.write(' Last image has reasonable DR. Hence made.'+'\n')
				DR_final=clean_as_necessary(want_deep,last_msname,final_image_poln,cell,imsize,weighting,deep_sigma,maskfile,final_img_name,final_model_name,BOX,header_dict,scales)
				exportuvfits(vis=last_msname,fitsfile=final_uvfits,datacolumn='data',overwrite=True)
				DRfil.write(fitstim+'\t'+str(DR_final)+'\n')
				os.system('rm -rf '+file_str+".ms*")
				os.system('rm -rf '+file_str+".image")
				os.system('rm -rf '+file_str+".residual")
				os.system('rm -rf '+file_str+".mask")
				os.system('rm -rf '+file_str+".model")
				os.system('rm -rf '+file_str+".psf")
				os.system('rm -rf '+file_str+".flux")
			else:
				os.system("touch "+final_img_name[:-5]+".junk")
				bad_finfil.write(scan_timeranges(part_time).split('~')[0]+'\t'+error_msgs(msg)+'\n')
				badtime_fil.write('\n')
				os.system('rm -rf '+file_str)
				num_imgs+=1	
		elif msg==6:# Bad calibration tables.
			num_imgs+=1
			os.system("touch "+final_img_name[:-5]+".junk")	
			log_file.write(gen_str+" "+error_msgs(msg)+"\n")
			os.system('rm -rf '+file_str)
			badtime_fil.write(scan_timeranges(part_time).split('~')[0]+'\t'+error_msgs(msg)+'\n')
		elif msg==7:# maximum iteration hit
			log_file.write(gen_str+" "+error_msgs(msg)+"\n")
			last_msname=glob.glob(file_str+'.ms')[0]
			badtime_fil.write(scan_timeranges(part_time).split('~')[0]+'\t'+error_msgs(msg))
			num_imgs+=1
			if success==True:	
				log_file.write(' Last image has reasonable DR. Hence made.'+'\n')
				badtime_fil.write(' Last image has reasonable DR. Hence made.'+'\n')
				DR_final=clean_as_necessary(want_deep,last_msname,final_image_poln,cell,imsize,weighting,deep_sigma,maskfile,final_img_name,final_model_name,BOX,header_dict,scales)
				exportuvfits(vis=last_msname,fitsfile=final_uvfits,datacolumn='data',overwrite=True)
				DRfil.write(fitstim+'\t'+str(DR_final)+'\n')
				os.system('rm -rf '+file_str+".ms*")
				os.system('rm -rf '+file_str+".image")
				os.system('rm -rf '+file_str+".residual")
				os.system('rm -rf '+file_str+".mask")
				os.system('rm -rf '+file_str+".model")
				os.system('rm -rf '+file_str+".psf")
				os.system('rm -rf '+file_str+".flux")
				
			else:
				os.system("touch "+final_img_name[:-5]+".junk")
				badtime_fil.write('\n')
				os.system('rm -rf '+file_str)	
		else: # Success = True				
			last_msname=glob.glob(file_str+'.ms')[0]
			print 'Successfull imaging..Final MS: ',last_msname
			os.system("mv clean.last clean_previous.last")
			DR_final=clean_as_necessary(want_deep,last_msname,final_image_poln,cell,imsize,weighting,deep_sigma,maskfile,final_img_name,final_model_name,BOX,header_dict,scales)
			exportuvfits(vis=last_msname,fitsfile=final_uvfits,datacolumn='data',overwrite=True)
			num_imgs+=1
			DRfil.write(fitstim+'\t'+str(DR_final)+'\n')
			os.system('rm -rf '+file_str+".ms*")
			os.system('rm -rf '+file_str+".image")
			os.system('rm -rf '+file_str+".residual")
			os.system('rm -rf '+file_str+".mask")
			os.system('rm -rf '+file_str+".model")
			os.system('rm -rf '+file_str+".psf")
			os.system('rm -rf '+file_str+".flux")
		## Finding location shifts or saving the ref location if its the ref time slice.
		if ref_timeslice_done==False and msg!=0:
			 return 8,num_imgs
		elif ref_timeslice_done==True:
			if os.path.isfile(final_img_name):
				max_pix_coords_final=imstat(imagename=final_img_name)['maxpos'][0:2] 
				shift_factor=max_pix_coords_final-max_pix_coords_before_selfcal
				testfil.write(fitstim+'\t ('+str(max_pix_coords_before_selfcal[0])+','+str(max_pix_coords_before_selfcal[1])+').'+'\t ('+str(max_pix_coords_final[0])+','+str(max_pix_coords_final[1])+').\n')
				print 'Time: ',fitstim
				print 'Final max pixel loc found: (',max_pix_coords_final[0],',',max_pix_coords_final[1],').'
				print 'Shift : ',shift_factor,' written out to dictionary..'			
				location_shift+=[shift_factor]
				time_keys+=[fitstim]
			else:
				print 'Time: ',fitstim
				print 'No final image made.. '
		end_time=time.time()-beg_time
		print '&&&&&&&&&&&&&&&&&&&&'
		print 'Time taken: ',str(round(end_time/60.,2))+' min.'
		print '&&&&&&&&&&&&&&&&&&&&'
		timfil.write(fitstim+'\t'+str(round(end_time/60.,2))+' min.\n')
		if ref_timeslice_done==False:
			ref_timeslice_done=True
		del final_uvfits
		del fitstim
		del final_img_name
		del msg
		del success
		bad_finfil.close()
		badtime_fil.close()
		timfil.close()
		bad_loc_fil.close()
		DRfil.close()
		testfil.close()
		if ref_timeslice_done==True:
			shift_dict=dict(zip(time_keys,location_shift))					
			pickle.dump(shift_dict,open(centroid_shift_dict_file,'wb')) # Write out new centroid shift array
			del shift_dict
		del time_keys
		del location_shift
	pickle.dump(num_imgs,open("number_of_images_generated.p","wb"))
	return 0,num_imgs

del Out

beyond_intelligence=False        ### gives authority to code to go beyond supplied intelligence if need
				### arises. set it to True if you really know the whole code like the palm
				### of your hand. 
in_unknown=False
tdt=dt.timedelta(seconds=tdt)
home	=	os.getcwd()
var,vals=np.genfromtxt('selfcal_input.py',dtype=str,delimiter='=',autostrip=True,usecols=(0,1),unpack=1,skip_header=3,skip_footer=True)
header_dict=dict(zip(var,vals))
del var
del vals
os.chdir(basedir)

#### entend flags of calibrator

if (extend_flag):
	print 'Extending Calibrator flags....'
	extend_calibrator_flags(calibrator_cal)
log_file=open('CASA_imaging.txt','a')

chan_beg=	chan
try:
	os.mkdir(out_dir+chan_beg)
	print 'Made directory: ',out_dir+chan_beg
except OSError:
	pass

ms_out	=	out_dir+chan_beg+"/"+out_dir+str(chan_beg)+".ms"
	
if (os.path.isdir(ms_out)==False):	
	split(vis=msname,outputvis=ms_out,keepmms=True,field="",spw="0:"+chan,scan="",antenna="",correlation=corr,timerange="",intent="",array="",uvrange="",observation="",feed="",datacolumn="data",keepflags=True,width=channel_to_avg,timebin="0s",combine="")

if os.path.isfile(out_dir+chan_beg+'/.calibrator_applied.check')==False:
	print "Applying calibrator solutions"
	applycal(vis=ms_out,field="",spw="",intent="",selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",msselect="",docallib=False,callib="",gaintable=calibrator_cal,gainfield=[],interp='nearest',spwmap=[],calwt=[True],parang=False,applymode="calflag",flagbackup=False)
	os.system("touch "+out_dir+chan_beg+"/.calibrator_applied.check")
	os.chdir(out_dir+chan_beg)
	split(vis="chan_"+str(chan_beg)+".ms",outputvis="dummy.ms",keepmms=True,field="",spw="",scan="",antenna="",correlation=corr,timerange="",intent="",array="",uvrange="",observation="",feed="",datacolumn="corrected",keepflags=True,width=1,timebin="0s",combine="")
	print "Remade chan_"+str(chan_beg)+".ms with data column as corrected data after application of Calibrator solns."
	os.system("rm -rf "+"chan_"+str(chan_beg)+".ms")
	os.system("mv dummy.ms chan_"+str(chan_beg)+".ms")
else:
	os.chdir(out_dir+chan_beg)

chan_msname='chan_'+str(chan_beg)+".ms"
log_file.close()

if os.path.isfile("number_of_images_generated.p"):
	os.system("rm -rf "+"number_of_images_generated.p")

print 'Calling control selfcal'
msg,num_imgs_made=control_selfcal(basedir=basedir,header_dict=header_dict,msname=chan_msname,intelligence_file=intelligence_file,times=times,starttime=starttime,endtime=endtime,tdt=tdt,imagedir=imagedir,modeldir=modeldir,true_loc_imgs_folder=true_loc_imgs_folder, out_dir=out_dir,corr=corr,need_both_poln=need_both_poln, ref_ant=ref_ant,max_iter=max_iter,num_ant=num_ant,region_file=region_file,cell=cell, calc_imsize=calc_imsize,imsize=imsize, calc_scale_flag=calc_scale_flag, scales=scales, maskfile=maskfile,sigma=sigma,phase_snr=phase_snr, beyond_intelligence=beyond_intelligence, DR_delta1=DR_delta1, DR_delta2=DR_delta2,calc_final_sigma=calc_final_sigma,limit_dyn=limit_dyn, final_image_poln=final_image_poln, worst_DR_allowed=worst_DR_allowed, weighting=weighting,ref_timeslice_done=ref_timeslice_done, ref_time=ref_time,img_number=img_number,safety_standard=safety_standard,want_deep=want_deep,deep_sigma=deep_sigma,final_sigma=final_sigma, need_uvsub=need_uvsub,want_check_dyn=want_check_dyn)


print error_msgs(msg)

os.chdir(home)

