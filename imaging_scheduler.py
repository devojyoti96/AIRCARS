import os,pexpect,sys,glob,numpy as np,time
import selfcal_input

reload(selfcal_input)

from selfcal_input import *

def run_casa_cmd(casatask):	
	prompt = [r'CASA <[0-9]+>:','libgomp',pexpect.EOF]
	cmd=['--nologger','--nogui','--colors=NoColor','--nologfile']
	child=pexpect.spawn(casapath,cmd)
	child.logfile =sys.stdout	
	child.expect(prompt,timeout=30)
	child.sendline(casatask)
	while(True):
		out=child.expect(prompt,timeout=None)
		if out==2:
			isalive=child.isalive()
			if (isalive):
				continue
			else:
				raise pexpect.EOF ("child died an unnatural death")
		elif out==1:
			isalive=child.isalive()
			if (isalive):
				child.sendline("exit")
			else:
				raise RuntimeError ("libgomp error. Probable resource crunch in threads")
		else:
			child.sendline("exit")
			break
	return

def flag_MWA(ncoarse_chan,minchan):
	'''
		A function to generate the list of coarse-channels edges + the 
		central channel in each coarse channel to be flagged.
		Assumes a 40 kHz channel width and hence 32 spectral channels per 
		coarse channel. If needed, can be made more flexible by reading 
		this information from MS itself. All the data we have so far, is 
		at 40 kHz spectral resolution.

		Divya 07Apr2016
	'''

	M = (minchan/4)-1 # No. of channels to be flagged at the start of the coarse channel
	N = (minchan/4)-1 # No. of channels to be flagged at the tail of the coarse channel
	a = 16	# The central channel which occassionally shows the DC spike
	b = minchan-N
	c = minchan+(M-1)
	CHAN_FLAG_STR='0:'
	i = 0
	ch0 = 0
	ch1 = minchan
	while i < ncoarse_chan:
		# The 0th coarse channel requires special treatment (one less ';')
		if ch0 == 0:
			CHAN_FLAG_STR=CHAN_FLAG_STR+str(ch0)+'~'+str(ch0+M)+';'+str(ch1-(N+1))+'~'+str(ch1-1)
	#			CHAN_FLAG_STR=CHAN_FLAG_STR+str(ch0)+'~'+str(ch0+M)+';'+str(a)+';'+str(ch1-(N+1))+'~'+str(ch1)			
		else:
			CHAN_FLAG_STR=CHAN_FLAG_STR+';'+str(ch0)+'~'+str(ch0+M)+';'+str(ch1-(N+1))+'~'+str(ch1-1)		
	#			CHAN_FLAG_STR=CHAN_FLAG_STR+';'+str(ch0)+'~'+str(ch0+M)+';'+str(a)+';'+str(ch1-(N+1))+'~'+str(ch1)
	#		a = a + 32
		ch0 = ch0 + minchan
		ch1 = ch1 + minchan
		i = i + 1

	return CHAN_FLAG_STR

avg_chn=channel_to_avg
reload_selfcal_LTS=True # Do you want to reload the selfcal_LTS.py in every chan_XX/ folder XX=8~11 , 12~15 etc..

##########################################################################################

home=os.getcwd()
spws=[]

minchan=int(1.28*10**3/freq_res) # 1.28MHz is a course channel width of MWA
M = (minchan/4)-1    # This is the last flagged channel at the start.Also the half-1 of the net allowed channel width.
a = (M+1)*2       # This is the middle of a coarse channel which is bad. 
flagwidth=a    # Width of the channels that are flagged in total at the edge of 2 adjoint channels. This is = minchan
DEF_STR='0:'
strt_chan=M+1        # 1st useful/unflagged channel
#print minchan,a,M,flagwidth
Total_channels=1.28*10**3*ncoarse/freq_res



while strt_chan<Total_channels:
	spw=DEF_STR+str(strt_chan)+'~'+str(a-1)+';'
	spw=spw+str(a+1)+'~'+str(strt_chan+flagwidth-1)+';'
	strt_chan=strt_chan+minchan  # Moving the start of good channels to the next block of coarse channel
	a+=minchan		# moving a to the middle of the next coarse channel
	DEF_STR=spw
spw=spw[:-1]
fqrs=spw.split(':')[1].split(';')


#Creating spw windows for channel averaging
indx=[]
count=0
for i in np.arange(len(fqrs)): 
	tmp=fqrs[i]
	ok=1
	ptr=int(tmp.split('~')[0])
	endptr=int(tmp.split('~')[1])
	while ok==1: 
		stf=ptr
		nxt=min(ptr+avg_chn-1,endptr)
		spws+=[str(stf)+'~'+str(nxt)]
		indx+=[count]
		count+=1
		ptr+=avg_chn
		if ptr>=endptr:
			ok=0
os.chdir(basedir)


#spws	=	spws[0:1]
#spws	=	['8~8']
spws	=	['8~11']

if flag_channels==True:
	minchan=int(1.28*10**3/freq_res)
	CHAN_FLAG_STR = flag_MWA(ncoarse,minchan)	
	print 'Flagging Data...',CHAN_FLAG_STR
	
	run_casa_cmd("flagdata(vis=\'"+msname+"\',spw=\'"+CHAN_FLAG_STR+"\',flagbackup=False)")

#run_casa_cmd("fixvis(vis=\'"+msname+"\',outputvis=\'"+msname+"\',phasecenter=\'J2000 16h45m00.119s -23d32m05.687s\')")	
os.system("cp "+home+"/decor.py "+basedir)
os.system("cp "+home+"/selfcal_input.py "+basedir)
run_casa_cmd("execfile(\'decor.py\')")

while (os.path.isfile('.decor_applied.check')==False):
	time.sleep(60)


if reload_selfcal_LTS ==True:
	for spw in spws:
		fold='chan_'+spw
		try:
			os.mkdir(fold)
		except OSError:
			pass
		os.system('cp '+selfcal_code_loc+'/selfcal_*.p* '+selfcal_code_loc+'/casa_runner.py  '+fold)
		ms_out	=	out_dir+spw+"/"+out_dir+str(spw)+".ms"
		if (os.path.isdir(ms_out)==False):  
			print  "split(vis=\'"+msname+"\',outputvis=\'"+ms_out+"\',datacolumn=\'data\',width="+str(channel_to_avg)+",spw=\'0:"+spw+"\')"
			run_casa_cmd("split(vis=\'"+msname+"\',outputvis=\'"+ms_out+"\',datacolumn=\'data\',width="+str(channel_to_avg)+",spw=\'0:"+spw+"\')")
		print 'Made ',ms_out,' from ',msname,'\n'
		os.chdir(fold)	
		fil=open('selfcal_input.py','rw+')
		lines=fil.readlines()
		for line_num,line in enumerate(lines):
			if 'chan\t\t=' in line:
				lines[line_num]="chan\t\t=\t'"+spw+"'\n"
			if 'extend_flag' in line:				
				if (spw!=spws[0]):
					lines[line_num]="extend_flag\t=\tFalse"
				else:
					if "True" in line:
						lines[line_num]="extend_flag\t=\tTrue"
					else:
						lines[line_num]="extend_flag\t=\tFalse"
				
		fil.seek(0)
		fil.writelines(lines)
		fil.close()
		os.chdir(basedir)


#spws=np.loadtxt('/Data1/atul/20141103/Solar_Event/1099030336/spws.dat',dtype='str')
#spws=spws[31:]
folds=['chan_'+i for i in spws]

os.chdir(basedir)
for i in range(len(folds)):
	os.chdir(folds[i])
	chunks=basedir.split('/')
	os.system('screen -S '+chunks[-2]+"_"+chunks[-1]+"_"+spws[i]+' -X quit')	
	os.system('screen -mdS '+chunks[-2]+"_"+chunks[-1]+"_"+spws[i])
	print 'Made Screen :',chunks[-2]+"_"+chunks[-1]+"_"+spws[i]
	os.system('screen -S '+chunks[-2]+"_"+chunks[-1]+"_"+spws[i]+' -X stuff "python casa_runner.py\n"')
	time.sleep(time_before_next_screen_opens)
	os.chdir(basedir)

print 'All screens are spawned..'
