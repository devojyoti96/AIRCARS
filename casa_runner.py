import pexpect
import pickle
import sys
import os
import selfcal_input
reload(selfcal_input)
del selfcal_input
from selfcal_input import *

filename='selfcal_LTS.py'
fold=os.getcwd()
chan=fold.split('/')[-1].split('_')[1]
max_img_one_run=2
num_img=max_img_one_run
max_iteration=500/max_img_one_run
iteration=0
while num_img>=max_img_one_run and iteration<max_iteration:
	flag=1
	if iteration!=0:
		if os.path.isfile("number_of_images_generated.p")==False:	
			flag=0
			os.system("echo 'Error in spw:  "+chan+"..'| mail -s 'Job done' atul@ncra.tifr.res.in")
			break
	prompt = [r'CASA <[0-9]+>:','libgomp',pexpect.EOF]
	cmd=['--nologger','--nogui','--colors=NoColor','--nologfile']
	child=pexpect.spawn(casapath,cmd)
	child.logfile =sys.stdout
	child.expect(prompt,timeout=30)
	child.sendline("execfile(\'"+filename+"\')")
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
	num_img=pickle.load(open("number_of_images_generated.p","rb"))
	os.system("sed -i 's/extend_flag\t=\tTrue/extend_flag\t=\tFalse/g' selfcal_input.py")
	iteration+=1

chan=os.getcwd().split('/')[-1].split('_')[-1]
if flag==1:
	os.system("echo 'Run imaging for images in "+chan+"..\nTotal images: \nFinal message: .'| mail -s 'Job done' surajit@ncra.tifr.res.in")
