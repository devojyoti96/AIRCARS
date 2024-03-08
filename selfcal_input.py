import datetime as dt

basedir		=	'/media/surajit/_data2/PhD/1130644136/b187-188'
msname		=	'1130644136-b187-188.MS'
calibrator_cal=	'/media/surajit/_data2/PhD/1130633816/1130633816_b187-188.cal'  ###give the full path if not in same path as ms
#calibrator_cal=		'/Data1/surajit/PhD/decor_data/div_data/time_32943_stokesI_chan_8_00/1130633816_113-114.cal'
ncoarse		=	2  # ncoarse_chan=2 for picket fence, 12 for harmonic and 24 for continuous mode observations
freq_res	=	40 # freq resolution of the data in kHz
flag_channels	=	False
metafits	=	'/media/surajit/_data2/PhD/MWA_test_data/1062131888.metafits'
max_time_delta	=	10
imagedir	=	'All_fits_05_160kHz_uncalib' # This should be inside basedir
modeldir	=	'All_model_fits'	# This should be inside basedir
true_loc_imgs_folder	= 'All_fits_05_160kHz_true_loc'
decor_applied	=	False

safety_standard= 1   ### takes 3 values 0/1/2 0: normal steps 1: medium safety 3: highest safety

corr		=	"XX,YY"
need_both_poln	=	True   #### right now the code only takes True. It doesn't matter iof you give false
ref_ant		=	0 ### give just the antenna number. not a string
max_iter	=	500
num_ant		=	128
region_file	=	""# Leave it "" to autoset
out_dir		=	"chan_"

calc_imsize	=	True	#### if True tells code to auto calculate the imsize and cellsize
cell		=	60 #### enter in arc-secs
imsize		=	[1620]

calc_scale_flag	=	True	#### If True will calculate scales which will be used in multiscale
scales		=	[]
casapath	=	"/media/surajit/_data2/casa-release-4.6.0-el6/bin/casa"

ref_timeslice_done	=	False

weighting	=	"natural"    #### weighting to be used in final image
generate_time_indices=True
times=['01:42:00.0~01:45:58.0']
time_indices=[0]
ref_time_index=0
ref_time 	=	'01:43:00.0'
limit_dyn	=	700 #### stop selfcal if you reach this dynamic range
want_deep	=	False
final_image_poln=	"I"
maskfile	=	""
sigma		=	10	### threshold=sigma*rms
phase_snr	=	3.0	#### snr of phase solution
calc_final_sigma=	True
threshold_check	=	False	
final_sigma	=	6.0	### threshold for final selfcal. if cal_final_sigma=True, calls the automatic threshold setter. Otherwise use what the user has given.
deep_sigma	=	5
increase_threshold_no_selfcal=0.2   ### The deep_sigma will be increased by this amount for cleaning datasets where no selfcalibration has been done

beyond_intelligence=False        ### gives authority to code to go beyond supplied intelligence if need arises
DR_delta1	= 	20	## The DR improvement factor that decides if the image quality is improvement.
DR_delta2	=	20
in_unknown	=	False   ### tells the code if it is already in unknown regime. If true, it will make some informed choices.
channel_to_avg	=	4 # No. of channels to average while imaging. 4 => 4 CASA channels to avg
chan		=	'8~11'
starttime	=	dt.datetime(2017,11,27,1,42,00) # start time of observations
endtime		=	dt.datetime(2017,11,27,1,46,00)
tdt		=	0.5# Time resolution of images in seconds
selfcal_code_loc=	'/Data1/surajit/PhD/AIRCARS_stable'
intelligence_file= selfcal_code_loc+"/selfcal_intelligence.p" #Give full path
time_before_next_screen_opens=10
worst_DR_allowed=	200 # The DR limit minimum necessary for the final image to be written out
img_number	=	2 #At a time how many images you need the pipeline to make. No guarentee.  
maskfile	=	''#circle[[1238pix,1289pix],80pix]' # Mask for cleaning and model develeopment in selfcal
need_uvsub	=	False
want_check_dyn	=	False
max_baseline	=	5000  ### in metres. Only used in calculation of taper. taper=0.95*max_baseline/wavelength
extend_flag	=	False     ##### This should be the last line of this file (for developers only)
################################################################################################################################

