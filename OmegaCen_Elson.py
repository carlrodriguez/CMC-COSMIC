import numpy as np
import os
import sys
from cosmic.sample import InitialCMCTable
from cosmic.sample.sampler import cmc
from astropy import units
from astropy import constants as const
import cosmic
print(cosmic.__version__)


RUNDIR = '/proj/rodriguezlab/users/elena/OmegaCen/CMC-COSMIC/OmegaCen/'
#tidal radius info, ~208.60 pc  : https://people.smp.uq.edu.au/HolgerBaumgardt/globular/parameter.html

#parameters commonly varied 
r_tidal = 208.60	# tidal radius in pc
r_v = 4 # virial radius in pc
mbh = 10000 # Mass of central BH in Msun
mbh_str = '1e4' # string for directory name
seed = 1234 
BHLC_flag = 1 # 0 (off) or 1 (on) 
alpha3 = 2.3
gamma = [3.0] #parameter for Elson Profile

if BHLC_flag == 0:
    Params_file = 'ElsonProfile_OmegaCen.ini'
    bhlc_flag = ''
else:
    Params_file = 'ElsonProfile_OmegaCen_bhlc.ini'
    bhlc_flag = '_bhlc'

for g in gamma:
	Singles, Binaries = InitialCMCTable.sampler('cmc', binfrac_model=.1, primary_model='custom',cluster_profile='elson', gamma=g,
											mcuts=[0.08,0.5,1.0,150.], alphas=[-1.3,-2.3,-alpha3],
											ecc_model='thermal', porb_model='log_uniform', qmin=-1.0,
												met=0.00017, size=12000000,
											params='./' + Params_file, virial_radius=r_v,
											tidal_radius=r_tidal/r_v,
											seed=seed, sample_porb_first=True, set_radii_with_BSE=True,
											central_bh=mbh)
	Mcl = Singles.mass_of_cluster                                                                             
	print('Mcl = %s Msun'%Mcl)                         

	#This normalizes to N-body units
	InitialCMCTable.ScaleToNBodyUnits(Singles,Binaries)
	InitialCMCTable.write(Singles, Binaries, filename="elson.hdf5")
				
	DIR=RUNDIR+ 'N12e6_MBH' + mbh_str + '_RV' + str(r_v) + '_gamma_' + str(g) + bhlc_flag
	os.system('mkdir -p %s'%DIR)
	print(DIR)
	# # COPY SUBMIT FILE AND CMC FILE
	os.system('scp submit_job_elson.sh %s/.'%DIR)
	os.system('scp elson.hdf5 %s/.'%DIR)
	os.system('scp OmegaCen_Elson.py %s/.'%DIR)
	os.system('scp %s %s/.'%(Params_file, DIR))
	os.system('rm elson.hdf5')
	os.system('scp /proj/rodriguezlab/users/elena/OmegaCen/CMC-COSMIC/CMC/bin/cmc %s/.'%DIR)

