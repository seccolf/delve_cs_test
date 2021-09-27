from astropy.io import fits
from ngmix.medsreaders import NGMixMEDS
import numpy as np
import galsim
import fitsio
import os
import sys
from _step import _run_metacal as run_metacal
sys.path.insert(0,'/Users/secco/Documents/projects/delve_shear/delve_cs_test/code/newish_metacal/tests')
sys.path.insert(0,'/Users/secco/Documents/projects/delve_shear/delve_cs_test/code/newish_metacal')
sys.path.insert(0,'/Users/secco/Documents/projects/delve_shear/delve_cs_test/code/newish_metacal/metacal')
import test_metacal
import Metacal_test_V1 

################################################
################################################

# SETUP ACTUAL RUN!
SEED = int(os.environ['SLURM_PROCID']) +50
meds_dir = '/home/secco/scratch-midway2/delve/simulated_meds/'
input_medsname = 'TESTBED_DES0954+0252_r4918p01_i_meds-DECADE.fits' #missing a '.fits' for now

#now replace all of the objects in this by simulated objects
replace_all_cutouts_by_sims(input_medsname=meds_dir+input_medsname,SEED, +0.02, 0.00, verbose=True)
#now metacalibrate the positive g1 output:
output_dir = '/home/secco/scratch-midway2/delve/simulated_mcal/'
output = run_metacal(meds_dir+input_medsname,8)
fitsio.write(output_dir+'metacal_output_seed'+str(SEED)+'_positive.fits', output, clobber=True)

#now replace all of the objects in this by simulated objects
replace_all_cutouts_by_sims(input_medsname=meds_dir+input_medsname,SEED, -0.02, 0.00, verbose=True)
#now metacalibrate the positive g1 output:
output_dir = '/home/secco/scratch-midway2/delve/simulated_mcal/'
output = run_metacal(meds_dir+input_medsname,8)
fitsio.write(output_dir+'metacal_output_seed'+str(SEED)+'_negative.fits', output, clobber=True)

#will take the SLURM process ID as the random seed - be careful here as followup runs WILL have repeated SEEDS


