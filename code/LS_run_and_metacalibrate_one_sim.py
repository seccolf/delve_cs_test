#from astropy.io import fits
#from ngmix.medsreaders import NGMixMEDS
import numpy as np
import galsim
import fitsio
import os
import sys
from _step import _run_metacal as run_metacal
sys.path.insert(0,'/home/secco/delve_cs_test/code/newish_metacal/tests')
#sys.path.insert(0,'/Users/secco/Documents/projects/delve_shear/delve_cs_test/code/newish_metacal')
#sys.path.insert(0,'/Users/secco/Documents/projects/delve_shear/delve_cs_test/code/newish_metacal/metacal')
#import test_metacal
from Metacal_test_V1_with_dither import replace_all_cutouts_by_sims 
import time
################################################
################################################

# SETUP ACTUAL RUN!
NTASKS = int(os.environ['SLURM_NTASKS'])
#start looping over possible seed numbers based on how many tasks are running:

for loop_counter in np.arange(100): #each process will likely not have time to do all of these, but it's ok
    SEED = int(os.environ['SLURM_PROCID']) + 50 + loop_counter*NTASKS
    PROCID = str(int(os.environ['SLURM_PROCID'])+1)
    print('PROCESS %s is starting calculation on SEED %d'%(PROCID,SEED))
    meds_dir = '/home/secco/scratch-midway2/delve/simulated_meds/'
    input_medsname = 'TESTBED_DES0954+0252_r4918p01_i_meds-DECADE_copy'+PROCID+'.fits'

    #now replace all of the objects in this by simulated objects
    time1 = time.time()
    replace_all_cutouts_by_sims(meds_dir+input_medsname, SEED, +0.02, 0.00, verbose=True)
    time2 = time.time()
    print('replace_cutout_by_sims took %1.2f minutes'%((time2-time1)/60.0))

    #now metacalibrate the positive g1 output:
    output_dir = '/home/secco/scratch-midway2/delve/simulated_mcal/'
    time3 = time.time()
    output = run_metacal([meds_dir+input_medsname],SEED)
    time4 = time.time()
    fitsio.write(output_dir+'metacal_output_seed'+str(SEED)+'_positive.fits', output, clobber=True)
    print('positive metacal run took %1.2f minutes'%((time4-time3)/60.0))


    #now replace all of the objects in this by simulated objects
    replace_all_cutouts_by_sims(meds_dir+input_medsname,SEED, -0.02, 0.00, verbose=True)
    #now metacalibrate the negative g1 output:
    output_dir = '/home/secco/scratch-midway2/delve/simulated_mcal/'
    output = run_metacal([meds_dir+input_medsname],SEED)
    fitsio.write(output_dir+'metacal_output_seed'+str(SEED)+'_negative.fits', output, clobber=True)




