from _step import _run_metacal as run_metacal
import fitsio

filename = ['/scratch/midway2/secco/delve/meds/DES0954+0252_r4918p01_i_meds-DECADE.fits.fz']
output = run_metacal(filename,8)
fitsio.write('./test_metacal_output.fits', output, clobber=True)
