from astropy.io import fits
from ngmix.medsreaders import NGMixMEDS
import numpy as np
import galsim
import fitsio

#Functions
def flux_to_mag(F):
    return -2.5*np.log10(F)+30.0

def mag_to_flux(m):
    return 10**(-0.4*(m-30))

def make_sim_one_cutout(*, seed, g1, g2, open_meds, obj_ind, cutout):
    m = open_meds
    noise_rng  = np.random.RandomState(seed=seed)
    dither_rng = np.random.RandomState(seed=seed)
    gal = galsim.Exponential(half_light_radius=0.5).shear(g1=g1, g2=g2) #let's do sheared exponentials again

    #fwhm=rng.uniform(low=0.8,high=1.5) #will use a random FWHM for each sim just because
    psf = galsim.Gaussian(fwhm=0.9)
    obj = galsim.Convolve([gal, psf])

    dim = m[1]['box_size'][obj_ind] #size of the image box

    cutout_row = m[1]['cutout_row'][obj_ind][cutout] + 1
    cutout_col = m[1]['cutout_col'][obj_ind][cutout] + 1
    dither = dither_rng.uniform(size=2, low=-1, high=1)
    #scale = 0.263 #the size of a DECam pixel in arcsec on the sky

    #get WCS:
    #jacobian: column in MEDS is x, row is y
    dudx = m[1]['dudcol'][obj_ind][cutout]
    dudy = m[1]['dudrow'][obj_ind][cutout]
    dvdx = m[1]['dvdcol'][obj_ind][cutout]
    dvdy = m[1]['dvdrow'][obj_ind][cutout]
    wcs = galsim.JacobianWCS(dudx, dudy, dvdx, dvdy)

    #from matt: drawImage(nx=stampsize, ny=stampsize, wcs=jac_wcs, center=galsim.PositionD(x=cutout_col, y=cutout_row))
    #flux = m[1].data['flux'][obj_ind]
    im = obj.drawImage(nx=dim, ny=dim,  wcs=wcs,
                       center = galsim.PositionD(x=cutout_col, y=cutout_row),
                       offset = galsim.PositionD(x=dither[0],  y=dither[1])).array * mag_to_flux(12.0) #multiplying by flux here, is that what I really want?

    im_start_row = m[1]['start_row'][obj_ind][cutout]
    im_end_row = im_start_row + dim**2
    weight=m[5][im_start_row:im_end_row].reshape(dim,dim) #get the weight of each cutout, will use this to add noise to the simulated image

    with np.errstate(divide='ignore', invalid = 'ignore'):
        nse=np.sqrt(1.0/weight)

    im += noise_rng.normal(size=im.shape, scale=nse)
    im[weight==0.0]=0.0 #setting the flux to zero where the weight is zero

    bmask_im = m['bmask_cutouts'][im_start_row:im_end_row]
    bmask_im[bmask_im!=2**30] = 0 #setting bmask to zero everywhere it is not 2^30


    psf_row = m[1]['psf_row_size'][obj_ind][cutout] #size of the PSF (rows) for each cutout
    psf_col = m[1]['psf_col_size'][obj_ind][cutout] #size of the PSF (cols) for each cutout
    assert psf_row==psf_col, "PSF is not a square, which they say is not cool"

    psf_start_row = m[1]['psf_start_row'][obj_ind][cutout]
    psf_end_row = psf_start_row + psf_row*psf_col
    psf_center_y = m[1]['psf_cutout_row'][obj_ind][cutout]+1
    psf_center_x = m[1]['psf_cutout_col'][obj_ind][cutout]+1
    psf_im = psf.drawImage(nx=psf_row, ny=psf_col, wcs=wcs,
                           center=galsim.PositionD(x=psf_center_x, y=psf_center_y)).array

    return im.ravel(), psf_im.ravel(), bmask_im.ravel(), im_start_row, im_end_row, psf_start_row, psf_end_row


def replace_all_cutouts_by_sims(input_medsname, seed, g1, g2, verbose=False):
    m = fitsio.FITS(input_medsname,mode='rw') #opens the input file
    Nobjects = m[1].get_nrows() #gets number of objects (rows) in the file
    for i_obj in range(Nobjects): #loops over objects in the file
        Ncutout = m[1]['ncutout'][i_obj]
        if verbose and i_obj%1000==0: print("Object %d has %d cutouts to be replaced"%(i_obj,Ncutout))
        for j_cut in range(Ncutout): #loops over cutout in each object in the file
            im, psf, bmask, im_start, im_end, psf_start, psf_end = make_sim_one_cutout(seed=seed, g1= g1, g2=g2,
                                                                                       open_meds=m, obj_ind=i_obj, cutout=j_cut)

            m['image_cutouts'].write(im, im_start)
            m['bmask_cutouts'].write(bmask, im_start)
            m['psf'].write(psf, psf_start)
    m.close()

if __name__ == "__main__":
    
  import sys
  sys.path.insert(0,'/home/dhayaa/Desktop/DECADE/delve_cs_test/code/newish_metacal/tests')
  sys.path.insert(0,'/home/dhayaa/Desktop/DECADE/delve_cs_test/code/newish_metacal')
  sys.path.insert(0,'/home/dhayaa/Desktop/DECADE/delve_cs_test/code/newish_metacal/metacal')

  import test_metacal
    
  seed = int(sys.argv[1])
  
  version = '1.1'

  input_medsname = r'/scratch/midway2/dhayaa/Metacal_test_V%s/gplus_V%d.fits'%(version, seed)
  replace_all_cutouts_by_sims(input_medsname, seed,  0.02, 0.00, verbose=True)

  input_medsname = r'/scratch/midway2/dhayaa/Metacal_test_V%s/gminus_V%d.fits'%(version, seed)
  replace_all_cutouts_by_sims(input_medsname, seed, -0.02, 0.00, verbose=True)

  #Generate output files

  from _step import _run_metacal as run_metacal

  filename = ['/scratch/midway2/dhayaa/Metacal_test_V%s/gplus_V%d.fits'%(version, seed)]
  output = run_metacal(filename, 8)
  fitsio.write('/scratch/midway2/dhayaa/Metacal_test_V%s/gplus_V%d_metacal_output.fits'%(version, seed), output, clobber=True)

  filename = ['/scratch/midway2/dhayaa/Metacal_test_V%s/gminus_V%d.fits'%(version, seed)]
  output = run_metacal(filename, 8)
  fitsio.write('/scratch/midway2/dhayaa/Metacal_test_V%s/gminus_V%d_metacal_output.fits'%(version, seed), output, clobber=True)
