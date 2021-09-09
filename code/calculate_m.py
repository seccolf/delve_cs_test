from astropy.io import fits
import numpy as np 

f = fits.open('/Users/secco/Documents/projects/delve_shear/delve_cs_test/data/simulated_metacal_output.fits')

mask = (f[1].data['mcal_s2n_noshear']>10) & (f[1].data['mcal_T_ratio_noshear']>0.5)
e1 = f[1].data['mcal_g_noshear'][:,0]
e2 = f[1].data['mcal_g_noshear'][:,1]
R11 = np.mean((f[1].data['mcal_g_1p'][:,0] - f[1].data['mcal_g_1m'][:,0])/0.02)
R22 = np.mean((f[1].data['mcal_g_2p'][:,1] - f[1].data['mcal_g_2m'][:,1])/0.02)

g1 = e1/R11
g2 = e2/R22

trueg1 = 0.02
trueg2 = 0.0

print('recovered g1= %1.3f'%np.mean(g1),'+/- %1.3f'%np.std(g1))
print('trueg1=%1.3f'%trueg1)
print('recovered g2= %1.3f'%np.mean(g2),'+/- %1.3f'%np.std(g2))
print('trueg2=%1.3f'%trueg2)
