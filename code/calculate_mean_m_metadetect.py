
from astropy.io import fits
import numpy as np

import sys
from tqdm import tqdm

start, stop, V = int(sys.argv[1]), int(sys.argv[2]), str(sys.argv[3])
N = stop - start

gplus, gminus = [0]*N, [0]*N

for i in tqdm(range(start, stop)):

	if (i < 50) | (i > 100):
		gplus[i - start]  = fits.open('/scratch/midway2/dhayaa/Metacal_test_V%s/gplus_V%d_metacal_output.fits'%(V, i))
		gminus[i - start] = fits.open('/scratch/midway2/dhayaa/Metacal_test_V%s/gminus_V%d_metacal_output.fits'%(V, i))

	elif (i >= 50) & (i <= 97):
		gplus[i - start]  = fits.open('/project2/chihway/secco/mcal/metacal_output_seed%d_positive.fits'%i)
		gminus[i - start] = fits.open('/project2/chihway/secco/mcal/metacal_output_seed%d_negative.fits'%i)

	#EDGE CASE: TEMPORARY
	elif (i > 97) & (i <= 100):
		gplus[i - start]  = fits.open('/project2/chihway/secco/mcal/metacal_output_seed96_positive.fits')
		gminus[i - start] = fits.open('/project2/chihway/secco/mcal/metacal_output_seed96_negative.fits')

	# print("Done with", i)

#mask = (gplus[1].data['mcal_s2n_noshear']>10) & (gplus[1].data['mcal_T_ratio_noshear']>0.5)

#print(mask.sum(), mask.size)

#mask = (gminus[1].data['mcal_s2n_noshear']>10) & (gminus[1].data['mcal_T_ratio_noshear']>0.5)

#print(mask.sum(), mask.size)


e1_plus_array = np.concatenate([gplus[i][1].data['mcal_g_noshear'][:,0] for i in range(N)])
e2_plus_array = np.concatenate([gplus[i][1].data['mcal_g_noshear'][:,1] for i in range(N)])

e1_minus_array = np.concatenate([gminus[i][1].data['mcal_g_noshear'][:,0] for i in range(N)])
e2_minus_array = np.concatenate([gminus[i][1].data['mcal_g_noshear'][:,1] for i in range(N)])


R11_plus_array = np.concatenate([(gplus[i][1].data['mcal_g_1p'][:,0] - gplus[i][1].data['mcal_g_1m'][:,0])/0.02 for i in range(N)])
R22_plus_array = np.concatenate([(gplus[i][1].data['mcal_g_2p'][:,1] - gplus[i][1].data['mcal_g_2m'][:,1])/0.02 for i in range(N)])


R11_minus_array = np.concatenate([(gminus[i][1].data['mcal_g_1p'][:,0] - gminus[i][1].data['mcal_g_1m'][:,0])/0.02 for i in range(N)])
R22_minus_array = np.concatenate([(gminus[i][1].data['mcal_g_2p'][:,1] - gminus[i][1].data['mcal_g_2m'][:,1])/0.02 for i in range(N)])

del gplus, gminus

nBoot = 1000

e1_plus  = np.zeros(nBoot)
e1_minus = np.zeros(nBoot)
e2_plus  = np.zeros(nBoot)
e2_minus = np.zeros(nBoot)

R11_plus  = np.zeros(nBoot)
R11_minus = np.zeros(nBoot)
R22_plus  = np.zeros(nBoot)
R22_minus = np.zeros(nBoot)

print(e1_plus_array.shape)
print(R11_plus_array.shape)
print(e1_minus_array.shape)
print(R11_minus_array.shape)


for i in tqdm(range(nBoot)):

	index_plus  = np.random.randint(0, len(e1_plus_array),  len(e1_plus_array))
	index_minus = np.random.randint(0, len(e1_minus_array), len(e1_minus_array))

	e1_plus[i]  = np.mean(e1_plus_array[index_plus])
	e1_minus[i] = np.mean(e1_minus_array[index_minus])
	e2_plus[i]  = np.mean(e2_plus_array[index_plus])
	e2_minus[i] = np.mean(e2_minus_array[index_minus])

	R11_plus[i]   = np.mean(R11_plus_array[index_plus])
	R11_minus[i]  = np.mean(R11_minus_array[index_minus])
	R22_plus[i]   = np.mean(R22_plus_array[index_plus])
	R22_minus[i]  = np.mean(R22_minus_array[index_minus])


#print(len(e1_plus))
R11 = (R11_plus + R11_minus)/2
R22 = (R22_plus + R22_minus)/2

g1 = (e1_plus - e1_minus)/(2*R11)
g2 = (e2_plus - e2_minus)/(2*R22)

trueg1 = 0.02
trueg2 = 0.0

print('recovered g1= %1.6f'%np.mean(g1),'+/- %1.6f'%np.std(g1))
print('trueg1=%1.3f'%trueg1)
print('recovered g2= %1.6f'%np.mean(g2),'+/- %1.6f'%np.std(g2))
print('trueg2=%1.3f'%trueg2)
print('---------------------------------')
print('recovered m= %1.6f'%np.mean(g1/trueg1 - 1),'+/- %1.6f'%np.std(g1/trueg1 - 1))
