from astropy.io import fits
import numpy as np

gplus  = fits.open('gplus_metacal_output.fits')
gminus = fits.open('gminus_metacal_output.fits')

#mask = (gplus[1].data['mcal_s2n_noshear']>10) & (gplus[1].data['mcal_T_ratio_noshear']>0.5)

#print(mask.sum(), mask.size)

#mask = (gminus[1].data['mcal_s2n_noshear']>10) & (gminus[1].data['mcal_T_ratio_noshear']>0.5)

#print(mask.sum(), mask.size)


e1_plus_array = gplus[1].data['mcal_g_noshear'][:,0]
e2_plus_array = gplus[1].data['mcal_g_noshear'][:,1]

e1_minus_array = gminus[1].data['mcal_g_noshear'][:,0]
e2_minus_array = gminus[1].data['mcal_g_noshear'][:,1]


R11_plus_array = (gplus[1].data['mcal_g_1p'][:,0] - gplus[1].data['mcal_g_1m'][:,0])/0.02
R22_plus_array = (gplus[1].data['mcal_g_2p'][:,1] - gplus[1].data['mcal_g_2m'][:,1])/0.02


R11_minus_array = (gminus[1].data['mcal_g_1p'][:,0] - gminus[1].data['mcal_g_1m'][:,0])/0.02
R22_minus_array = (gminus[1].data['mcal_g_2p'][:,1] - gminus[1].data['mcal_g_2m'][:,1])/0.02


nBoot = 1000

e1_plus  = np.zeros(nBoot)
e1_minus = np.zeros(nBoot)
e2_plus  = np.zeros(nBoot)
e2_minus = np.zeros(nBoot)
for i in range(nBoot):

        #Need two different indexes because len(e1_plus_array) != len(e2_plus_array)
        #We miss ~8 galaxies (0.02%) of the total sample of ~29,000 galaxies
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

R11_plus  = np.zeros(nBoot)
R11_minus = np.zeros(nBoot)
R22_plus  = np.zeros(nBoot)
R22_minus = np.zeros(nBoot)

print(e1_plus_array.shape)
#print(R11_plus_array.shape)
print(e1_minus_array.shape)
#print(R22_plus_array.shape)


