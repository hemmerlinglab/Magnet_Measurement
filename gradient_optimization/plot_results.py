import numpy as np
import pickle
import matplotlib.pyplot as plt

from helper_functions import *

f = open('store.pckl', 'rb')
(ylocs, zlocs, B0, Bg, L0, shift_z_arr, B_arr, dBdz_arr, zlocs_arr, xoffset, oris) = pickle.load(f)

ys = np.linspace(0, 300.0, 100.0) # in mm

Bfield = get_ideal_field(ys, B0, Bg, L0)

(M, D, L) = get_bfield(ys, zlocs, ylocs, xoffset, zoffset, oris, return_magnets = True)

plt.figure()
plt.subplot(2,2,1)
plt.plot(ylocs, B_arr)
plt.plot(ys, Bfield, 'r--')

plt.subplot(2,2,2)
plt.plot(ys, Bfield, 'r--')
plt.plot(ys, get_bfield(ys, zlocs, ylocs, xoffset, zoffset, oris), 'b')

plt.subplot(2,2,3)
plt.bar(np.arange(0,len(zlocs)), zlocs_arr[:, -1])

for k in range(len(zlocs)):
    plt.text(k - 0.5, 5.0 + zlocs_arr[k, -1], "{0:3.2f}".format(zlocs_arr[k, -1]))

plt.ylim(0, 60)

plt.ylabel('Magnet distances')


plt.subplot(2,2,4)

#pos = np.arange(0, len(zlocs))

pos = ylocs/10.0
for k in range(len(zlocs)):
    if k in [6,7,8]:
        mycolor = 'r'        
    else:        
        mycolor = 'b'
    if oris[k] > 0:
        leg = 'N'
    else:
        leg = 'S'

    
    plt.bar(pos[k], L[k], bottom = 0.0, width = D[k]/10.0, color = mycolor)
    
    my_offset = 60.0

    plt.bar(pos[k], L[k], bottom = my_offset, width = D[k]/10.0, color = mycolor)
    
    plt.text(pos[k] - 0.125, L[k] + 1.0, leg)

plt.show()


