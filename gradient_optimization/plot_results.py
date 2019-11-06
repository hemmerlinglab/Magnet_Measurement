import numpy as np
import pickle
import matplotlib.pyplot as plt

from helper_functions import *

f = open('store.pckl', 'rb')
(ylocs, zlocs, B0, Bg, L0, shift_z_arr, B_arr, dBdz_arr, zlocs_arr, xoffset, zoffset, oris) = pickle.load(f)

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

pos = np.arange(0, len(zlocs))
plt.bar(pos, D)

plt.show()


