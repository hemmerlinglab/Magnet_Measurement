import numpy as np
import matplotlib.pyplot as plt
import magpylib as magpy
import scipy.optimize
import pickle

import lmfit
from lmfit import Minimizer, Parameters, report_fit

from functools import partial

from helper_functions import *



# ideal field
B0 = 537.0 # Gs
L0 = 26.6 # cm
Bg = B0/2 # Gs

zoffset = (2.5/2+.375)*25.4 # mm, radius of nipple plus half the magnet
xoffset = .25*25.4 # mm, magnet radius for double stacked

oris = [1,1,1,1,1,1,1,1,-1,-1,-1, -1] # orientation of magnets



zlocs = np.zeros(len(oris)) 
ylocs = np.arange(0,len(zlocs))*25.4 # mm, beam line locations

shift_z = 10.0 * np.ones(len(zlocs))
no_of_iterations = 30

Bideal = get_ideal_field(ylocs, B0, Bg, L0)

shift_z_arr = np.zeros([len(zlocs), no_of_iterations])
B_arr = np.zeros([len(zlocs), no_of_iterations])
dBdz_arr = np.zeros([len(zlocs), no_of_iterations])
zlocs_arr = np.zeros([len(zlocs), no_of_iterations])

for k in range(no_of_iterations):
    print(k)

    (delta, shift_z, B, dBdz) = do_iter(Bideal, ylocs, zlocs, xoffset, zoffset, oris, shift_z)
    zlocs += delta

    ind = np.where(zlocs < 0.0)
    if len(ind[0])>0:
        zlocs[ind[0][0]] = 0.0

    shift_z_arr[:, k] = shift_z
    B_arr[:, k] = B
    dBdz_arr[:, k] = dBdz
    zlocs_arr[:, k] = zlocs



f = open('store.pckl', 'wb')
pickle.dump([ylocs, zlocs, B0, Bg, L0, shift_z_arr, B_arr, dBdz_arr, zlocs_arr, xoffset, zoffset, oris], f)
f.close()


