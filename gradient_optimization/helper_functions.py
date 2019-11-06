import numpy as np
import matplotlib.pyplot as plt
import magpylib as magpy
import scipy.optimize

import lmfit
from lmfit import Minimizer, Parameters, report_fit

from functools import partial


def Cylmag(loc,ori,d,l,m):
    return magpy.source.magnet.Cylinder(mag=[0,0,ori*m],dim=[d,l],pos=loc)


def get_bfield(ys, zlocs, ylocs, xoffset, oris, return_magnets = False):
    sources = {}

    M = 6335.0/(4.0*1.26) * np.ones(len(zlocs))
    D = 0.50 * 25.4        * np.ones(len(zlocs))
    L = 0.75 * 25.4        * np.ones(len(zlocs))

    D[[0,1,11]] = 0.8 * 25.4
    L[[0,1,11]] = 0.8 * 25.4

    # needs fixing
    zoffset = 2.5/2 * 25.4 + L/2.0 # mm, radius of nipple plus half the magnet

    if return_magnets:
        return (M, D, L)

    for i in range(len(zlocs)):
        if i in range(6,9):
            sources['{}tr'.format(i)] = Cylmag(loc=[xoffset,ylocs[i],zlocs[i]+zoffset[i]],ori=oris[i], d = D[i], l = L[i], m = M[i])
            sources['{}br'.format(i)] = Cylmag(loc=[-xoffset,ylocs[i],zlocs[i]+zoffset[i]],ori=oris[i], d = D[i], l = L[i], m = M[i])
            sources['{}tl'.format(i)] = Cylmag(loc=[xoffset,ylocs[i],-zlocs[i]-zoffset[i]],ori=oris[i], d = D[i], l = L[i], m = M[i])
            sources['{}bl'.format(i)] = Cylmag(loc=[-xoffset,ylocs[i],-zlocs[i]-zoffset[i]],ori=oris[i], d = D[i], l = L[i], m = M[i])
        else:
            sources['{}r'.format(i)] = Cylmag(loc=[0,ylocs[i],zlocs[i]+zoffset[i]],ori=oris[i], d = D[i], l = L[i], m = M[i])
            sources['{}l'.format(i)] = Cylmag(loc=[0,ylocs[i],-zlocs[i]-zoffset[i]],ori=oris[i], d = D[i], l = L[i], m = M[i])

    C = magpy.Collection()

    for sc in sources:
        C.addSources(sources[sc])

    Y = [[0,y,0] for y in ys]
    BZy = C.getBsweep(Y) # mT

    return BZy[:, 2]*10.0


def get_grad(ylocs, zpos, xoffset, oris, dz = 0.1):
    dBdz = np.zeros(len(zpos))
    for k in range(len(zpos)):

        hlp = np.copy(zpos)
        B1 = get_bfield(ylocs, hlp, ylocs, xoffset, oris)

        hlp[k] += dz
        B2 = get_bfield(ylocs, hlp, ylocs, xoffset, oris)
    
        dBdz[k] = B2[k] - B1[k]

    return dBdz


def do_iter(Bideal, ylocs, zlocs, xoffset, oris, shift_z, dz = 0.1):
    B = get_bfield(ylocs, zlocs, ylocs, xoffset, oris)
    dBdz = get_grad(ylocs, zlocs, xoffset, oris, dz = dz)
    diffs = Bideal - B

    delta = np.sign(diffs) * dBdz * shift_z

    Bnew = get_bfield(ylocs, zlocs + delta, ylocs, xoffset, oris)
    new_diffs = Bideal - Bnew

    # check if we need to decrease the step size
    for k in range(len(zlocs)):
        if np.abs(diffs[k]) < np.abs(new_diffs[k]):            
            print('Reducing step size')
            shift_z[k] *= 0.5

    delta_final = np.sign(diffs) * dBdz * shift_z
    B_final = get_bfield(ylocs, zlocs + delta_final, ylocs, xoffset, oris)
    dBdz_final = get_grad(ylocs, zlocs + delta_final, xoffset, oris, dz = dz)

    return (delta, shift_z, B_final, dBdz_final)


def get_ideal_field(ys, B0, Bg, L0):

    plot_ys = np.copy(ys)
    ind = np.where(plot_ys >= L0*10.0)[0]
    if len(ind) > 0:
        plot_ys[ind] = L0*10.0

    Bfield = B0*(1-(plot_ys)/(L0*10.0))**(0.5) - Bg

    return Bfield


