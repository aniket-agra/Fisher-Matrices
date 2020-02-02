#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 17:12:28 2018
LN PDF and LN convolved with Gaussian PDF. Auxiliary code to "smooth" P(k) for 
use in variance of convergence
@author: aagrawal
"""
import numpy as np
from scipy import interpolate
from scipy.integrate import quad
import scipy as sp
def lnpdf(k, kmin, sigln):
       kmin = abs(kmin)
       if k+kmin > 1e-9:
           dpdk = 1/np.sqrt(2*np.pi*sigln)*np.exp(-1/2/sigln* \
                           (np.log(1+k/kmin)+sigln/2)**2)/(k+kmin)
       else:
           dpdk = 0
       return dpdk

def convpdf(kmin, signl, sig, xtot):
        k_xlens = lambda x : 1-np.exp(np.log(10)/5*(xtot-x))
        I = lambda x : (1-k_xlens(x))*lnpdf(k_xlens(x), kmin, signl)* \
                                    np.exp(-x**2/2/sig**2)/np.sqrt(2*np.pi)/sig
        return quad(I, -10*sig, 10*sig)[0]

def gausspdf(sig, x):
       dpdk = np.exp(-x**2/2/sig**2)/np.sqrt(2*np.pi)/sig
       return dpdk
    
def pkint(karr, zarr, pkarr, kcarr):
    pksmooth = np.zeros(np.size(zarr))
    for i in np.arange(np.size(zarr)):
        if kcarr[i] > 0 :
            pksmooth[i] = np.trapz(pkarr[:,i]*karr*np.exp(-(karr/kcarr[i])**2), karr) 
        else : 
            pksmooth[i] =  np.trapz(pkarr[:,i]*karr, karr)
    pksmooth = pksmooth/2/np.pi
    pks = interpolate.InterpolatedUnivariateSpline(zarr, pksmooth)
    return pks

def pkint2(karr, zarr, pkarr, darr, kcarr):
    pksmooth = np.zeros(np.size(zarr))
    for i in np.arange(np.size(zarr)):
        if kcarr[i] > 0 :
            pksmooth[i] = np.trapz(karr*pkarr[:,i]*(darr[:,i]*(1+zarr[i])-1)**2*np.exp(-(karr/kcarr[i])**2), karr) 
        else : 
            pksmooth[i] =  np.trapz(karr*pkarr[:,i]*(darr[:,i]*(1+zarr[i])-1)**2, karr)
    pksmooth = pksmooth/2/np.pi
    pks = interpolate.InterpolatedUnivariateSpline(zarr, pksmooth)
    return pks

def covmatgen(karr, zarr, pkarr, hfs, chis, r, z1, z2) : 
    zarrsize = np.size(zarr)
    covmat = np.zeros((zarrsize, zarrsize))
    karr2 = karr**2
    j0kr = sp.special.spherical_jn(0, karr*r)
    for i in np.arange(0, zarrsize) :
        for j in np.arange(i, zarrsize) : 
            pks = np.trapz(np.sqrt(pkarr[:, i]*pkarr[:, j])*karr2*j0kr, karr)/(2*np.pi**2) 
            covmat[i,j] = chis(zarr[i])*chis(zarr[j])*(1+zarr[i])*(1+zarr[j])*pks/hfs(zarr[i])/hfs(zarr[j])
            if i != j :
                covmat[j,i] = covmat[i,j]
    return covmat