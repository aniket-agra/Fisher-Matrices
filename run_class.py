#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  7 15:48:50 2019

@author: aagrawal
"""
import numpy as np
from classy import Class
import pdf_and_aux as pdf
import convergence as conv
from scipy import optimize, interpolate

zmax = 4
nz = 100
nk = 1000
sig = 0.12
Mc = 1e11
norm = (5/np.log(10))**2

def constraints(params, zs):
    cosmo = Class()
    cosmo.set(params)
    cosmo.compute()
    h = cosmo.Hubble(0)*299792.458
    om0 = cosmo.Omega0_m()
    #print(cosmo.pars)
    zarr2 = cosmo.get_background().get('z')
    hz = cosmo.get_background().get('H [1/Mpc]')
    hf = interpolate.InterpolatedUnivariateSpline(zarr2[-1:0:-1], hz[-1:0:-1])
    chiz = cosmo.get_background().get('comov. dist.')
    chif = interpolate.InterpolatedUnivariateSpline(zarr2[-1:0:-1], chiz[-1:0:-1])
    
    pkz = np.zeros((nk,nz))
    pklz2 = np.zeros((nk,nz))
    dprime = np.zeros((nk,1))
    zarr = np.linspace(0, zmax, nz)
    karr = np.logspace(-3, np.log10(20), nk)
    rcollarr = np.zeros(nz)
    kcarr = np.zeros(nz)
    delz = 0.01
    
    for i in np.arange(nz):
        Dz = (cosmo.scale_independent_growth_factor(zarr[i])/cosmo.scale_independent_growth_factor(0))
        sigz = lambda x : Dz*cosmo.sigma(x, zarr[i])-1.192182033080519
        if (sigz(1e-5) > 0) & (sigz(10) < 0) : 
            rcollarr[i] = optimize.brentq(sigz, 1e-5, 10)
        else :
            rcollarr[i] = 0
        for j in np.arange(nk):
            pkz[j, i] = cosmo.pk(karr[j], zarr[i])
    for i in np.arange(nk) : 
        pklz0 = np.log(cosmo.pk(karr[i], zs-delz)/cosmo.pk(karr[i], 0))
        pklz1 = np.log(cosmo.pk(karr[i], zs+delz)/cosmo.pk(karr[i], 0))
        pklz2[i] = cosmo.pk(karr[i], 0)
        dprime[i] = -hf(zs)*np.sqrt(cosmo.pk(karr[i],zs)/pklz2[i,0])*(pklz1-pklz0)/4/delz 
        #divided by 2 for step size, another for defining D'   
    
    w0 = params.get('w0_fld')
    wa = params.get('wa_fld')
    mt = 5*np.log10(cosmo.luminosity_distance(zs))  
    #mt = 5*np.log10(fanal.dlatz(zs, om0, og0, w0, wa))  
    Rc = (2*4.302e-9*Mc/h**2/om0)**(1/3)
    mask = (0 < rcollarr) & (rcollarr < Rc)
    kcarr[mask] = 2*np.pi/rcollarr[mask]
    mask = (rcollarr >= Rc)
    kcarr[mask] = 2*np.pi/Rc
    #plt.semilogy(zarr, kcarr)
    pksmooth = pdf.pkint(karr, zarr, pkz, kcarr)
    
    par2 = {'w0' : w0, 'wa' : wa, 'Omega_m' : om0}  
    print(par2)
    #kmin = conv.kmin(zs, chif, hf)*(-3./2.*hf(0)**2*om0)
    kvar = conv.kvar(zs, pksmooth, chif, hf)*(3./2.*hf(0)**2*om0)**2
    #sigln = np.log(1+kvar/np.abs(kmin)**2)
    #L = pdf.convpdf(kmin, sigln, sig, mfid-mt)
    #sigln = np.sqrt(sig**2+(5/np.log(10))**2*kvar)
    #L = pdf.gausspdf(sig, mfid-mt)
    #lnL = mt/sig
    vvar = np.trapz(pklz2[:,0]*dprime[:,0]**2, karr)/6/np.pi**2*(1-(1+zs)/hf(zs)/chif(zs))**2
    #var_tot = norm*kvar+sig**2
    lnL = -mt**2/2
    print('Sigmasq = {}, {}, Likelihood = {}'.format(kvar, vvar, lnL))
    cosmo.struct_cleanup()                                                             
    cosmo.empty()                                                                      
    return [mt, kvar, vvar]
#[mt, var_tot]