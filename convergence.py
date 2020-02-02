#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 27 17:18:52 2018

@author: aagrawal
"""

from scipy.integrate import quad

def kmin(zs, chif, hf):
    chis = chif(zs)
    integrand = lambda z : (1+z)*chif(z)*(chis-chif(z))/(chis*hf(z))
    return quad(integrand, 0, zs)[0]
    
def kvar(zs, pkf, chif, hf):
    chis = chif(zs)
    integrand = lambda z : ((1+z)*chif(z)*(chis-chif(z))/chis)**2*pkf(z)/hf(z)
    return quad(integrand, 0, zs)[0]

def deltat(zs, pkf, chif, hf):
    chis = chif(zs)
    integrand = lambda z : (chif(z)*(1-chif(z)/chis))**2*pkf(z)*hf(z)
    return quad(integrand, 0, zs)[0]