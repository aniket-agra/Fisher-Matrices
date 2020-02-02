#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 26 08:20:22 2018

@author: aagrawal
"""
import numpy as np
import copy 
import run_class as rc
import time
#import matplotlib.pyplot as plt

keys_params = np.array(['Omega_b', 'Omega_cdm', 'm_ncdm', 'Omega_k', 'w0_fld', 'wa_fld'])   #, 'w0_fld', 'wa_fld'
fid_params = np.array([0.049, 0.267, 0.06, 0.1, -1, 0.0])                     #, -1, 0.0
eps_params = np.array([0.001, 0.0001, 0.01, 0.0001, 0.01, 0.01])                       #, 0.01, 0.01
num_params = len(fid_params)
base_cosmo = dict(zip(keys_params, fid_params))
base_cosmo.update({'output' : 'mPk', 'z_max_pk' : 4, 'P_k_max_1/Mpc' : 40, \
                    'ln10^{10}A_s' : 3.094, 'n_s' : 0.9645, \
                    'Omega_scf' : 0.0, 'H0' : 67.27, 'tau_reio' : 0.079, \
                    'N_ur' : 2.0328, 'N_ncdm' :1, 'Omega_Lambda' : 0.0})#, 'input_verbose' : 1, \
                    #'background_verbose' : 1, 'thermodynamics_verbose' : 1, 'Omega_Lambda' : 0.0,
                    #'non linear' : 'halofit', })
fish = np.zeros((num_params, num_params))

zsarr = np.arange(0.2, 1.8, 0.1)
numsn = np.array([60, 200, 400, 220, 320, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140, 140])
dlarr = np.array([1017.5255, 1610.3278, 2250.4729, 2932.2815, 3650.7709, 4401.6222, 5181.1216, 5986.0891, \
                6813.8064, 7661.9492, 8528.5258, 9411.8244, 10310.3680, 11222.8763, 12148.2346, 13085.4670] )
mfidarr = 5*np.log10(dlarr)+30
start = time.time()
indx = np.array([1, 3])
nvars = 2

for r in np.arange(10, np.size(zsarr), 19) :
    fish1 = np.zeros((num_params, num_params))
    base_cosmo.update({'zs' : zsarr[r]})
    for s in np.linspace(0, 1, 1) : #np.linspace(mfidarr[r]-0.12, mfidarr[r]+0.12, 1)
        base_cosmo.update({'mobs' : mfidarr[r]+s*0.12})
        fish_temp = np.zeros((num_params, num_params))
        for p in np.arange(0, nvars) :                         #np.arange(indx, indx+1)
            i = indx[p]
            for q in np.arange(p, nvars) :                      #np.arange(i, indx+1)
                j = indx[q]
                update_dict1 = copy.deepcopy(base_cosmo)
                update_dict2 = copy.deepcopy(base_cosmo)
                update_dict3 = copy.deepcopy(base_cosmo)
                update_dict4 = copy.deepcopy(base_cosmo)
                if j == i :
                    update_dict1.update({keys_params[i] : fid_params[i]+eps_params[i]})
                    update_dict2.update({keys_params[i] : fid_params[i]-eps_params[i]})
                    fish_temp[i,j] = -1/eps_params[i]**2*(rc.constraints(update_dict1)+rc.constraints(update_dict2) \
                             -2*rc.constraints(update_dict3))
                else :
                    update_dict1.update({keys_params[i] : fid_params[i]+eps_params[i], \
                                         keys_params[j] : fid_params[j]+eps_params[j]})
                    update_dict2.update({keys_params[i] : fid_params[i]-eps_params[i], \
                            keys_params[j] : fid_params[j]-eps_params[j]})
                    update_dict3.update({keys_params[i] : fid_params[i]+eps_params[i], \
                            keys_params[j] : fid_params[j]-eps_params[j]})
                    update_dict4.update({keys_params[i] : fid_params[i]-eps_params[i], \
                            keys_params[j] : fid_params[j]+eps_params[j]})
                    fish_temp[i,j] = -1/(4*eps_params[i]*eps_params[j])*(rc.constraints(update_dict1)+\
                         rc.constraints(update_dict2)-rc.constraints(update_dict3)-rc.constraints(update_dict4))
        #print(fish_temp)    
        fish1 = fish1+fish_temp
    fish = fish+fish1*numsn[r]

end = time.time()
print(end-start)