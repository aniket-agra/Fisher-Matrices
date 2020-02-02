#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 16:08:08 2019

@author: aagrawal
"""
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

fig_width_pt = 510                     # Suited for LaTeX
inches_per_pt = 1.0/72.27              # Convert pt to inch
golden_mean = (math.sqrt(5)-1.0)/2.0   # Aesthetic ratio
fig_width = fig_width_pt*inches_per_pt # width in inches
fig_height = fig_width*golden_mean     # height in inches
fig_size =  [fig_width,fig_height]
params = {'backend': 'ps', 'axes.linewidth' : 2, 'legend.fontsize': 14, 'text.usetex' : False, \
          'figure.figsize' : fig_size, 'font.family' : 'mono'}
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True

plt.rcParams.update(params)
nplot = 3
indx_c = np.array([0.267, 0.06, -1, 0.0])
indx_l = np.array([0.1, 0.0, -2, -2.7])
indx_u = np.array([0.5, 1, -0.9, 2.7])
indx_lab = np.array([r'$\Omega_{CDM}$', r'$m_{\nu}$'+' (eV)', r'$w_0$', r'$w_a$']) #, 

fig, axs = plt.subplots(figsize = (7, 7), sharex = 'col', sharey = 'row', \
                        ncols = nplot, nrows = nplot)

plt.subplots_adjust(wspace = 0., hspace = 0.)
for i in np.arange(0, nplot) :
    for j in np.arange(0, nplot) : 
        if i < j :
            axs[i, j].axis('off')
        else : 
            submat = np.ix_([j, i+1], [j, i+1])
            C1ii = np.linalg.inv(C3i[submat])
            C2ii = np.linalg.inv(C4i[submat])
#            C1ii = fish3[submat]
#            C2ii = fish2[submat]
            a = C1ii[0, 0]#
            b = 2*C1ii[0, 1]#A[1]
            c = C1ii[1, 1]#A[2]
            aux = np.sqrt((a-c)**2+b**2)
            aux2 = b**2-4*a*c
            ang = np.arctan((c-a-aux)/b)*180/np.pi
            w = -1/aux2*np.sqrt(-2*aux2*(a+c+aux))
            h = -1/aux2*np.sqrt(-2*aux2*(a+c-aux))
            ell1 = Ellipse((indx_c[j], indx_c[i+1]), w*1.52, h*1.52, ang, color = 'r', lw = 2., fill = False, linestyle = '--')
            ell2 = Ellipse((indx_c[j], indx_c[i+1]), w*2.48, h*2.48, ang, color = 'r', lw = 2., fill = False, linestyle = '--')
            a = C2ii[0, 0]#
            b = 2*C2ii[0, 1]#A[1]
            c = C2ii[1, 1]#A[2]
            aux = np.sqrt((a-c)**2+b**2)
            aux2 = b**2-4*a*c
            ang = np.arctan((c-a-aux)/b)*180/np.pi
            w = -1/aux2*np.sqrt(-2*aux2*(a+c+aux))
            h = -1/aux2*np.sqrt(-2*aux2*(a+c-aux))
            ell3 = Ellipse((indx_c[j], indx_c[i+1]), w*1.52, h*1.52, ang, alpha = 1)
            ell4 = Ellipse((indx_c[j], indx_c[i+1]), w*2.48, h*2.48, ang, alpha = 0.40) #linestyle = '--',                         
            axs[i, j].add_artist(ell3)
            axs[i, j].add_artist(ell4)
            axs[i, j].add_artist(ell1)
            axs[i, j].add_artist(ell2)
            axs[i, j].set_xlim(indx_l[j], indx_u[j])
            axs[i, j].xaxis.set_major_locator(plt.MaxNLocator(3))
            axs[i, j].yaxis.set_major_locator(plt.MaxNLocator(3))
            axs[i, j].xaxis.set_minor_locator(AutoMinorLocator())
            axs[i, j].yaxis.set_minor_locator(AutoMinorLocator())
            axs[i, j].tick_params(which='minor', length=2)
            axs[i, j].tick_params(which='major', length=5)
            axs[i, j].set_ylim(indx_l[i+1], indx_u[i+1])
            if j == 0 : 
                axs[i, j].set_ylabel(indx_lab[i+1], fontsize = 14)
                axs[i, j].tick_params(labelsize = 10)
            if i == nplot-1 :
                axs[i, j].set_xlabel(indx_lab[j], fontsize = 14)
#            if j > 0 : 
#                axs[i, j].yaxis.set_visible(False)
            axs[i, j].tick_params(labelsize = 12)
#axs[-1, -1].set_yticks([-1.03, -1.00, -0.97])
plt.subplots_adjust(left=0.01, bottom=0.01)
#fig.tight_layout()
#plt.savefig('plot_4pars_wwom_ztf_new.pdf', format = 'pdf')
plt.show()