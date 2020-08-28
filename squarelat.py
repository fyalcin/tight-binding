#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 18:26:29 2018

@author: firaty
"""

import numpy as np
#import matplotlib.pyplot as plt
from pylab import *
from mpl_toolkits.mplot3d import Axes3D
import copy
from tools import *
from solveH import *
from objects import *

def energysq(t,kx,ky):
    return 4 + 2*t*(np.cos(kx) + np.cos(ky))

def energyhex(t,kx,ky):
    fk = 2*np.cos(np.sqrt(3)*ky*latConst) + 4*np.cos(np.sqrt(3)*ky*latConst/2)*np.cos(3*latConst*kx/2)
    return -t*np.sqrt(3+fk), +t*np.sqrt(3+fk)

latConst = 1.42

a = kpointsBZ(60)

#xcent = 2*np.pi/(3*latConst)
#ycent = 2*np.pi/(3*np.sqrt(3)*latConst)
#ax = np.linspace(xcent - 0.05,xcent + 0.05, 51)
#ay = np.linspace(ycent - 0.05,ycent + 0.05, 51)
#
#X, Y = meshgrid(ax,ay)
#Z1, Z2 = energyhex(-2.7,X,Y)

x = np.linspace(-2*np.pi, 2*np.pi, 300)
y = np.linspace(-2*np.pi, 2*np.pi, 300)
#xcent = 2*np.pi/(3*latConst)
#ycent = 2*np.pi/(3*np.sqrt(3)*latConst)
#x = np.linspace(-2*np.pi/(3*latConst), 2*np.pi/(3*latConst), 100)
#y = np.linspace(-2*np.pi/(3*latConst), 2*np.pi/(3*latConst), 100)
#x = np.linspace(xcent-0.1*xcent, xcent+0.1*xcent, 250)
#y = np.linspace(ycent-0.1*ycent, ycent+0.1*ycent, 250)
X, Y = np.meshgrid(x,y)
Z = energysq(-1,X,Y)


#cell = SimulationCell(1,1,Lattice="PG")

#H,S,bands,eigvecs = solveSecularSO(cell, kpath = a)

#bandz1 = np.asarray([item[0] for item in bands])
#bandz2 = np.asarray([item[1] for item in bands])

#bandzzz = meshgrid(bandz1,bandz2)
#fig2, ax2 = subplots()
#
#im = imshow(Z, cmap=cm.RdBu)
#im.set_interpolation('bilinear')
#
#cb = fig2.colorbar(im)

#fig1, ax1 = subplots()
#
#cnt = contour(X,Y, Z, cmap=cm.RdBu, extent=[0, 1, 0, 1])

fig = plt.figure(figsize=(8,6))

#ax = Axes3D(fig)
ax = fig.add_subplot(111, projection='3d')

#ax = fig.add_subplot(1,1,1, projection='3d')

ax.plot_surface(X, Y, Z, cmap="RdBu", rstride=1, cstride=1, alpha=1)
#ax.plot_surface(X, Y, Z2, cmap=cm.coolwarm, rstride=1, cstride=1, alpha=1)
#ax.view_init(0,60)
#ax.axis('off')
#ax.view_init(0,60)
#cset = ax.contour(X, Y, Z, zdir='z', offset=-pi, cmap=cm.coolwarm)
#cset = ax.contour(X, Y, Z, zdir='x', offset=-pi, cmap=cm.coolwarm)
#cset = ax.contour(X, Y, Z, zdir='y', offset=3*pi, cmap=cm.coolwarm)
#
#ax.set_xlim3d(min(x), max(x));
#ax.set_ylim3d(min(y), max(y));
#ax.set_zlim3d(-8, 8);
#
##cb = fig.colorbar()
#
#norm = matplotlib.colors.Normalize(vmin=-3, vmax=3)
#m = cm.ScalarMappable(cmap=cm.coolwarm, norm=norm)
#m.set_array([])
#plt.colorbar(m, fraction=0.030,pad=0.005)
#b1 = [item[0] for item in bands]
#b2 = [item[1] for item in bands]
#a1 = [item/latConst for item in a1]
#a2 = [item/latConst for item in a2]
##t=np.linspace(0,30,5640)
#ax.scatter(a1,a2,b1)
#ax.scatter(a1,a2,b2)
#    
#plt.show()
##ax.view_init(0,60)
#
ax.set_xlabel("$k_x(1/a)$")
ax.set_ylabel("$k_y(1/a)$")
ax.set_zlabel("$E(eV)$")
#plt.savefig('sqlat-bz.png', format='png', bbox_inches = 'tight', transparent =True)