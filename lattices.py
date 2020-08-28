#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from constants import *
from tools import intMatrix

__all__ = ["getlatinfo"]

Lattices = {}        
Lattices["square"] = {"basisVec":np.array([[0, 0, 0]]), \
                    "latVec":latConst*np.array([[1, 0, 0],[0, 1, 0]])}
Lattices["PG"] = {"basisVec":np.array([[0.0, 0.0, 0.0], [latConst, 0.0, 0.0]]), \
                    "latVec":latConst*0.5*np.array([[3,np.sqrt(3),0.0],[3,-np.sqrt(3),0.0]])}
Lattices["AG"] = {"basisVec":np.array([[0, 0, 0], [latConst/2, latConst*np.sqrt(3)/2, 0.0],\
                    [0.0, latConst*np.sqrt(3), 0.0], [latConst/2, 3*latConst*np.sqrt(3)/2, 0.0]]),\
                    "latVec":latConst*0.5*np.array([[3, np.sqrt(3), 0.0], [0.0, 4*np.sqrt(3), 0.0]])}
Lattices["PAG"] = {"basisVec":np.array([[0.0, 0.0, 0.0], [latConst, 0.0, 0.0]]),\
                    "latVec":latConst*0.5*np.array([[3,np.sqrt(3),0.0],[0,-2*np.sqrt(3),0.0]])}
Lattices["flake"] = {"basisVec":np.array([[0.0, 0.0, 0.0], [latConst, 0.0, 0.0]]), \
                    "latVec":latConst*0.5*np.array([[3,np.sqrt(3),0.0],[3,-np.sqrt(3),0.0]])}
Lattices["armchair"] = {"basisVec":np.array([[0.0, 0.0, 0.0], [latConst, 0.0, 0.0],[3*latConst/2,\
                        latConst*np.sqrt(3)/2,0],[5*latConst/2,latConst*np.sqrt(3)/2,0]]),\
                        "latVec":latConst*np.array([[3,0,0],[0,np.sqrt(3),0]])}

def getlatinfo(lattice):
    if lattice in Lattices.keys():
        return Lattices[lattice]["basisVec"], Lattices[lattice]["latVec"]
    raise Exception('Please input a correct lattice type')

IntMatGrapS = [Vpi, Spi]

IntMatGrapAO = []
for i in range(6):
    IntMatGrapAO.append(intMatrix(np.pi*i/3, latConst, orbitals="all"))
    