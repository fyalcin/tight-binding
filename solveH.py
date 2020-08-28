#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import scipy.linalg as la
from tools import intMatrix, GetPeierlsPhase, fluxtubes
from lattices import IntMatGrapAO
from constants import *
from tools import dot
#import line_profiler
import time

__all__ = ["solveSecularSO", "solveSecularAO"]
#@profile
def solveSecularSO(lattice, kpath, orthogonal = True, alpha=0, Bcenter = np.array([0,0,0]), sigma = 0.1):
    """ For each kpoint in kpath, creates Hamiltonian and overlap matrix elements
    using the 4x4 interaction matrix and the neighbor list for each atom in lattice """
    s=time.time()
    Natoms = len(lattice.atoms)
    bands = np.empty((0,Natoms), float)
    for kpoint in kpath:
        H = np.zeros((Natoms,Natoms), dtype = complex)
        S = np.zeros((Natoms,Natoms), dtype = complex)
        for atom in lattice.atoms:
            for neighbor in atom.neighbors:
#                blochPhase = np.exp(np.complex(0, np.dot(kpoint,neighbor.sepVec)))
                blochPhase = np.exp(np.complex(0, np.dot(kpoint,neighbor.Rvec-atom.Rvec)))
                if lattice.Lattice is "custom":
                    PP = 1
                else:
                    PP = GetPeierlsPhase(lattice, atom, neighbor, alpha)
#                PP = GetPeierlsPhase(lattice, atom, neighbor, alpha)
                H[atom.MUCindex[1]][neighbor.MUCindex[1]] += Vpi*blochPhase*PP/(latConst**2)
                S[atom.MUCindex[1]][neighbor.MUCindex[1]] += Spi*blochPhase*PP/(latConst**2)
        S += np.eye(Natoms)
        H += -V2s*np.eye(Natoms)
        if lattice.Lattice is "custom" and alpha != 0:
            tubes,indices = fluxtubes(lattice,sigma,alpha,Bcenter)
            for index in indices:
                PeierlsFactor = 2*np.pi*(index[2])*np.sign(neighbor.sepVec[0])
                H[index[0]][index[1]] *= np.exp(np.complex(0, PeierlsFactor))
                H[index[1]][index[0]] *= np.exp(np.complex(0, -PeierlsFactor))
        if orthogonal:
            S = np.identity(Natoms)
            energies, eigvecs = la.eigh(H,S)
        else:
            energies, eigvecs = la.eigh(H,S)
#        energies = np.sort(np.real(energies))
        energies = np.real(energies)
        bands = np.vstack([bands, energies])
    e=time.time()
    print("Eigendecomposition took",e-s,"seconds")
    return H, S, bands, eigvecs

def solveSecularAO(lattice, kpath, alpha=0):
    """ For each kpoint in kpath, creates Hamiltonian and overlap matrix elements
    using the 4x4 interaction matrix and the neighbor list for each atom in lattice """
    Natoms = len(lattice.atoms)
    bands = np.empty((0,4*Natoms), float)
    for kpoint in kpath:
        H = np.zeros((4*Natoms,4*Natoms), dtype = complex)
        S = np.zeros((4*Natoms,4*Natoms), dtype = complex)
        for atom in lattice.atoms:
            for neighbor in atom.neighbors:
                blochPhase = np.exp(np.complex(0, np.dot(kpoint,neighbor.Rvec-atom.Rvec)))
                PP = GetPeierlsPhase(lattice, atom, neighbor, alpha)
                hm, sm = intMatrix(neighbor.sepAngle, neighbor.sepDist, orbitals = "all")
                H[4*atom.MUCindex[1]:4*atom.MUCindex[1]+4, 4*neighbor.MUCindex[1]:4*neighbor.MUCindex[1]+4] += hm*blochPhase*PP
                S[4*atom.MUCindex[1]:4*atom.MUCindex[1]+4, 4*neighbor.MUCindex[1]:4*neighbor.MUCindex[1]+4] += sm*blochPhase*PP
        S += np.eye(4*Natoms)
        H += np.diag(np.tile([V2s, V2p, V2p, V2p], Natoms))
        energies, eigvecs = la.eigh(H,S)
        energies = np.sort(np.real(energies))
        bands = np.vstack([bands, energies])
    return H, S, bands, eigvecs

#def solveSecularAO(lattice, kpath, alpha=0):
#    """ For each kpoint in kpath, creates Hamiltonian and overlap matrix elements
#    using the 4x4 interaction matrix and the neighbor list for each atom in lattice """
#    Natoms = len(lattice.atoms)
#    bands = np.empty((0,4*Natoms), float)
#    for kpoint in kpath:
#        H = np.zeros((4*Natoms,4*Natoms), dtype = complex)
#        S = np.zeros((4*Natoms,4*Natoms), dtype = complex)
#        for atom in lattice.atoms:
#            for neighbor in atom.neighbors:
#                blochPhase = np.exp(np.complex(0, np.dot(kpoint,neighbor.Rvec-atom.Rvec)))
#                PP = GetPeierlsPhase(lattice, atom, neighbor, alpha)
#                index = int(np.round(neighbor.sepAngle/(np.pi/3)))
#                H[4*atom.MUCindex[1]:4*atom.MUCindex[1]+4, 4*neighbor.MUCindex[1]:4*neighbor.MUCindex[1]+4] += IntMatGrapAO[index][0]*blochPhase*PP
#                S[4*atom.MUCindex[1]:4*atom.MUCindex[1]+4, 4*neighbor.MUCindex[1]:4*neighbor.MUCindex[1]+4] += IntMatGrapAO[index][1]*blochPhase*PP
#        S += np.eye(4*Natoms)
#        H += np.diag(np.tile([V2s, V2p, V2p, V2p], Natoms))
#        energies, eigvecs = la.eigh(H,S)
#        energies = np.sort(np.real(energies))
#        bands = np.vstack([bands, energies])
##        print(sum(sum(H**2 - (H.transpose().conjugate())**2)))
#    return H, bands, eigvecs
    
#def solveSecular(lattice, kpath, orbitals = "s", alpha=0):
#    """ For each kpoint in kpath, creates Hamiltonian and overlap matrix elements
#    using the 4x4 interaction matrix and the neighbor list for each atom in lattice """
#    M = 4-3*int(orbitals == "s")
#    N = M*len(lattice.atoms)
#    bands = np.empty((0,N), float)
#    for kpoint in kpath:
#        H = np.zeros((N,N), dtype = complex)
#        S = np.zeros((N,N), dtype = complex)
#        for atom in lattice.atoms:
#            for neighbor in atom.neighbors:
#                blochPhase = np.exp(np.complex(0, np.dot(kpoint,neighbor.Rvec-atom.Rvec)))
#                PP = GetPeierlsPhase(lattice, atom, neighbor, alpha)
#                hm, sm = intMatrix(neighbor.sepAngle, neighbor.sepDist, orbitals)
#                H[M*atom.MUCindex[1]:M*atom.MUCindex[1]+M, M*neighbor.MUCindex[1]:M*neighbor.MUCindex[1]+M] += hm*blochPhase*PP
#                S[M*atom.MUCindex[1]:M*atom.MUCindex[1]+M, M*neighbor.MUCindex[1]:M*neighbor.MUCindex[1]+M] += sm*blochPhase*PP
#        S += np.eye(N)
#        if orbitals == "all":
#            H += np.diag(np.tile([V2s, V2p, V2p, V2p], int(N/4)))
#        else:
#            H += -V2s*np.eye(N)
#        energies, eigvecs = la.eigh(H,S)
#        energies = np.sort(np.real(energies))
#        bands = np.vstack([bands, energies])
##        print(sum(sum(H**2 - (H.transpose().conjugate())**2)))
#    return H, bands, eigvecs