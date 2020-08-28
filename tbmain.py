#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import copy
import scipy.linalg as la
#from mpl_toolkits.mplot3d import Axes3D
from pandas import DataFrame
#matplotlib.interactive(True)
import time
from matplotlib.patches import RegularPolygon
from lattices import *
from tools import *
from objects import SimulationCell  
from solveH import *
from constants import *
from plotcell import plotcell
import datetime
import pylab
from mpl_toolkits.mplot3d import Axes3D
#from concurrent.futures import ProcessPoolExecutor, as_completed
#import multiprocessing as mp
#from sisl import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

gamma = [np.array([0,0,0])]

def butterfly(latticetype,N1,N2,qfactor,orbitals):
    SimCell = SimulationCell(N1, N2, Lattice = latticetype, q = qfactor)
    totBands = []
    Htot = []
    eigfuncs = []
    Stot = []
    for i in range(0,SimCell.q+1):
        if orbitals == "s":
            H,S,band,eigs = solveSecularSO(SimCell, gamma, alpha = i/SimCell.q)
        else:
            H,S,band,eigs = solveSecularAO(SimCell, gamma, alpha = i/SimCell.q)
        totBands.append(band)
        Htot.append(H)
        Stot.append(S)
        eigfuncs.append(eigs)
        print(i)
    for i in range(len(totBands[0][0])-1):
        xx = [item/SimCell.q for item in list(range(SimCell.q+1))]
        plt.plot(xx,[bands[0][i] for bands in totBands],'r.',markersize='1')
    plt.xlabel("p/q (unitless, magnetic field strength)")
    plt.ylabel("Energy (eV)")
    plt.xlim((0,1))
    figname = "bands"+str(N1)+"-"+str(N2)+"-"+str(SimCell.q)+".png"
    plt.savefig(figname,bbox_inches='tight', dpi = 300)
    return Htot, totBands, SimCell, eigfuncs



def plotdosandpdos(cell, Emin, Emax, Ebars, Esigma, PDcenter, PDdelta, PDsigma, PDslices, B, \
                   center, Bsigma, latstr ="cell", size = 1):
#    now = datetime.datetime.now()
#    datestr = str(now.hour) + "-" + str(now.minute) + "-" + str(now.second)
    E = np.linspace(Emin,Emax,Ebars)
    PDE = np.linspace(PDcenter-PDdelta,PDcenter+PDdelta,PDslices)
    
    H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True, Bcenter = center)
#    fig1 = plt.figure(1)
    plt.figure(figsize=(12,9))
#    fig1.tight_layout()
    
    plt.subplot(2,3,1)
    dos1 = DOS(E,eigvals[0],Esigma)
    dos1norm = max(dos1)-min(dos1)
    dos1 = [item/dos1norm for item in dos1]
    plt.plot(E,dos1)
    plt.axhline(y=0, linestyle = "dashed", c ="black")
    plt.ylabel("DOS(E)")
#    plt.ylim((0,2))
#    plt.ylabel("without B")
#    plt.yticks([], [])
    label1 = plt.title("without B")
    label1.set_color("red")
    dosbef = dos1
    
    plt.subplot(2,3,4)
    pdos1 = sumPDOSorthogonal(PDE, eigvals[0], eigvecs, S, PDsigma)
    plotPDOS(cell, pdos1, size)
    plt.ylabel("Real space PDOS")
    plt.xticks([], [])
    plt.yticks([], [])
#    plt.xlabel("with B")
    pdosbef = pdos1
    
    """ Magnetic Field On """
    H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = B, orthogonal = True, Bcenter = center, sigma = Bsigma)
    
    plt.subplot(2,3,2)
    dos1 = DOS(E,eigvals[0],Esigma)
    dos1 = [item/dos1norm for item in dos1]
    plt.plot(E, dos1)
    plt.axhline(y=0, linestyle = "dashed", c ="black")
#    plt.ylim((0,2))
#    plt.title("after")
#    plt.ylabel("with B")
    label2 = plt.title("with B")
    label2.set_color("red")
#    plt.yticks([], [])
    dosaft = dos1
    #
    plt.subplot(2,3,5)
    pdos1 = sumPDOSorthogonal(PDE, eigvals[0], eigvecs, S, PDsigma)
    plotPDOS(cell, pdos1, size)
    plt.xticks([], [])
    plt.yticks([], [])
#    plt.title("after")
    pdosaft = pdos1
    
    """ Difference """
    plt.subplot(2,3,3)
    label3 = plt.title("difference")
    label3.set_color("red")
#    plt.yticks([], [])
    dosdiff = [dosaft[i]-dosbef[i] for i in range(len(dosaft))]
    plt.plot(E,dosdiff)
    plt.axhline(y=0, linestyle = "dashed", c ="black")
    print("Max Difference in DOS occurs at: ",np.argmax(dosdiff),E[np.argmax(dosdiff)])
    
    plt.subplot(2,3,6)
    pdosdiff = [pdosaft[i]-pdosbef[i] for i in range(len(pdosaft))]
    plotPDOS(cell,pdosdiff,size)
    plt.xticks([], [])
    plt.yticks([], [])
    
    plt.savefig(latstr+"_"+str(len(cell.atoms))+ "_"+str(np.round(B*394,2))+"_"+str(Bsigma)+"lowres.png" \
                ,format='png', bbox_inches="tight", transparent = True)
    return H,dosbef,dosaft,pdosbef,pdosaft,pdosdiff

#start = time.time()
#Htot, ads, cell, eigfuncs = butterfly('square',1,1,251,orbitals="s")
#

#cell = SimulationCell(Lattice='custom', M = flake(3,"z"))
#plotflux(cell, 0.5, 10, np.array([0,-2.4595,0]))

#cell = SimulationCell(Lattice = "custom", M = flake(35,"z"))
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 2)
#sE = np.linspace(0,0.3,300)
#dos = DOS(sE,eigvals[0],0.005)
#pE = np.linspace(0.03,0.065,30)
#pdos = sumPDOSorthogonal(pE,eigvals[0],eigvecs, S, pE[1]-pE[0])
#plt.figure()
#plotPDOS(cell,pdos,1000)
#
#plt.figure()
#plotDOS(sE,dos)
#
#pE = np.linspace(0.24,0.27,30)
#pdos = sumPDOSorthogonal(pE,eigvals[0],eigvecs, S, pE[1]-pE[0])
#plt.figure()
#plotPDOS(cell,pdos,1000)
#plotcell(cell)

#cell = SimulationCell(Lattice="custom", M = flake(3,"z"))
#plotflux(cell,1e10,1)



''' DOS BY SIGMA '''
#cell = SimulationCell(Lattice="custom", M = flake(15,"a"))
#
#fig = plt.figure(figsize=(12,6))
#fig.subplots_adjust(hspace=0.05, wspace=0.05)
#for i in range(1,9):  
#    a = plt.subplot(2,4,i)
##    if i == 1:
##        AA = flake(3,"z")
##        aa = inset_axes(a,1,1, loc="lower left")
##        aa.scatter(AA[:,0],AA[:,1], c='grey')
##        a.axis('off')
##        aa.axis('off')
##    elif i > 1 and i < 5:
##        a.axis('off')
##    else:
#    plt.xticks([], [])
#    plt.yticks([], [])
#    alpha = 0.05
#    H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = alpha, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 5*i)
#    E = np.linspace(-8,8,2000)
#    dos = DOS(E,eigvals[0],0.1)
#    plotDOS(E,dos)
##    a.text(0.05, 0.9,"N="+str(len(M)),
##    horizontalalignment='left',
##    verticalalignment='center',
##    transform = a.transAxes)
###        pylab.legend(loc='upper left')
#plt.savefig("dosbysigma-15-a.png", transparent = True, bbox_inches = 'tight', dpi = 200)
''' DOS BY SIGMA '''



''' CENTER STATES BY SIGMA FLAKE '''
#cell = SimulationCell(Lattice="custom", M = flake(15,"a"))
#N = len(cell.atoms)
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 2)
#sE = np.linspace(-8.0,8.0,800)
#dosall = DOS(sE,eigvals[0],0.01)
#norm = sum(dosall)
#E = np.linspace(-0.5,0.5,50)
#ref = DOS(E,eigvals[0],0.01)
#ref = [item*N/norm for item in ref]
#ref = sum(ref)
#
#for j in range(1,5):
#    dosbysigma = []
#    B = 0.0250*j
#    sigma = []
#    for i in range(1,31):
#        H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = B, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 2*i)
#        dos = DOS(E,eigvals[0],0.01)
#        dos = [item*N/norm for item in dos]
#        dosbysigma.append(sum(dos)-ref)
#        sigma.append(2*i)
#    lab = str(int(400*B))+"T"
#    pylab.plot(sigma, dosbysigma, label=lab)
#    pylab.xlabel(r'$\sigma(A)$')
#    pylab.ylabel(r'$DOS(E \in [-0.5,0.5]eV)$')
#    pylab.legend(loc='upper left')
#    pylab.xlim((0,60))
#    pylab.ylim((0,80))
#pylab.savefig("centerstatesbysigma-flake15a-magfieldeffect.png", bbox_inches='tight', transparent = True)
''' CENTER STATES BY SIGMA FLAKE '''


''' CENTER STATES BY SIGMA JUNCTION A'''
#cell = SimulationCell(Lattice="custom", M = junction(15,5,"a"))
#N = len(junction(15,5,"a"))
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 2)
#sE = np.linspace(-8.0,8.0,800)
#dosall = DOS(sE,eigvals[0],0.01)
#norm = sum(dosall)
#E = np.linspace(-0.5,0.5,50)
#
#for j in range(1,5):
#    dosbysigma = []
#    B = 0.0250*j
#    sigma = []
#    for i in range(1,30):
#        H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = B, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 2*i)
#        dos = DOS(E,eigvals[0],0.01)
#        dos = [item*N/norm for item in dos]
#        dosbysigma.append(sum(dos))
#        sigma.append(2*i)
#    lab = str(int(400*B))+"T"
#    pylab.plot(sigma, dosbysigma, label=lab)
#    pylab.xlabel(r'$\sigma(A)$')
#    pylab.ylabel(r'$DOS(E \in [-0.5,0.5]eV)$')
#    pylab.legend(loc='upper left')
#    pylab.xlim((0,60))
#    pylab.ylim((0,100))
#pylab.savefig("centerstatesbysigma-junction15_5_a.png", bbox_inches='tight', dpi = 300, transparent = True)
''' CENTER STATES BY SIGMA JUNCTION A'''


''' CENTER STATES BY SIGMA JUNCTION Z'''
#cell = SimulationCell(Lattice="custom", M = junction(17,5,"z"))
#N = len(junction(17,5,"z"))
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 2)
#sE = np.linspace(-8.0,8.0,800)
#dosall = DOS(sE,eigvals[0],0.01)
#norm = sum(dosall)
#E = np.linspace(-0.5,0.5,50)
#
#for j in range(1,5):
#    dosbysigma = []
#    B = 0.0250*j
#    sigma = []
#    for i in range(1,30):
#        H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = B, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 2*i)
#        dos = DOS(E,eigvals[0],0.01)
#        dos = [item*N/norm for item in dos]
#        dosbysigma.append(sum(dos))
#        sigma.append(2*i)
#    lab = str(int(400*B))+"T"
#    pylab.plot(sigma, dosbysigma, label=lab)
#    pylab.xlabel(r'$\sigma(A)$')
#    pylab.ylabel(r'$DOS(E \in [-0.5,0.5]eV)$')
#    pylab.legend(loc='upper left')
#    pylab.xlim((0,60))
#    pylab.ylim((0,100))
#pylab.savefig("centerstatesbysigma-junction17_5_z.png", bbox_inches='tight', dpi = 300, transparent = True)
''' CENTER STATES BY SIGMA JUNCTION A'''




''' CENTER STATES BY SIZE '''

#E = np.linspace(-0.5,0.5,50)
#sE = np.linspace(-8.0,8.0,800)
#for i in range(0,3):
#    B = 0.025*i
#    dosbyB = []
#    size = []
#    for j in range(3,15):
#        cell = SimulationCell(Lattice="custom", M = flake(2*j-3,"z"))
#        N = len(cell.atoms)
#        H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 0)
#        dosall = DOS(sE,eigvals[0],0.01)
#        norm = sum(dosall)
#        H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = B, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 1e10)
#        dos = DOS(E,eigvals[0],0.01)
#        dos = [item*N/norm for item in dos]
#        dosbyB.append(sum(dos))
#        size.append(N)
#    lab = str(int(400*B))+"T"
#    pylab.plot(size, dosbyB, label = lab)
#    pylab.xlabel('Number of atoms')
#    pylab.ylabel(r'$DOS(E \in [-0.5,0.5]eV)$')
#    pylab.legend(loc='upper left')
#pylab.savefig("centerstatesbysize-uniformfield-armchair.png", bbox_inches='tight', dpi = 300, transparent = True)

''' CENTER STATES BY SIZE '''






''' DOS BY SIZE - ZIGZAG FLAKE '''

#fig = plt.figure(figsize=(12,9))
#fig.subplots_adjust(hspace=0.05, wspace=0.05)
#for i in range(1,13):  
#    a = plt.subplot(3,4,i)
#    if i == 1:
#        AA = flake(3,"z")
##        aa = inset_axes(a,1,1, loc="lower left")
##        aa.scatter(AA[:,0],AA[:,1], c='grey')
#        a.axis('off')
##        aa.axis('off')
#    elif i > 1 and i < 5:
#        a.axis('off')
#    else:
#        plt.xticks([], [])
#        plt.yticks([], [])
#        alpha = 0
#        M = flake(6*i-25,"z")
#        cell = SimulationCell(Lattice="custom", M = M)
#        H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = alpha, orthogonal = True, Bcenter=np.array([-5.68,0,0]), sigma = 1)
#        E = np.linspace(-8,8,1000)
#        dos = DOS(E,eigvals[0],0.05)
#        plotDOS(E,dos)
#        a.text(0.05, 0.9,"N="+str(len(M)),
#        horizontalalignment='left',
#        verticalalignment='center',
#        transform = a.transAxes)
##        pylab.legend(loc='upper left')
#    
#plt.savefig("dossizeflakez.ps", format='ps', bbox_inches = 'tight', transparent = True)

''' DOS BY SIZE - ZIGZAG FLAKE '''
    



''' DOS BY SIZE - ARMCHAIR FLAKE '''

#fig = plt.figure(figsize=(12,9))
#fig.subplots_adjust(hspace=0.05, wspace=0.05)
#for i in range(1,13):  
#    a = plt.subplot(3,4,i)
#    if i == 1:
#        AA = flake(3,"a")
#        aa = inset_axes(a,1,1, loc="lower left")
#        aa.scatter(AA[:,0],AA[:,1], c='grey')
#        a.axis('off')
#        aa.axis('off')
#    elif i > 1 and i < 5:
#        a.axis('off')
#    else:
#        plt.xticks([], [])
#        plt.yticks([], [])
#        alpha = 0
#        M = flake(4*i-15,"a")
#        cell = SimulationCell(Lattice="custom", M = M)
#        H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = alpha, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 1)
#        E = np.linspace(-8,8,1000)
#        dos = DOS(E,eigvals[0],0.05)
#        plotDOS(E,dos)
#        a.text(0.05, 0.9,"N="+str(len(M)),
#        horizontalalignment='left',
#        verticalalignment='center',
#        transform = a.transAxes)
##        pylab.legend(loc='upper left')
#    
#plt.savefig("dossizeflakea.png", format='png', bbox_inches = 'tight', transparent = True)

''' DOS BY SIZE - ARMCHAIR FLAKE '''



''' DOS BY SIZE - ZIGZAG TRIANGLE '''

#fig = plt.figure(figsize=(12,9))
#fig.subplots_adjust(hspace=0.05, wspace=0.05)
#for i in range(1,13):  
#    a = plt.subplot(3,4,i)
#    if i == 1:
#        AA = triangle(3,"z")
#        for index,item in enumerate(AA):
#            AA[index] = rotate(item, -np.pi/2)
#        aa = inset_axes(a,1,1, loc="lower left")
#        aa.scatter(AA[:,0],AA[:,1], c='grey')
#        a.axis('off')
#        aa.axis('off')
#    elif i > 1 and i < 5:
#        a.axis('off')
#    else:
#        plt.xticks([], [])
#        plt.yticks([], [])
#        alpha = 0
#        M = triangle(8*i-35,"z")
#        cell = SimulationCell(Lattice="custom", M = M)
#        H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = alpha, orthogonal = True, Bcenter=np.array([-5.68,0,0]), sigma = 1)
#        E = np.linspace(-8,8,1000)
#        dos = DOS(E,eigvals[0],0.1)
#        plotDOS(E,dos)
#        a.text(0.05, 0.9,"N="+str(len(M)),
#        horizontalalignment='left',
#        verticalalignment='center',
#        transform = a.transAxes)
#    
#plt.savefig("dosbysize-triangle-z.png", bbox_inches = 'tight', dpi = 300, transparent = True)

''' DOS BY SIZE - ZIGZAG TRIANGLE '''

#time.sleep(20)

''' DOS BY SIZE - ARMCHAIR TRIANGLE '''
#
#fig = plt.figure(figsize=(12,9))
#fig.subplots_adjust(hspace=0.05, wspace=0.05)
#for i in range(1,13):  
#    a = plt.subplot(3,4,i)
#    if i == 1:
#        AA = triangle(3,"a")
#        aa = inset_axes(a,1,1, loc="lower left")
#        aa.scatter(AA[:,0],AA[:,1], c='grey')
#        a.axis('off')
#        aa.axis('off')
#    elif i > 1 and i < 5:
#        a.axis('off')
#    else:
#        plt.xticks([], [])
#        plt.yticks([], [])
#        alpha = 0
#        M = triangle(4*i-15,"a")
#        cell = SimulationCell(Lattice="custom", M = M)
#        H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = alpha, orthogonal = True, Bcenter=np.array([-5.68,0,0]), sigma = 1)
#        E = np.linspace(-8,8,1000)
#        dos = DOS(E,eigvals[0],0.1)
#        plotDOS(E,dos)
#        a.text(0.05, 0.9,"N="+str(len(M)),
#        horizontalalignment='left',
#        verticalalignment='center',
#        transform = a.transAxes)
#    
#plt.savefig("dosbysize-triangle-a.png", bbox_inches = 'tight', transparent = True)

''' DOS BY SIZE - ARMCHAIR TRIANGLE '''






''' DOS BY SIZE - ARMCHAIR/ZIGZAG JUNCTION '''

#fig = plt.figure(figsize=(12,9))
#fig.subplots_adjust(hspace=0.05, wspace=0.05)
#for i in range(1,13):  
#    a = plt.subplot(3,4,i)
#    if i == 1:
#        AA = junction(5,2,"z")
##        for index,item in enumerate(AA):
##            AA[index] = rotate(item, -np.pi/2)
#        aa = inset_axes(a,1,1, loc="lower left")
#        aa.scatter(AA[:,0],AA[:,1], c='grey',s=10)
#        a.axis('off')
#        aa.axis('off')
#    elif i > 1 and i < 5:
#        a.axis('off')
#    else:
#        plt.xticks([], [])
#        plt.yticks([], [])
#        alpha = 0
#        M = junction(3*i-12,int(4*i/3)-4,"a")
#        cell = SimulationCell(Lattice="custom", M = M)
#        H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = alpha, orthogonal = True, Bcenter=np.array([-5.68,0,0]), sigma = 1)
#        E = np.linspace(-8,8,1000)
#        dos = DOS(E,eigvals[0],0.05)
#        plotDOS(E,dos)
#        a.text(0.05, 0.9,"N="+str(len(M)),
#        horizontalalignment='left',
#        verticalalignment='center',
#        transform = a.transAxes)
#    
#plt.savefig("dossizejunctionz.png", bbox_inches = 'tight', transparent = True)

''' DOS BY SIZE - ARMCHAIR/ZIGZAG JUNCTION '''





#dos = DOS(E,eigvals[0],0.1)
#plotDOS(E,dos)

#dosg = []
#spec = []
#for en in E:
#    G = greensclosed(H,en)
#    spec.append(spectral(G))
#    dosg.append(-np.trace(spectral(G))/(2*np.pi))

#plotDOS(E,dosg)
#plotflux(cell, 3, 1)


#for atom in cell.atoms:
#    print(np.asarray(atom.coord)+np.array([10,10,0]))
#plotcell(cell,"all")
###E = 0.5
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True, Bcenter=np.array([-5.68,0,0]), sigma = 1)
#
#mu = 0
#kT = 1e8
#rho = 1+np.exp((mu*np.eye(len(H))-np.diag(eigvals[0]))/kT)
#rho = np.dot(eigvecs,np.dot(rho,eigvecs.T.conj()))
#rho = np.round(np.real(np.diag(rho)),2)
#minrho = min(rho)
#rangerho = max(rho)-min(rho)
#rho = [(item -minrho)/rangerho for item in rho]

#E = np.linspace(-8,8,100)
#dos = DOS(E,eigvals[0],0.1)
#plotDOS(E,dos)


#cell = triangle(3,"z")
#a = open("mag.txt")
#b = a.read().split("\n")
#del b[-1]
#
#mag = [float(item) for item in b]
#minmag = min(mag)
#mag = [500*(item - minmag) for item in mag]
#x = [atom[0] for atom in cell]
#y = [atom[1] for atom in cell]
#plt.scatter(x,y,c=mag, cmap = "coolwarm", s=mag)
##plt.annotate(str(mag[i]),(j.coord[0],j.coord[1]))



#E = np.linspace(-0.5,0.5,200)
#pdos = sumPDOSorthogonal(E,eigvals[0],eigvecs,S,0.01)
#print(np.round(eigvals[0],2))

#lenatoms = []
#for i in range(1,30):
#    lenatoms.append(len(triangle(2*i+1,"z")))
#plt.plot(list(range(1,30)),lenatoms)
##plotPDOS(cell,pdos,100)
#kT=0.025
#mu= 0.25
##rho = la.logm(1+np.exp((mu*np.eye(len(H))-H)/kT))
##rho = np.diag(rho)/(1.42*1e-9)
#rho = np.log(1+np.exp((mu*np.eye(len(H))-np.diag(eigvals[0]))/kT))
#rho = np.dot(eigvecs,np.dot(rho,eigvecs.T.conj()))
#rho = np.round(np.real(np.diag(rho)),2)
#
#bef = rho
#
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0.5, orthogonal = True, Bcenter=np.array([-5.68,0,0]), sigma = 1)
#kT=0.025
#mu= 0.25
#rho = np.log(1+np.exp((mu*np.eye(len(H))-np.diag(eigvals[0]))/kT))
#rho = np.dot(eigvecs,np.dot(rho,eigvecs.T.conj()))
#rho = np.round(np.real(np.diag(rho)),2) 
#
#aft = rho
#aa =aft-bef
#aa = np.round(aa,2)
#
#for i,j in enumerate(aa):
#    plt.scatter(cell.atoms[i].coord[0],cell.atoms[i].coord[1],c="r")
#    plt.annotate(str(j), (cell.atoms[i].coord[0],cell.atoms[i].coord[1]))
    



''' START : RATIO OF EDGE STATES with B OFF and ON for ZIGZAG type '''
#width = 13
#edgetype = "z"
#latstr = "triangle"+str(width)+edgetype
#
#M = triangle(width,edgetype)
#cell = SimulationCell(Lattice = "custom", M= M)
##plotcell(cell)
#pdosbydos = []
#dE = np.linspace(-8,8,500)
#pE = np.linspace(-0.5,0.5,5)
#for alp in np.linspace(0,1.5,32):
#    H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = alp, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 1e7)
#    dos = DOS(dE,eigvals[0],0.01)
#    pdos = sumPDOSorthogonal(pE,eigvals[0],eigvecs,S,0.01)
#    dosnorm = DOS(pE,eigvals[0],0.01)
#    dosnormsum = sum(dosnorm)
#    pdossum = sum(pdos)
#    pdos = [item*dosnormsum/pdossum for item in pdos]
#    pdosbydos.append(sum(pdos)/sum(dos))
#
#plt.plot(np.linspace(0,1.5,32),pdosbydos)
''' END : RATIO OF EDGE STATES with B OFF and ON for ZIGZAG type '''




''' START : RADIO OF EDGE STATES FOR VARYING LATTICE SIZE with CONSTANT B '''
#ratio0 = []
#ratio1 = []
#dE = np.linspace(-8,8,500)
#pE = np.linspace(-0.5,0.5,5)
#
#for i in range(4,15):
#    print(i)
#    cell = SimulationCell(Lattice="custom", M = triangle(2*i+1,"a"))
#    H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 1)
#    dos = DOS(dE,eigvals[0],0.01)
#    pdos = sumPDOSorthogonal(pE,eigvals[0],eigvecs,S,0.01)
#    dosnorm = DOS(pE,eigvals[0],0.01)
#    dosnormsum = sum(dosnorm)
#    pdossum = sum(pdos)
#    pdos = [item*dosnormsum/pdossum for item in pdos]
#    ratio0.append(sum(pdos)/sum(dos))
#    
#    H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0.05, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = i)
#    dos = DOS(dE,eigvals[0],0.01)
#    pdos = sumPDOSorthogonal(pE,eigvals[0],eigvecs,S,0.01)
#    dosnorm = DOS(pE,eigvals[0],0.01)
#    dosnormsum = sum(dosnorm)
#    pdossum = sum(pdos)
#    pdos = [item*dosnormsum/pdossum for item in pdos]
#    ratio1.append(sum(pdos)/sum(dos))
#
#   
#pylab.plot(list(range(4,15)),ratio0, label = "zero B")
#pylab.plot(list(range(4,15)),ratio1, label = "localized B") 
#pylab.legend(loc='upper right')
##plt.savefig("ratioofedgestates-localizedB.png", dpi=300, bbox_inches="tight", transparent = True)
''' END : RADIO OF EDGE STATES FOR VARYING LATTICE SIZE with CONSTANT B '''










''' START : JUNCTION PDOS/DOS DIFF '''
#length = 15
#width = 7
#edgetype = "z"
#latstr = "junction"+str(length)+"_"+str(width)+edgetype
#
#M = junction(length,width,edgetype)
#
#cell = SimulationCell(Lattice = "custom", M= M)
##plotcell(cell)
#H,dosbef,dosaft,pdosbef,pdosaft,pdosdiff = plotdosandpdos(cell,Emin=-8,Emax=8,Ebars=500,Esigma=0.1,PDcenter=0,\
#            PDdelta=0.50, PDsigma=0.05,PDslices=3, B = 0.2, center=np.array([0,0,0]), Bsigma = 10, latstr = latstr)
''' END : JUNCTION PDOS/DOS DIFF '''



''' START : PDOS by changing B - junction '''
#length = 15
#width = 6
#edgetype = "a"
#latstr = "junction"+str(length)+"_"+str(width)+edgetype
#
#M = junction(length,width,edgetype)
#
#cell = SimulationCell(Lattice = "custom", M= M)
#for ind,alp in enumerate(np.linspace(0,1.5,16)):
#    plt.subplot(4,4,ind+1)
#    H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, orthogonal = True, alpha = alp, sigma = 3)
#    E = np.linspace(-0.5,0.5,3)
#    pdos = sumPDOSorthogonal(E,eigvals[0],eigvecs,S,0.05)
#    plotPDOS(cell,pdos, 10)
#plt.savefig(latstr + "_" +"PDOSbyB_"+ "sigma3.png", dpi = 300)
''' END : PDOS by changing B - junction '''




#plt.figure(3,figsize=(11,1))
#for index, alpha in enumerate(np.linspace(0,1,11)):
#    H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, orthogonal = True, alpha = alpha, sigma = 0.2)
#    plt.subplot(1,11,index+1)
#    e = np.linspace(-0.5,0.5,3)
#    pdos = sumPDOSorthogonal(e,eigvals[0],eigvecs,S,0.05)
#    plotPDOS(cell,pdos,1)
  
#cell = SimulationCell(Lattice = "custom", M = triangle(11,"a"))
#plotcell(cell)


''' START : EDGE EFFECTS STUDIED '''
#for i in range(0,9):
#    size = 2*i+1
#    cell = SimulationCell(Lattice = "custom", M = triangle(size,"z"))
#    H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0.2, orthogonal = True,\
#                                            Bcenter=np.array([0,0,0]), sigma = 2)
#    plt.subplot(3,3,i+1)
#    en = np.linspace(-8,8,500)
#    pdos = sumPDOSorthogonal(en,eigvals[0],eigvecs,S,0.05)
##    dos = DOS(en,eigvals[0],0.1)
#    plotPDOS(cell,pdos)
##    plotDOS(en,dos)
#    plt.savefig("triangle_zigzag_gaussB_pdos.png", dpi=300)
''' END : EDGE EFFECTS STUDIED '''


''' START : BASIC TRANSMISSION PROFILE PLOT  -- ZIGZAG RIBBON '''
cell = SimulationCell(Lattice = "custom", M = ribbon(2,4,"z"))    
cell.remove([32,33,34,35,36,37,38,39])
#cell.remove([26])
#cell.rotate(np.pi/2)
#plotcell(cell)
#
H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 10000)
H = H.conj()
#plotflux(cell,1,0.5,Bcenter=np.array([-8.41,0,0]))

T = []
dosg = []
N = len(H[0])
spect = np.zeros((N,N),dtype=complex)
for en in np.linspace(-8,8,1000):
    G, sigma, subsigma = greensopen(leaddict["ribbon_2_4_z"], H, en)   ### change this line
#    subG = gsubmat("ribbon_4_1_z",0,1,G)
    gamma1, subgamma1 = broadening(sigma[0],subsigma[0])
    gamma2, subgamma2 = broadening(sigma[1],subsigma[1])
    spect += spectral(G)
#    gamma1 = np.complex(0,1)*(sigma[0]-sigma[0].T.conj())
#    gamma2 = np.complex(0,1)*(sigma[1]-sigma[1].T.conj())
#    Tepoint = np.trace(np.dot(subbroad1,np.dot(subG, np.dot(subbroad2,subG.T.conj()))))
    t = transmission(G,en,gamma1,gamma2)
    dosg.append(-np.imag(np.trace(G))/np.pi)
    T.append(t)
plt.plot(np.linspace(-8,8,1000),T) 
plt.ylabel(r'$G(2e^2/h)$')
plt.xlabel(r'$E(eV)$')
plt.title('Conductance vs. E for Zigzag Nanoribbon')
#plt.savefig("zigzagtrans.png", bbox_inches='tight', transparent = True)

#G, sigma, subsigma = greensopen(leaddict["ribbon_2_4_z"], H, 0.02)
#H += sigma[0] + sigma[1]
#kT=0.025
#mu= 0
#rho = np.log(1+np.exp((mu*np.eye(len(H))-np.diag(eigvals[0]))/kT))
#rho = np.dot(eigvecs,np.dot(rho,eigvecs.T.conj()))
#rho = np.round(np.real(np.diag(rho)),2)
#
#plotPDOS(cell,rho,1)
#pdosopen = np.diag(spect)

##plt.plot(np.linspace(-8,8,1000),dosg)
#print(sum(T))
''' END : BASIC TRANSMISSION PROFILE PLOT  -- ZIGZAG RIBBON '''





''' START : COMPARISON OF ARMCHAIR RIBBON WITH ARMCHAIR JUNCTION - TRANSMISSION '''
#cell = SimulationCell(Lattice = "custom", M = junction(4,2,"a"))    
##cell.remove([32,33,34,35,36,37,38,39])
##cell.remove([26])
##cell.rotate(np.pi/2)
#pylab.figure(figsize=(12,9))
##plotcell(cell)
#ax1 = pylab.subplot(2,2,1)
#ax1.axis('off')
#for atom in cell.atoms:
#    ax1.scatter(atom.coord[0],atom.coord[1],c="r", s=3000/len(cell.atoms))
#ax1.set_ylim((-16.25,16.25))
#ax1.set_xlim((-21,21))
#pylab.xticks([],[])
#pylab.yticks([],[])
##
#
#''' MAG FIELD 1 '''
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0.05, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 5)
#
#T1 = []
#dosg = []
#N = len(H[0])
#spect = np.zeros((N,N),dtype=complex)
#for en in np.linspace(-8,8,1000):
#    G, sigma, subsigma = greensopen(leaddict["junction_4_2_a"], H, en)   ### change this line
##    subG = gsubmat("ribbon_4_1_z",0,1,G)
#    gamma1, subgamma1 = broadening(sigma[0],subsigma[0])
#    gamma2, subgamma2 = broadening(sigma[1],subsigma[1])
#    spect += spectral(G)
#    t = transmission(G,en,gamma1,gamma2)
#    dosg.append(-np.imag(np.trace(G))/np.pi)
#    T1.append(t)
#    
#T2 = []
#for en in np.linspace(-8,8,1000):
#    G, sigma, subsigma = greensopen(leaddict["junction_4_2_a"], H, en)   ### change this line
##    subG = gsubmat("ribbon_4_1_z",0,1,G)
#    gamma1, subgamma1 = broadening(sigma[0],subsigma[0])
#    gamma2, subgamma2 = broadening(sigma[2],subsigma[2])
#    spect += spectral(G)
#    t = transmission(G,en,gamma1,gamma2)
#    dosg.append(-np.imag(np.trace(G))/np.pi)
#    T2.append(t)    
#T = [T1[i]+T2[i] for i in range(len(T1))]
#
#pdosopen = np.diag(spect)
#ax2 = pylab.subplot(2,2,2)
#ax2.plot(np.linspace(-8,8,1000),T, label = "B = 20T", linestyle = ':')
#
#''' MAG FIELD 2 '''
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0.1, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 5)
#
#T1 = []
#dosg = []
#N = len(H[0])
#spect = np.zeros((N,N),dtype=complex)
#for en in np.linspace(-8,8,1000):
#    G, sigma, subsigma = greensopen(leaddict["junction_4_2_a"], H, en)   ### change this line
##    subG = gsubmat("ribbon_4_1_z",0,1,G)
#    gamma1, subgamma1 = broadening(sigma[0],subsigma[0])
#    gamma2, subgamma2 = broadening(sigma[1],subsigma[1])
#    spect += spectral(G)
#    t = transmission(G,en,gamma1,gamma2)
#    dosg.append(-np.imag(np.trace(G))/np.pi)
#    T1.append(t)
#    
#T2 = []
#for en in np.linspace(-8,8,1000):
#    G, sigma, subsigma = greensopen(leaddict["junction_4_2_a"], H, en)   ### change this line
##    subG = gsubmat("ribbon_4_1_z",0,1,G)
#    gamma1, subgamma1 = broadening(sigma[0],subsigma[0])
#    gamma2, subgamma2 = broadening(sigma[2],subsigma[2])
#    spect += spectral(G)
#    t = transmission(G,en,gamma1,gamma2)
#    dosg.append(-np.imag(np.trace(G))/np.pi)
#    T2.append(t)    
#T = [T1[i]+T2[i] for i in range(len(T1))]
#
#pdosopen = np.diag(spect)
#ax2.plot(np.linspace(-8,8,1000),T, label = "B = 40T", linestyle = '--')
#
#
#''' ZERO FIELD '''
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True, Bcenter=np.array([0,0,0]), sigma = 5)
#
#T1 = []
#dosg = []
#N = len(H[0])
#spect = np.zeros((N,N),dtype=complex)
#for en in np.linspace(-8,8,1000):
#    G, sigma, subsigma = greensopen(leaddict["junction_4_2_a"], H, en)   ### change this line
##    subG = gsubmat("ribbon_4_1_z",0,1,G)
#    gamma1, subgamma1 = broadening(sigma[0],subsigma[0])
#    gamma2, subgamma2 = broadening(sigma[1],subsigma[1])
#    spect += spectral(G)
#    t = transmission(G,en,gamma1,gamma2)
#    dosg.append(-np.imag(np.trace(G))/np.pi)
#    T1.append(t)
#    
#T2 = []
#for en in np.linspace(-8,8,1000):
#    G, sigma, subsigma = greensopen(leaddict["junction_4_2_a"], H, en)   ### change this line
##    subG = gsubmat("ribbon_4_1_z",0,1,G)
#    gamma1, subgamma1 = broadening(sigma[0],subsigma[0])
#    gamma2, subgamma2 = broadening(sigma[2],subsigma[2])
#    spect += spectral(G)
#    t = transmission(G,en,gamma1,gamma2)
#    dosg.append(-np.imag(np.trace(G))/np.pi)
#    T2.append(t)    
#T = [T1[i]+T2[i] for i in range(len(T1))]
#
#pdosopen = np.diag(spect)
#ax2.plot(np.linspace(-8,8,1000),T, label = "B = 40T")
#
#ax2.legend(loc='upper right')
#
#
#''' COMPARE WITH RIBBON '''
#cell = SimulationCell(Lattice = "custom", M = ribbon(6,2,"a"))    
#ax3 = pylab.subplot(2,2,3)
#ax3.axis('off')
#for atom in cell.atoms:
#    ax3.scatter(atom.coord[0],atom.coord[1],c="r")
#ax3.set_ylim((-7.5,15))
#ax3.set_xlim((-4,26))
#pylab.xticks([],[])
#pylab.yticks([],[])
#xcm = sum(ribbon(6,2,"a")[:,0])/len(ribbon(6,2,"a"))
#ycm = sum(ribbon(6,2,"a")[:,1])/len(ribbon(6,2,"a"))
#
#''' ZERO FIELD '''
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True, Bcenter=np.array([xcm,ycm,0]), sigma = 3)
#
#T = []
#dosg = []
#N = len(H[0])
#spect = np.zeros((N,N),dtype=complex)
#for en in np.linspace(-8,8,1000):
#    G, sigma, subsigma = greensopen(leaddict["ribbon_6_2_a"], H, en)   ### change this line
##    subG = gsubmat("ribbon_4_1_z",0,1,G)
#    gamma1, subgamma1 = broadening(sigma[0],subsigma[0])
#    gamma2, subgamma2 = broadening(sigma[1],subsigma[1])
#    spect += spectral(G)
##    gamma1 = np.complex(0,1)*(sigma[0]-sigma[0].T.conj())
##    gamma2 = np.complex(0,1)*(sigma[1]-sigma[1].T.conj())
##    Tepoint = np.trace(np.dot(subbroad1,np.dot(subG, np.dot(subbroad2,subG.T.conj()))))
#    t = transmission(G,en,gamma1,gamma2)
#    dosg.append(-np.imag(np.trace(G))/np.pi)
#    T.append(t)
#    
#pdosopen = np.diag(spect)
#ax4 = pylab.subplot(2,2,4)
#ax4.plot(np.linspace(-8,8,1000),T, label = "ZERO B")
#
#''' NONZERO FIELD '''
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0.05, orthogonal = True, Bcenter=np.array([xcm,ycm,0]), sigma = 3)
#
#T = []
#dosg = []
#N = len(H[0])
#spect = np.zeros((N,N),dtype=complex)
#for en in np.linspace(-8,8,1000):
#    G, sigma, subsigma = greensopen(leaddict["ribbon_6_2_a"], H, en)   ### change this line
##    subG = gsubmat("ribbon_4_1_z",0,1,G)
#    gamma1, subgamma1 = broadening(sigma[0],subsigma[0])
#    gamma2, subgamma2 = broadening(sigma[1],subsigma[1])
#    spect += spectral(G)
##    gamma1 = np.complex(0,1)*(sigma[0]-sigma[0].T.conj())
##    gamma2 = np.complex(0,1)*(sigma[1]-sigma[1].T.conj())
##    Tepoint = np.trace(np.dot(subbroad1,np.dot(subG, np.dot(subbroad2,subG.T.conj()))))
#    t = transmission(G,en,gamma1,gamma2)
#    dosg.append(-np.imag(np.trace(G))/np.pi)
#    T.append(t)
#    
#pdosopen = np.diag(spect)
#ax4.plot(np.linspace(-8,8,1000),T, label = "B = 10T", linestyle = ':')
#
#
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0.1, orthogonal = True, Bcenter=np.array([xcm,ycm,0]), sigma = 3)
#
#T = []
#dosg = []
#N = len(H[0])
#spect = np.zeros((N,N),dtype=complex)
#for en in np.linspace(-8,8,1000):
#    G, sigma, subsigma = greensopen(leaddict["ribbon_6_2_a"], H, en)   ### change this line
##    subG = gsubmat("ribbon_4_1_z",0,1,G)
#    gamma1, subgamma1 = broadening(sigma[0],subsigma[0])
#    gamma2, subgamma2 = broadening(sigma[1],subsigma[1])
#    spect += spectral(G)
##    gamma1 = np.complex(0,1)*(sigma[0]-sigma[0].T.conj())
##    gamma2 = np.complex(0,1)*(sigma[1]-sigma[1].T.conj())
##    Tepoint = np.trace(np.dot(subbroad1,np.dot(subG, np.dot(subbroad2,subG.T.conj()))))
#    t = transmission(G,en,gamma1,gamma2)
#    dosg.append(-np.imag(np.trace(G))/np.pi)
#    T.append(t)
#    
#pdosopen = np.diag(spect)
#ax4.legend(loc='upper right')
#ax4.plot(np.linspace(-8,8,1000),T, label = "B = 20T", linestyle ='--')
#
#pylab.title('Conductance vs. E for different B')
##
#plt.savefig("transmission-armchair-ribbonvsjunction.png", format='png', bbox_inches='tight', transparent = True)
''' END : COMPARISON OF RIBOON VS JUNCTION '''



''' START : TRANSMISSION VS ALPHA '''
#cell = SimulationCell(Lattice = "custom", M = ribbon(6,2,"a"))    
#tbyalpha = []
#for alp in np.linspace(0,6,100):
#    H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = alp, orthogonal = True, Bcenter=np.array([10.67,1.19,0]))
#    T = []
#    dosg = []
#    N = len(H[0])
#    spect = np.zeros((N,N),dtype=complex)
#    for en in np.linspace(-0.2,0.2,100):
#        G, sigma, subsigma = greensopen(leaddict["ribbon_6_2_a"], H, en)
##    subG = gsubmat("ribbon_4_1_z",0,1,G)
#        gamma1, subgamma1 = broadening(sigma[0],subsigma[0])
#        gamma2, subgamma2 = broadening(sigma[1],subsigma[1])
#        spect += spectral(G)
##    gamma1 = np.complex(0,1)*(sigma[0]-sigma[0].T.conj())
##    gamma2 = np.complex(0,1)*(sigma[1]-sigma[1].T.conj())
##    Tepoint = np.trace(np.dot(subbroad1,np.dot(subG, np.dot(subbroad2,subG.T.conj()))))
#        t = transmission(G,en,gamma1,gamma2)
#        dosg.append(-np.imag(np.trace(G))/np.pi)
#        T.append(t)
#    tbyalpha.append(sum(T))
#
#plt.plot(np.linspace(0,6,100),tbyalpha)
#
#def transmission(celltype,lead1,lead2,H,en):
#    G, sigma, subsigma = greensopen(leaddict[celltype], H, en)
#    subG = gsubmat(celltype,0,1,G)
#    broad1, subbroad1 = broadening(sigma[0],subsigma[0])
#    broad2, subbroad2 = broadening(sigma[1],subsigma[1])
#    Tepoint = np.trace(np.dot(subbroad1,np.dot(subG, np.dot(subbroad2,subG.T.conj()))))
#    Tepoint = np.real(Tepoint)
''' END : TRANSMISSION BY ALPHA '''



#normdosg = sum(dosg)
#dosg = [item/normdosg for item in dosg]
#
#dos = DOS(np.linspace(-8,8,100),eigvals[0],0.01)
#normdos = sum(dos)
#dos = [item/normdos for item in dos]
#plt.subplot(2,2,2)
#plotDOS(np.linspace(-8,8,100), dos)
#plt.subplot(2,2,3)
#plt.plot(np.linspace(-8,8,100),dosg, c="g")

#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 1, orthogonal = True)
#plt.plot(np.linspace(-8,8,100),dosG)

#cell=SimulationCell(1,1,Lattice="PG")
#kpath1 = kpath(100)
#H,S,eigvals,eigvecs = solveSecularAO(cell,kpath1)
#for i in range(len(eigvals[0])):
#    plt.plot(list(range(len(kpath1))),eigvals[:,i],c="b")
    
#atoms = flake(15,"z")
#end = time.time()
#print("Cell creation",end-start)
#
#start=time.time()
#cell.populateWlatticepts(atoms)
#plotcell(cell)
#end = time.time()
#print("Lattice population",end-start)
#
#plotdosandpdos(cell,bars = 300, gamma=gamma, delta=0.05,center=0)
#E = np.linspace(0,0.2,500)
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0, orthogonal = True)
#
#leadz = [[1, 'armchair', 1, [8, 10, 6, 11, 7, 9]],[2, 'armchair', 1, [3, 1, 5, 0, 4, 2]]]
#leadz = [[1, 'z', 1, [7,9,8,6]],[2, 'z', 1, [2,0,1,3]]]
#H, G, selfens = greens(leadz, H, 0.5)  
#magoff = DOS(E,eigvals[0],E[4]-E[0])
#H, S, eigvals, eigvecs = solveSecularSO(cell, gamma, alpha = 0.4, orthogonal = True)
#magon = DOS(E,eigvals[0],E[4]-E[0])
#diff = [magon[i]-magoff[i] for i in range(len(magon))]
#plotDOS(E,diff)

#plotdosandpdos(cell,Emin=0,Emax=0.2,Ebars=300,Esigma=0.01,PDcenter=0.05,PDdelta=0.01,PDsigma=0.01,PDslices=1,center=np.array([0,0,0]))



#tube,index = fluxtubes(cell,sigma=0.5,B=0.5)
#for atom in cell.atoms:
#    plt.scatter(atom.coord[0], atom.coord[1], s = 10, c="m")
####    plt.annotate(str((atom.Spaceindex[0],atom.Spaceindex[1])), (atom.coord[0],atom.coord[1]))
#    plt.annotate(str(atom.MUCindex[1]),(atom.coord[0],atom.coord[1]))
#for i in range(len(index)):
#    plt.annotate(str(np.round(tube[i][1],2)),(cell.atoms[index[i][0]].coord[0]+latConst/5-0.1,cell.atoms[index[i][0]].coord[1]-0.2), size=8)
#
#start = time.time()
#H,S,eigvals,eigvecs = solveSecularSO(cell,gamma,alpha=0.2)

#E = np.linspace(-10,10,1000)  
#dos = DOS(E,eigvals[0],E[2]-E[0])
#plt.plot(E,dos)
#plt.show()
#plotdosandpdos(cell,gamma,)
#plotdosandpdos(cell,gamma,delta=0.05,sigma=0.01,center=0,slices=100)
#center = 0
#delta = 0.05
#slices = 5
#sigma = 2*delta/slices
#E = np.linspace(center-delta,center+delta,slices)
#pdos1 = sumPDOS(E, eigvals[0], eigvecs, S, sigma)
#plotPDOS(cell, pdos1, 50)

#end = time.time()
#print(end-start)

#start = time.time()
#a,b = la.eigh(H)
#end = time.time()
#print(end-start)














