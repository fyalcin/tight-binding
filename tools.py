#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import inspect
from functools import wraps
import numpy as np
from constants import latConst, Vpp, Vss, Vsp, Vpi, Spp, Sss, Ssp, Spi, leaddict
import copy
import matplotlib.pyplot as plt
import time
#from objects import SimulationCell

__all__ = ["NeighborStr", "initializer", "gaussianize", "intMatrix", "kpath", \
           "rotate", "kpointsBZ", "dot", "GetPeierlsPhase","BuildJunction", "PDOS", \
           "plotPDOS", "PDOS", "sumPDOS", "DOS", "hexagon", "ribbon", "triangle", "flake",\
           "crosssheet","junction","fluxtubes","plotDOS","PDOSorthogonal","sumPDOSorthogonal",\
           "h_zigzag","h_armchair","self_energy","greensopen","gsubmat","plotflux",\
           "broadening","greensclosed","self_en_iter","spectral","transmission"]

def NeighborStr(atom):
    neighborstr = ''
    for neighbor in atom.neighbors:
        separationVec = neighbor.coord - atom.coord
#        separationSepDist = (separationVec.dot(separationVec))**0.5
        neighborstr += 'Carbon ' + str(neighbor.MUCindex[1]) + ' at coord. ' + \
            str(np.round(neighbor.coord,2)) + ', Rvec ' + str(np.round(neighbor.Rvec,2)) + \
           " and SepDistance "+ str(1.42) + '\n'
    return neighborstr
    
def initializer(fun):
    names, varargs, keywords, defaults = inspect.getargspec(fun)
    @wraps(fun)
    def wrapper(self, *args, **kargs):
        for name, arg in list(zip(names[1:], args)) + list(kargs.items()):
            setattr(self, name, arg)
        for i in range(len(defaults)):
            index = -(i + 1)
            if not hasattr(self, names[index]):
                setattr(self, names[index], defaults[index])
        fun(self, *args, **kargs)
    return wrapper

def check(v1,v2):
    if len(v1)!=len(v2):
        raise ValueError("the length of both arrays must be the same")
    pass

def dot(v1,v2):
    """                                                                                                     
    d0 is Nominal approach:                                                                                 
    multiply/add in a loop                                                                                  
    """
    check(v1,v2)
    out = 0
    for k in range(len(v1)):
        out += v1[k] * v2[k]
    return out

def gaussianize(values, xarray, sigma):
    yarray = []
    total = 0
    for x in xarray:
        total = 0
        for value in values:
            total += np.exp(-(x-value)**2/(2*sigma**2)) 
        yarray.append(total)
    return xarray, yarray

def dist(a,b,sigma):
    return np.exp(-(a-b)**2/(2*sigma**2))

def DOS(E, eigvals, sigma):
    s=time.time()
    dos = []
    for i in range(len(E)):
        total = 0
        for j in range(len(eigvals)):
            total += dist(E[i],eigvals[j],sigma)
        dos.append(total)
    e=time.time()
    print("DOS took",e-s,"seconds")
    return dos
        
def PDOS(E, eigvals, eigvecs, S, sigma):
    eigvecs = eigvecs.T
    pdos = eigvecs[0].conj()*np.dot(S,eigvecs[0])*dist(E,eigvals[0],sigma)
    for i in range(1,len(eigvals)):
        pdos += eigvecs[i].conj()*np.dot(S,eigvecs[i])*dist(E,eigvals[i],sigma)
    return pdos.real

def PDOSorthogonal(E, eigvals, eigvecs, S, sigma):
    eigvecs = eigvecs.T
    pdos = eigvecs[0].conj()*eigvecs[0]*dist(E,eigvals[0],sigma)
    for i in range(1,len(eigvals)):
        pdos += eigvecs[i].conj()*eigvecs[i]*dist(E,eigvals[i],sigma)
    return pdos.real

def sumPDOS(E, eigvals, eigvecs, S, sigma):
    s=time.time()
    totpdos = PDOS(E[0], eigvals, eigvecs, S, sigma)
    for i in range(1, len(E)):
        totpdos += PDOS(E[i], eigvals, eigvecs, S, sigma)
    e=time.time()
    print("sumPDOS took",e-s,"seconds")
    return totpdos

def sumPDOSorthogonal(E, eigvals, eigvecs, S, sigma):
    s=time.time()
    totpdos = PDOSorthogonal(E[0], eigvals, eigvecs, S, sigma)
    for i in range(1, len(E)):
        totpdos += PDOSorthogonal(E[i], eigvals, eigvecs, S, sigma)
    e=time.time()
    print("sumPDOS took",e-s,"seconds")
    return totpdos

def plotDOS(E,dos):
    plt.plot(E,dos, c="r", linewidth = 1)
    plt.axhline(y=0, linestyle = "dashed", c ="black")


def plotPDOS(cell, pdos, size):
    xcoord = [atom.coord[0] for atom in cell.atoms]
    ycoord = [atom.coord[1] for atom in cell.atoms]
    minpdos = 0
    pdos = [(item-minpdos)*size for item in pdos]
    plt.scatter(xcoord, ycoord, c = pdos, cmap= "coolwarm", s = pdos)
#    plt.scatter(xcoord, ycoord, c = pdos, cmap = "RdBu", s = pdos)
#    for index,item in enumerate(pdos):
#        if item < 0:
#            plt.scatter(xcoord[index], ycoord[index], s= size*pdos, c = "b")
#        else:
#            plt.scatter(xcoord[index], ycoord[index], s= size*pdos, c = "r")
    
def GetPeierlsPhase(lattice, atom, neighbor, alpha): 
#    if lattice.Lattice is "custom":
#        if atom.Spaceindex[0] == 0 and neighbor.Spaceindex[0] == 0:
#            PeierlsFactor = -2*np.pi*alpha*int(atom.Spaceindex[1]>5)*np.sign(neighbor.sepVec[0])
#            PeierlsPhase = np.exp(np.complex(0, PeierlsFactor))
#        else:
#            PeierlsPhase = 1
#    if lattice.Lattice is "flake":
#        if atom.Spaceindex[0]+neighbor.Spaceindex[0] == 0 and atom.Spaceindex[1] > 0:
#            PeierlsFactor = 2*np.pi*alpha*np.sign(neighbor.sepVec[0])
#            PeierlsPhase = np.exp(np.complex(0, PeierlsFactor))
#        else:
#            PeierlsPhase = 1
    if lattice.Lattice is "AG" or lattice.Lattice is "PAG":
        if np.abs(neighbor.sepAngle-2*np.pi/3) < 0.1 or np.abs(neighbor.sepAngle+np.pi/3) < 0.1:
            PeierlsFactor = 2*np.pi*alpha*(atom.Spaceindex[0])*np.sign(neighbor.sepAngle)
            PeierlsPhase = np.exp(np.complex(0, PeierlsFactor))
        else:
            PeierlsPhase = 1
    elif lattice.Lattice is "square":
        if np.abs(neighbor.sepAngle-np.pi/2) < 0.1 or np.abs(neighbor.sepAngle+np.pi/2) < 0.1:
            PeierlsFactor = 2*np.pi*alpha*(atom.UCindex[0][0])*np.sign(neighbor.sepAngle)
            PeierlsPhase = np.exp(np.complex(0, PeierlsFactor))
        else:
            PeierlsPhase = 1
    elif lattice.Lattice is "PG":
        if np.abs(neighbor.sepAngle-2*np.pi/3) < 0.1 or np.abs(neighbor.sepAngle+np.pi/3) < 0.1:
            PeierlsFactor = 2*np.pi*alpha*(atom.Spaceindex[0])*np.sign(neighbor.sepAngle)
            PeierlsPhase = np.exp(np.complex(0, PeierlsFactor))
        else:
            PeierlsPhase = 1
    else:
        PeierlsPhase = 1
    return PeierlsPhase

def intMatrix(SepAngle, SepDist, orbitals = "s"):
    """ Creates 4x4 interaction matrices between possible orbitals,
    hopping parameters are scaled to SepDistance by an inverse SepDistance squared 
    relation """
    if orbitals is "s":
        return Vpi/SepDist**2, Spi/SepDist**2
    n11 = Vss
    n12 = -Vsp*np.cos(SepAngle)
    n13 = -Vsp*np.sin(SepAngle)
    n14 = n24 = n34 = n41 = n42 = n43 = 0
    n21 = -n12
    n22 = (-Vpp*np.cos(SepAngle)**2) + (Vpi*np.sin(SepAngle)**2)
    n23 = -(Vpp+Vpi)*np.sin(SepAngle)*np.cos(SepAngle)
    n31 = -n13
    n32 = n23
    n33 = (-Vpp*np.sin(SepAngle)**2) + (Vpi*np.cos(SepAngle)**2)
    n44 = Vpi
      
    m1 = np.array([[n11, n12, n13, n14],
                   [n21, n22, n23, n24],
                   [n31, n32, n33, n34],
                   [n41, n42, n43, n44]])
    m1 = m1/SepDist**2;
               
    s11 = Sss
    s12 = -Ssp*np.cos(SepAngle);
    s13 = -Ssp*np.sin(SepAngle);
    s14 = s24 = s34 = s41 = s42 = s43 = 0
    s21 = -s12
    s22 = (-Spp*np.cos(SepAngle)**2 + Spi*np.sin(SepAngle)**2)
    s23 = -(Spp+Spi)*np.sin(SepAngle)*np.cos(SepAngle)
    s31 = -s13
    s32 = s23
    s33 = (-Spp*np.sin(SepAngle)**2 + Spi*np.cos(SepAngle)**2)
    s44 = Spi
            
    m2 = np.array([[s11, s12, s13, s14],
                   [s21, s22, s23, s24],
                   [s31, s32, s33, s34],
                   [s41, s42, s43, s44]])
    m2 = m2/SepDist**2
    return m1, m2

def kpath(N):
    """ Creates kpoints between the high symmetry points in the
    first BZ of graphene """
#    kpath = []
    kpath = np.empty((0,3), float)
    for n in range(N+1):
        ky = (1-n/N)*2*np.pi/(latConst*3*(3**0.5))
        kx = ky*(3**0.5)
        kpoint = np.array([kx ,ky ,0])
        kpath = np.vstack((kpath, kpoint))
    for n in range(N+1):
        kx = (n/N)*2*np.pi/(latConst*3)
        ky = 0
        kpoint = np.array([kx, ky, 0])
        kpath = np.vstack((kpath, kpoint))
    for n in range(N+1):
        kx = 2*np.pi/(3*latConst)
        ky = (n/N)*2*np.pi/(latConst*3*(3**0.5));
        kpoint = np.array([kx, ky, 0])
        kpath = np.vstack((kpath, kpoint))
    return kpath  

def rotate(point, angle):
    """ This is just a rotation function that rotates the given point
    counterclockwise by <N*angle> degrees """
    rotMatrix = np.array([[np.cos(angle), -np.sin(angle), 0], \
        [np.sin(angle), np.cos(angle), 0],[0, 0, 1]])
    point = np.dot(rotMatrix, point)
    return point

def kpointsBZ(N):
    ''' Generate kpoints inside the irreducible triangle in the first BZ of graphene '''
    prop = (2*np.pi/(3*latConst))
    irred = [np.array([prop*nx/(N-1), prop*ny/(N-1), 0]) for nx in range(N) \
             for ny in range(N) if nx > ny*np.sqrt(3)]
    tmp = copy.deepcopy(irred)
    for point in tmp:
        point[1] = -point[1]
        irred.append(point)
    tmp = copy.copy(irred)
    for n in range(1,6):
        for point in tmp:
            point = rotate(point, n*np.pi/3)
            irred.append(point)
    return irred

def BuildJunction(width, length):
    angle = np.pi/3
    rotMatrix = np.array([[np.cos(angle), -np.sin(angle), 0], \
            [np.sin(angle), np.cos(angle), 0],[0, 0, 1]])
    basis = []
    for i in range(6):
        rotMatrix = np.array([[np.cos(i*angle), -np.sin(i*angle), 0], \
                        [np.sin(i*angle), np.cos(i*angle), 0],[0, 0, 1]])
        rotatedcoord = np.dot(rotMatrix, np.array([latConst,0,0]))
        basis.append(rotatedcoord)
    arms = np.empty((0,3))
    for w in range(width):
        for l in range(length):
            newcoords = basis + np.array([l*(3*latConst/2), -(2*w+l)*(latConst*np.sqrt(3))/2, 0])
            for coord in newcoords: 
                arms = np.vstack((arms, coord))
    mirrorCell = -arms
    arms[:,1] *= -1
    arms = np.vstack((arms, mirrorCell))
    arms = np.unique(arms, axis=0)
    
    body = np.empty((0,3))
    bodyl = length
    for w in range(width):
        for bl in range(bodyl):
            newcoords = basis + np.array([w*(3*latConst/2), -(2*bl+w%2)*(latConst*np.sqrt(3))/2, 0])
            for coord in newcoords:
                body = np.vstack((body,coord))
    body += np.array([-(width//2)*(3*latConst/2), ((width//2)%2)*latConst*np.sqrt(3)/2, 0])
    junction = np.empty((0,3))
    junction = np.vstack((arms, body))
    junction = np.round(junction, 5)
    junction = np.unique(junction, axis=0)
    return junction

def hexagon(orientation="a", center = np.array([0,0,0])):
    angle = np.pi/3
    rotMatrix = np.array([[np.cos(angle), -np.sin(angle), 0], \
            [np.sin(angle), np.cos(angle), 0],[0, 0, 1]])
    basis = []
    for i in range(6):
        rotMatrix = np.array([[np.cos(i*angle), -np.sin(i*angle), 0], \
                        [np.sin(i*angle), np.cos(i*angle), 0],[0, 0, 1]])
        rotatedcoord = np.dot(rotMatrix, np.array([latConst,0,0]))
        basis.append(rotatedcoord)
    if orientation == "z":
        for i,j in enumerate(basis):
            basis[i] = rotate(basis[i], np.pi/2)
    basis += center
    return np.asarray(basis)

def ribbon(length,width,orientation = "a", start = np.array([0,0,0]),build = False):
    basis = hexagon(orientation,start)
    if orientation == "a":
        transVec = np.array([3*latConst,np.sqrt(3)*latConst,0])
    else:
        transVec = np.array([np.sqrt(3)*latConst, 3*latConst, 0])
    base = np.empty((0,3))
    for i in range(length):
        for j in range(width):
            newcoords = np.array([i*transVec[0], j*transVec[1],0]) + basis
            base = np.vstack([base,newcoords])
    if build is False:
        base = np.round(base,4)
    base = np.unique(base, axis=0)
    if build is False and orientation is "z":
        for i,j in enumerate(base):
            base[i] = rotate(base[i],np.pi/2)
    return base
    
def triangle(side, orientation="a"):
    if orientation == "a":
        transVec = np.array([3*latConst/2,3*np.sqrt(3)*latConst/2,0])
    else:
        transVec = np.array([np.sqrt(3)*latConst/2, 3*latConst/2, 0])
    base = ribbon(side,1,orientation, build = True)
    atoms = base
    for i in range(1,side):
        newlayer = ribbon(side-i,1,orientation, build=True) + np.array([i*transVec[0], i*transVec[1], 0])
        atoms = np.vstack([atoms,newlayer])
    xcm = 0
    ycm = 0
    offset = np.floor((side+1)/4)
    for i in atoms:
        xcm += i[0]/len(atoms)
        ycm += i[1]/len(atoms)
    cm = np.array([xcm,0,0])
    atoms += -cm
    if orientation is "z":
        offset = np.floor((side-1)/4) -1
        atoms += np.array([0,-(offset)*3*latConst,0])
    elif orientation is "a":
        offset = (side+1)/2 -1 
        atoms += np.array([0,-(offset)*np.sqrt(3)*latConst,0])
    if orientation == "z":
        for i,j in enumerate(atoms):
            atoms[i] = rotate(atoms[i],np.pi/2)
    atoms = np.round(atoms, 2)
    atoms = np.unique(atoms, axis=0)
    return atoms

#@profile
def flake(side, orientation="a"):
    if orientation == "a":
        transVec = np.array([[3*latConst/2,3*np.sqrt(3)*latConst/2,0],[3*latConst/2,-3*np.sqrt(3)*latConst/2,0]])
    else:
        transVec = np.array([[np.sqrt(3)*latConst/2, 3*latConst/2, 0],[np.sqrt(3)*latConst/2, -3*latConst/2, 0]])
    base = hexagon(orientation)
    atoms = np.empty((0,3))
    for i in range(side):
        for j in range(side):
            tmp = base + i*transVec[0] + j*transVec[1]
            atoms = np.vstack([atoms,tmp])
    removeList = []
    for index,atom in enumerate(atoms):
        cutoff = transVec[0][1]*((side+1*int(orientation == "z"))/2)
        if atom[1] > cutoff or atom[1] < -cutoff:
            removeList.append(index)
    for index in removeList[::-1]:
        atoms = np.delete(atoms,index, axis=0)
    xcm=0
    for atom in atoms:
        xcm += atom[0]
    xcm = xcm/len(atoms)
    cm = np.array([xcm,0,0])
    atoms -= cm
    if orientation == "z":
        for i,j in enumerate(atoms):
            atoms[i] = rotate(atoms[i],np.pi/2)
    atoms = np.round(atoms,2)
    atoms = np.unique(atoms,axis=0)
    return atoms

def crosssheet(length,width):
    base1 = ribbon(length,width,orientation="a", build=True)
    base2 = ribbon(length,width,orientation="z", build=True)
    xcm1 = (length-1)*3*latConst/2
    ycm1 = (width-1)*np.sqrt(3)*latConst/2
    cm1 = np.array([xcm1,ycm1,0])
    base1 -= cm1
    xcm2 = (length-1)*np.sqrt(3)*latConst/2
    ycm2 = (width-1)*3*latConst/2
    cm2 = np.array([xcm2,ycm2,0])
    base2 -= cm2
    for i,j in enumerate(base2):
        base2[i] = rotate(base2[i],np.pi/2)
    cross = np.concatenate((base1,base2),axis=0)
    cross = np.round(cross,2)
    cross = np.unique(cross,axis=0)
    return cross
  
def junction(length,width,orientation="a"):
    base = ribbon(length,width,orientation, build = True)
    if orientation is "a":
        base -= np.array([(width%2-1)*3*latConst/2,(width-1)*np.sqrt(3)*latConst/2,0])
    elif orientation is "z":
        base -= np.array([(width%2-1)*np.sqrt(3)*latConst/2,(width-1)*3*latConst/2,0])
    arm1 = base*np.array([-1,1,1])
    for i,j in enumerate(arm1):
        arm1[i] = rotate(arm1[i], -np.pi/3)
    arm2 = arm1*np.array([1,-1,1])
    body = np.concatenate((base,arm1,arm2),axis=0)
#    body = base
    if orientation is "z":
        for i,j in enumerate(body):
            body[i] = rotate(body[i],-np.pi/2)
    body = np.round(body,2)
    body = np.unique(body,axis=0)
    return body
    
def indexspace(cell):
    indices = [(int(np.round(atom[0]/(3*latConst/2))), int(np.round(atom[1]/(0.5*latConst*np.sqrt(3)))))\
               for atom in cell]
    return indices

def fluxtubes(cell,sigma,B,center=np.array([0,0,0])):
    tubes = []
    indices = []
    for atom in cell.atoms:
        for neighbor in atom.neighbors:
            if abs(atom.coord[1]-neighbor.coord[1]) < 0.05:
                x = (atom.coord[0]+neighbor.coord[0])/2
                y = (atom.coord[1]+neighbor.coord[1])/2 -latConst*np.sqrt(3)/2
                rsq = (x-center[0])**2 + (y-center[1])**2
                tmp = [(atom.Spaceindex,neighbor.Spaceindex),np.round(B*np.exp(-rsq/(2*sigma**2)),2)]
                if tmp not in tubes:
                    tubes.append(tmp)
                    indices.append([atom.MUCindex[1],neighbor.MUCindex[1]])
    tubez = copy.deepcopy(tubes)
    for i in range(len(tubez)):
        for j in range(i):
            if tubez[i][0][0][0] == tubes[j][0][0][0] and tubez[i][0][0][1] > tubes[j][0][0][1]:
                tubez[i][1] += tubes[j][1]
        indices[i].append(tubez[i][1])            
    return tubez, indices

def h_zigzag(width):
    a = Vpi/(latConst**2)
    h00 = a*(np.eye(4*width,k=1)+np.eye(4*width,k=-1))
    baseblock = np.array([[0,a,0,0],[0,0,0,0],[0,0,0,0],[0,0,a,0]])
    h01 = np.kron(np.eye(width),baseblock)
    h01t = h01.T
    return h00,h01,h01t

def h_armchair(width):
    dim = 2*(2*width+1)
    a = Vpi/(latConst**2)
    h00 = a*(np.eye(dim,k=-2)+np.eye(dim,k=2))
    ind1 = [[2*(2*i),2*(2*i)+1] for i in range(width+1)]
    for item in ind1:
        h00[item[0]][item[1]] = a
        h00[item[1]][item[0]] = a
    ind2 = [[2*(2*i-1),2*(2*i)-1] for i in range(1,width+1)]
    h01 = np.zeros((dim,dim))
    h01t = np.zeros((dim,dim))
    for item in ind2:
        h01[item[0]][item[1]] = a
    h01t = h01.T
    return h00,h01,h01t

def self_energy(h00,h01,h01t,E):
    N = len(h00[0])
    maxiter = 1000
    convergence = 1e-6
    epsilon = 1e-3*np.complex(0,1)
    temp = np.linalg.inv((E+epsilon)*np.eye(N)-h00)
    u = np.dot(temp,h01t)
    v = np.dot(temp,h01)
    pv = v
    S = u
    r = 0
    conv = 1
    while conv > convergence:
        r+=1
        temp = np.linalg.inv(np.eye(N)-np.dot(u,v)-np.dot(v,u))
        u = np.dot(temp,np.dot(u,u))
        v = np.dot(temp,np.dot(v,v))
        S = S + np.dot(pv,u)
        pv = np.dot(pv,v)
        conv = abs(np.trace(pv))
        if r > maxiter:
            break
    sigma = np.dot(h01,S)
    return sigma

def self_en_iter(v,vp,vm,E):
    N = len(vp[0])
    conv = 1e-2
    err = 1
    eps = 1e-03*np.complex(0,1)
    g = np.linalg.inv((E+eps)-v)
    itermax = 1000
    while err > conv:
        r += 1
        gnew = np.linalg.inv((E+eps)*np.eye(N)-np.dot(vp,npdot(g,vp.T.conj())))
        err = sum(sum(abs(gnew-g)))
        g = gnew
        if r > itermax:
            break
    return g

def greensopen(leads, Hsave, E):
    subsigma = []
    sigma = []
    H = copy.deepcopy(Hsave)
    for lead in leads:
        if lead[1] is "a":
            if lead[0] == "left":
                a,b,c = h_armchair(lead[2])
            elif lead[0] == "right":
                a,c,b = h_armchair(lead[2])
        elif lead[1] is "z":
            if lead[0] == "left":
                a,b,c = h_zigzag(lead[2])
            elif lead[0] == "right":
                a,c,b = h_zigzag(lead[2])
            
        selfen = self_energy(a,b,c,E)
        subsigma.append(selfen)
        sig = np.zeros((len(H),len(H)),dtype=complex)
        for i,j in enumerate(lead[3]):
            for k,l in enumerate(lead[3]):
                sig[j][l] += selfen[i][k]
                H[j][l] += selfen[i][k]
        sigma.append(sig)
    epsilon = 1e-3*np.complex(0,1)
    G = np.linalg.inv((E)*np.eye(len(H[0])) - H)
    return G, sigma, subsigma

def greensclosed(H,E):
    epsilon = 1e-01*np.complex(0,1)
    G = np.linalg.inv((E-epsilon)*np.eye(len(H[0])) - H)
    return G

def gsubmat(celltype,leadindex1,leadindex2,G):
    leads = leaddict[celltype]
    lead1 = leads[leadindex1]
    lead2 = leads[leadindex2]
    dims = len(lead1[3])
    submat = np.zeros((dims,dims), dtype=complex)
    for i,j in enumerate(lead1[3]):
        for k,l in enumerate(lead2[3]):
            submat[i][k] = G[j][l]
    return submat

def plotflux(cell,sigma,B,Bcenter=np.array([0,0,0])):
    plt.figure(figsize=(6.6,6))
    ax = plt.axes()
    ax.axis("off")
    tube,index = fluxtubes(cell,sigma, B, Bcenter)
    for ind in index:
        x = (cell.atoms[ind[0]].coord[0]+cell.atoms[ind[1]].coord[0])/2
        y = cell.atoms[ind[0]].coord[1]
#        plt.scatter(x, y, s = 5*ind[2], c="m")
        plt.annotate(str(np.round(ind[2],2))+r'$\Phi_0$',(x-0.3,y+0.1))
        plt.annotate("", xy=(cell.atoms[ind[1]].coord[0],cell.atoms[ind[1]].coord[1]), xytext =(cell.atoms[ind[0]].coord[0],cell.atoms[ind[0]].coord[1]), arrowprops = dict(arrowstyle="->"),size=15)
    for atom in cell.atoms:
        plt.scatter(atom.coord[0],atom.coord[1],c="black", s=100)
#    plt.savefig("fluxonflake.png",dpi=300, bbox_inches='tight', transparent =True)
    # update ax.viewLim using the new dataLim


def broadening(sigma, subsigma):
    broad = np.complex(0,1)*(sigma-sigma.T.conj())
    subbroad = np.complex(0,1)*(subsigma-subsigma.T.conj())
    return broad, subbroad

def spectral(G):
    return np.complex(0,1)*(G-G.T.conj())

def transmission(G,E,gamma1,gamma2):
    t = np.trace(np.dot(gamma1,np.dot(G, np.dot(gamma2,G.T.conj()))))
    return np.real(t)

#def rho(mu,H, eigvals, eigvecs):
    

#H, G = greens(leads,H,0.5)
#a,b,c = h_zigzag(1)
#n = self_energy(a,b,c,0.5)

#import time
#
#start = time.time()
#cell = SimulationCell(1,1,"PG")


#cell = SimulationCell(1, 1, Lattice="custom", M = ribbon(5,1,"z"))
#atoms = flake(3,"z")
#cell.populateWlatticepts(atoms)
#tube,index = fluxtubes(cell,sigma=1,B=0.5, center = np.array([0,4.9,0]))
#for atom in cell.atoms:
#    plt.scatter(atom.coord[0], atom.coord[1], s = 10, c="m")
#####    plt.annotate(str((atom.Spaceindex[0],atom.Spaceindex[1])), (atom.coord[0],atom.coord[1]))
##    plt.annotate(str(atom.MUCindex[1]),(atom.coord[0],atom.coord[1]))
#for i in range(len(index)):
#    plt.annotate(str(np.round(tube[i][1],2)),(cell.atoms[index[i][0]].coord[0]+latConst/5-0.1,cell.atoms[index[i][0]].coord[1]-0.2), size=8)
#plotflux(cell,0.5,1)
#H,S,bands,eigvecs = solveSecularSO(cell,gamma)
#ind = indexspace(atoms)
#for i,j in enumerate(atoms):
#    plt.scatter(j[0],j[1], c="r")
#    plt.annotate(str(ind[i]),(j[0],j[1]), size=5)

#tubezz = copy.deepcopy(tube)
#for i in tubezz:
#    for j in tube:
#        if i[0][0][0] == j[0][0][0] and i[0][0][1] > j[0][0][1]:
#            i[1] += j[1]

#for i in range(len(tubezz)):
#    for j in range(i):
#        if tubezz[i][0][0][0] == tube[j][0][0][0] and tubezz[i][0][0][1] > tube[j][0][0][1]:
#            tubezz[i][1] += tube[j][1]
        
#end = time.time()
#print((end-start))

#for index,atom in enumerate(a):
#    plt.scatter(atom[0],atom[1], c="r")
#    plt.annotate(str(index), (atom[0], atom[1]))
    
#fig1 = plt.figure(1, figsize = (9,9))
#
#plt.subplot(2,2,1)
#lattice = triangle(13,"z")
#plt.title(str(len(lattice)) + " atoms")
#for atom in lattice:
#    plt.scatter(atom[0],atom[1], c="r", s=7)
#    
#plt.subplot(2,2,2)
#lattice = triangle(11,"a")
#plt.title(str(len(lattice)) + " atoms")
#for atom in lattice:
#    plt.scatter(atom[0],atom[1], c="r", s=5)
#
#plt.subplot(2,2,3)
#lattice = flake(13,"z")
#plt.title(str(len(lattice)) + " atoms")
#for atom in lattice:
#    plt.scatter(atom[0],atom[1], c="r", s=7)
#
#plt.subplot(2,2,4)
#lattice = flake(11,"a")
#plt.title(str(len(lattice)) + " atoms")
#for atom in lattice:
#    plt.scatter(atom[0],atom[1], c="r", s=5)
#
#plt.savefig("nanostructures.png", dpi=300)
#a = triangle(100,"z")