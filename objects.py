#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from lattices import *
from tools import initializer, NeighborStr, dot, check, flake
import numpy as np
import copy
from constants import latConst
import cProfile as profile
import time

__all__ = ["SimulationCell", "Carbon"]

""" Main object we build our lattice with """
class Carbon(object):
    @initializer
    def __init__(self, coord, Rvec=np.zeros(3), UCindex=[[0,0],0], SCindex=[[0,0],0], MUCindex=[[0,0],0], neighbors=[]):
        pass
    def dump(self):
        print(repr(self.__dict__))
    def translate(self, translationVec):
        self.coord += translationVec
    def __str__(self):
        return 'Carbon ' + str(self.MUCindex[1]) + ' at ' + str(np.round(self.coord,2))+ " with UCindex " +\
        str(self.UCindex) + " and SCindex " +str(self.SCindex) + '\n' + 'with ' \
        + 'neighbors ' + '\n' + NeighborStr(self)

'''  Our main simulation cell composed of carbon objects '''
class SimulationCell:
    @initializer
    def __init__(self, N1=0, N2=0, n1=0, n2=0, Lattice="square", q=1, M=flake(5,"z"), periodic = True):
        """ Initialization parameters determine the size of the 
        simulation cell that will be used """
        if Lattice is "custom":
            self.periodic = False
            self.populateWlatticepts(M)
        else:
            self.basisVec, self.UClatVec = getlatinfo(Lattice)
            self.populateCell(N1, N2, n1, n2)
            if periodic:
                self.duplicate(q,1)
            else:
                for index, atom in enumerate(self.atoms):
                    atom.MUCindex = [[0,0], index]
            self.populateNeighbors(1.1*latConst, 1, 1)
        
    def bar(self):
        profile.runctx('self.populateNeighbors(1.5,1,1)', globals(), locals())
    def baralt(self):
        profile.runctx('self.populateNeighborsAlt(1.5,1,1)', globals(), locals())
    
    def dump(self):
        print(repr(self.__dict__))
    
    def populateCell(self, N1, N2, n1, n2):
        cellRvec = n1*self.UClatVec[0] + n2*self.UClatVec[1]
        self.SClatVec = np.array([(N1-n1)*self.UClatVec[0], (N2-n2)*self.UClatVec[1]])
        self.atoms = [Carbon(self.basisVec[k]+i*self.UClatVec[0]+j*self.UClatVec[1], Rvec = cellRvec, \
        UCindex = [np.array([i,j]),k]) for i in range(N1) for j in range(N2)\
        for k in range(len(self.basisVec))]
        for index, atom in enumerate(self.atoms):
            atom.SCindex = [np.array([0,0]), index]
            
    def populateHexagons(self, M):
        angle = np.pi/3
        basis = []
        for i in range(6):
            rotMatrix = np.array([[np.cos(i*angle), -np.sin(i*angle), 0], \
                        [np.sin(i*angle), np.cos(i*angle), 0],[0, 0, 1]])
            rotatedcoord = np.dot(rotMatrix, np.array([latConst,0,0]))
            basis.append(rotatedcoord)
        rows = len(M)
        cols = len(M[0])
        tmpCell = []
        self.gridxy = []
        for row in range(rows):
            for col in range(cols):
                transVec = np.array([col*(3*latConst/2), (rows-row+(col%2)/2)*(latConst*np.sqrt(3)), 0])
                if M[row][col] == 1:
                    newcoords = basis + transVec
                    for coord in np.round(newcoords,4):
                        tmpCell.append(coord)
                self.gridxy.append([transVec,[row,col],M[row][col]])
        hexagons = np.unique(tmpCell, axis=0)
        self.atoms = [Carbon(coord, MUCindex = [[0,0],index]) for index, coord in enumerate(hexagons)]
        
    def populateWlatticepts(self, array):
        self.periodic = False
        self.atoms = []
        self.atoms = [Carbon(coord, MUCindex = [[0,0],index]) for index, coord in enumerate(array)]
        self.populateNeighborsAlt(1.1*latConst, 1, 1)
        
    def duplicate(self, M1, M2):
        newCell = []
        index = 0
        for m1 in range(M1):
            for m2 in range(M2):
                for atom in self.atoms:
                    tmpAtom = copy.deepcopy(atom)
                    tmpAtom.translate(m1*self.SClatVec[0]+m2*self.SClatVec[1])
                    tmpAtom.UCindex[0] += np.array([m1*self.N1, m2*self.N2])
                    tmpAtom.MUCindex = [np.array([0,0]),index]
                    tmpAtom.SCindex[0] = [m1,m2]
                    newCell.append(tmpAtom)
                    index += 1
        self.MUClatVec = np.array([M1*self.SClatVec[0], M2*self.SClatVec[1]])
        self.atoms = newCell
        
#    def fluxtubes(self, magnitude, center, sigma):
#        tubes = []
#        for atom in self.atoms:
#            for neighbor in atom.neighbors:
#                if abs(atom.coord[1]-neighbor.coord[1]) < 0.05:
#                    fluxindex = [(atom.)]
#                    tubes.append()
                
    def translate(self, translationVec):
        """ Simple translation method to translate all the objects in the cell """
        for atom in self.atoms:
            atom.translate(translationVec)
            atom.Rvec += translationVec
            
    def remove(self, removearray):
        removearray = np.sort(removearray)
        for index in removearray[::-1]:
            del self.atoms[index]
        for index, atom in enumerate(self.atoms):
            atom.MUCindex[1] = index
        self.populateNeighbors(1.1*latConst, 1, 1)
            
    def makeFlake(self, N1, N2):
        xlowercut = (N1/2-1)*self.UClatVec[0][0] + latConst/2
        xuppercut = (3*N1/2)*self.UClatVec[0][0] - latConst
        removeList = []
        for atom in self.atoms:
            if atom.coord[0] < xlowercut or atom.coord[0] > xuppercut:           
                removeList.append(atom.SCindex[1])
        for item in removeList[::-1]:
            del self.atoms[item]
        for index, atom in enumerate(self.atoms):
            atom.SCindex[1] = index
            atom.MUCindex = atom.SCindex
            
    def rotate(self, angle, pivot = np.array([0,0,0])):
        """ Rotates the atoms in the cell around a given point in space
        by an axis parallel to z axis """
        rotMatrix = np.array([[np.cos(angle), -np.sin(angle), 0], \
            [np.sin(angle), np.cos(angle), 0],[0, 0, 1]])
        self.translate(-pivot)
        for atom in self.atoms:
            atom.coord = np.dot(rotMatrix, atom.coord)
        self.translate(pivot)
    
    def indexspace(self):
        for atom in self.atoms:
            atom.Spaceindex = [int(np.round(atom.coord[0]/(3*latConst/2))), int(np.round(atom.coord[1]/(0.5*latConst*np.sqrt(3))))]

    def indexspacealt(self):
        for atom in self.atoms:
            atom.Spaceindexalt = [int(np.round((atom.coord[0]+3*latConst/4)/(3*latConst/2))), int(np.round(atom.coord[1]/(0.5*latConst*np.sqrt(3))))]

    def resetNeighbors(self):
        for atom in self.atoms:
            atom.neighbors = []

    def populateNeighbors(self, rCut, N1, N2):
        """ Replicates the main cell and for every object in the main cell,
        creates neighbors based on a distance condition , rCut : cutoff radius"""
        s=time.time()
        self.resetNeighbors()
        superCell = []
        self.indexspace()
        if not self.periodic:
            superCell = copy.deepcopy(self.atoms)
        else:
            for i in range(-N1, N1+1):
                for j in range(-N2, N2+1):
                    tmpCell = copy.deepcopy(self)
                    tmpCell.translate(i*self.MUClatVec[0] + j*self.MUClatVec[1])
                    for atom in tmpCell.atoms:
                        atom.UCindex[0] += np.array([i*self.N1*self.q, j*self.N2]) ## our assumption is that the MUC is generated by duplicating in the first lattice vector direction
                        atom.MUCindex[0] = [i,j]
                        superCell.append(atom)
        for atom in self.atoms:
            atom.neighbors = []
            for neighbor in superCell:
                separationVec = neighbor.coord - atom.coord
                separationDist = dot(separationVec,separationVec)**0.5
#                separationDist = np.linalg.norm(separationVec)
                if separationDist < rCut and separationDist > 0.1:
                    tmpAtom = copy.copy(neighbor)
                    tmpAtom.sepDist = separationDist
                    tmpAtom.sepVec = separationVec
                    tmpAtom.sepAngle = np.arctan2(separationVec[1], separationVec[0])
                    atom.neighbors.append(tmpAtom)
        e=time.time()
        print("populateNeighbors took",e-s,"seconds")
        self.neighborCell = superCell
        
    def populateNeighborsAlt(self, rCut, N1, N2):
        """ Replicates the main cell and for every object in the main cell,
        creates neighbors based on a distance condition , rCut : cutoff radius"""
        s=time.time()
        self.resetNeighbors()
        superCell = []
        self.indexspace()
        self.indexspacealt()
        if not self.periodic:
            superCell = copy.deepcopy(self.atoms)
        else:
            for i in range(-N1, N1+1):
                for j in range(-N2, N2+1):
                    tmpCell = copy.deepcopy(self)
                    tmpCell.translate(i*self.MUClatVec[0] + j*self.MUClatVec[1])
                    for atom in tmpCell.atoms:
                        atom.UCindex[0] += np.array([i*self.N1*self.q, j*self.N2]) ## our assumption is that the MUC is generated by duplicating in the first lattice vector direction
                        atom.MUCindex[0] = [i,j]
                        superCell.append(atom)
        for atom in self.atoms:
            atom.neighbors = []
            for neighbor in superCell:
                sepVec = neighbor.coord - atom.coord
                m1,n1 = atom.Spaceindexalt
                m2,n2 = neighbor.Spaceindexalt
                if ((m1 == m2 and (n2 == n1+1 or n2 == n1-1))) or \
                ((sepVec[1] < 0.05 and sepVec[1]>-0.05) and ((sepVec[0]>-1.45 and sepVec[0]<-1.40) \
                  or (sepVec[0]<1.45 and sepVec[0]>1.40))):
                    neighbor.sepVec = sepVec
                    neighbor.sepAngle = np.arctan2(sepVec[1], sepVec[0])
                    atom.neighbors.append(neighbor)
        e=time.time()
        print("populateNeighbors took",e-s,"seconds")
        self.neighborCell = superCell
        
    def __str__(self):
        atoms = ''
        for atom in self.atoms:
            atoms += atom.__str__() + '\n'*(self.atoms.index(atom) != len(self.atoms)-1)
        return atoms

#cell = SimulationCell(1,1, Lattice = "PG")
#cell.populateWlatticepts(flake(25,"z"))
#cell.bar()