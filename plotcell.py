#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
import numpy as np
from constants import *
import datetime

__all__ = ["plotcell"]

def plotcell(cell, window = "neighbors", neighbors = False):
    xcoord = [atom.coord[0] for atom in cell.atoms]
    ycoord = [atom.coord[1] for atom in cell.atoms]
    margin = latConst
    xmax = max(xcoord) + margin
    ymax = max(ycoord) + margin
    xmin = min(xcoord) - margin
    ymin = min(ycoord) - margin
    now = datetime.datetime.now()
    datestr = str(now.hour) + "-" + str(now.minute) + "-" + str(now.second)
#    if cell.Lattice is "custom":
#        plt.figure(1,figsize=(8,6))
#        ax1 = plt.axes()
#        plt.scatter(xcoord, ycoord, s=30, c="r")
#        for item in cell.gridxy:
#            plt.annotate("["+str(item[1][0])+","+str(item[1][1])+"]", (item[0][0]-0.1,item[0][1]-0.1), size = 6)
#            Hex = RegularPolygon((item[0][0],item[0][1]), numVertices=6, radius=latConst,orientation=np.radians(30), facecolor=(1,0,0,0.1*item[2]), edgecolor = (0,0,0,0.1))
#            ax1.add_patch(Hex)    
#        plt.savefig('cell_1_' + datestr + '.pdf', format='pdf')
        
#    plt.figure(2,figsize=(8,6))
    
    if window is "all":
        plt.scatter(xcoord, ycoord, s=30, c="b")
        plt.xlim((xmin,xmax))
        plt.ylim((ymin,ymax))
        for atom in cell.atoms:
            plt.annotate("["+str(atom.Spaceindexalt[0])+","+str(atom.Spaceindexalt[1])+"]", (atom.coord[0]-0.3, atom.coord[1]+0.2), size=6)
#    plt.savefig('cell_2_' + datestr + '.pdf', format='pdf')
        
    plt.figure(3,figsize=(8,6))
    plt.scatter(xcoord, ycoord, s=30, c="grey")
    plt.xlim((xmin,xmax))
    plt.ylim((ymin,ymax))
    for atom in cell.atoms:
        plt.annotate(str(len(atom.neighbors)), (atom.coord[0]-0.4, atom.coord[1]+0.2), size=6)
#    plt.savefig('cell_3_' + datestr + '.pdf', format='pdf')
    
    if window is "all":
        plt.figure(4,figsize=(8,6))
        plt.scatter(xcoord, ycoord, s=30, c="g")
        plt.xlim((xmin,xmax))
        plt.ylim((ymin,ymax))
        for atom in cell.atoms:
            plt.annotate(str(atom.MUCindex[1]), (atom.coord[0]-0.4, atom.coord[1]+0.2), size=6)
#    plt.savefig('cell_5_' + datestr + '.pdf', format='pdf')
        
    if neighbors:
        plt.figure(5,figsize=(8,6))
        plt.scatter(xcoord, ycoord, s=30, c="r")
        nxcoord = [atom.coord[0] for atom in cell.neighborCell]
        nycoord = [atom.coord[1] for atom in cell.neighborCell]
        plt.scatter(nxcoord, nycoord, s=2, c="b")
        for atom in cell.neighborCell:
            plt.annotate(str(atom.Spaceindex), (atom.coord[0]-0.2, atom.coord[1]+0.2), size=6)
#        plt.savefig('cell_4_' + datestr + '.pdf', format='pdf')
