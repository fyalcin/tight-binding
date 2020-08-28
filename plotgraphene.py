#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 16:14:15 2019

@author: firaty
"""

import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib.patches import Rectangle
from matplotlib.patches import Circle
import numpy as np
latConst = 1.42

from matplotlib import rc

rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

def rotate(point, angle):
    """ This is just a rotation function that rotates the given point
    counterclockwise by <N*angle> degrees """
    rotMatrix = np.array([[np.cos(angle), -np.sin(angle), 0], \
        [np.sin(angle), np.cos(angle), 0],[0, 0, 1]])
    point = np.dot(rotMatrix, point)
    return point

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

def plotgraphene():
    fig1 = plt.figure(1,figsize=(12,8))
    plt.subplot(1,2,1)
    ax1 = plt.axes()
    plt.xlim((-5,10))
    plt.ylim((-5,5))
    centers = np.array([[0,0,0],[3*latConst/2,latConst*np.sqrt(3)/2,0],[3*latConst/2,-latConst*np.sqrt(3)/2,0],\
                        [3*latConst,0,0]])
    hexs = np.empty((0,3))
    for i in centers:
        hexs = np.vstack([hexs,hexagon(center = i)])
    hexs = np.round(hexs,2)
    hexs = np.unique(hexs,axis=0)
    centers -= np.array([latConst,0,0])
    hexs -= np.array([latConst,0,0])
    for item in centers:
#        plt.annotate("["+str(item[1][0])+","+str(item[1][1])+"]", (item[0][0]-0.1,item[0][1]-0.1), size = 6)
        Hex = RegularPolygon((item[0],item[1]), numVertices=6, radius=latConst, orientation=np.radians(30), \
                             facecolor=(0,0,0,0), edgecolor = (0,0,0,1), linewidth = 2)
        ax1.add_patch(Hex)
#    plt.scatter(hexs[:,0],hexs[:,1], c="black", s=6)
    plt.arrow(0,0,-latConst/2,latConst*np.sqrt(3)/2)
    plt.annotate("", xy=(-latConst/3,latConst*np.sqrt(3)/3), xytext =(0,0), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(-latConst/3,-latConst*np.sqrt(3)/3), xytext =(0,0), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(latConst/(1.5),0), xytext =(0,0), arrowprops = dict(arrowstyle="->"),size=15)
    
    plt.annotate("", xy=(0.1+5*latConst/2,0.05+latConst*np.sqrt(3)/2), xytext =(1*latConst,0), arrowprops = dict(arrowstyle="->"), size = 15)
    plt.annotate("", xy=(0.1+5*latConst/2,-0.05-latConst*np.sqrt(3)/2), xytext =(1*latConst,0), arrowprops = dict(arrowstyle="->"), size=15)
    plt.annotate(r'$\vec{a_1}$',(2*latConst,-0.05+latConst/3),size=12)
    plt.annotate(r'$\vec{a_2}$',(2*latConst,-0.03-latConst/2),size=12)
    plt.annotate(r'$\vec{\delta_2}$', (-0.25,0.65), size=12)
    plt.annotate(r'$\vec{\delta_1}$', (0.5,0.15), size=12)
    plt.annotate(r'$\vec{\delta_3}$', (-0.25,-1), size=12)
    plt.annotate(r'$A$',(-0.15,2.6),size=12)
    plt.annotate(r'$B$',(latConst-0.15,2.60),size=12)
#    plt.yticks([], [])
#    plt.xticks([], [])
#    for item in [fig1, ax1]:
#        item.patch.set_visible(False)
    ax1.axis('off')
    fig1.savefig('honeycombehres.png', dpi = 150, transparent=True, bbox_inches= 'tight')
    return hexs

plotgraphene()
import matplotlib.lines as mlines



def plotBZ():
#    plt.close()
    fig1 = plt.figure(1,figsize=(8,8))
    ax1 = plt.axes()
    hex1 = hexagon()
    Hex = RegularPolygon((0,0), numVertices=6, radius=4*np.pi/(3*np.sqrt(3)*latConst), \
                             facecolor=(0,0,0,0), edgecolor = (0,0,0,1), linewidth = 2)
    ax1.add_patch(Hex)
    mag = 2*np.pi/(3*latConst)
    plt.annotate("", xy=(mag,mag*np.sqrt(3)), xytext =(0,-0.03), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(mag,-mag*np.sqrt(3)), xytext =(0,0.03), arrowprops = dict(arrowstyle="->"),size=15)
    plt.xlim((-4,4))
    plt.ylim((-4,4))
    ax1.axis('off')
    l1 = mlines.Line2D([0.02,2*np.pi/(3*latConst)], [0,0])
    l2 = mlines.Line2D([2*np.pi/(3*latConst),2*np.pi/(3*latConst)], [0,2*np.pi/(3*np.sqrt(3)*latConst)])
    l3 = mlines.Line2D([2*np.pi/(3*latConst),0.02], [2*np.pi/(3*np.sqrt(3)*latConst),0])
    ax1.add_line(l1)
    ax1.add_line(l2)
    ax1.add_line(l3)
    plt.annotate(r'$\Gamma$', (-0.15,0), size=13)
    plt.annotate("M", (1.55,0), size=12)
    plt.annotate("K", (1.55,0.80), size=12)
    plt.annotate("K'", (1.55,-0.90), size=12)
    plt.annotate(r'$\vec{b_1}$', (0.9,2), size=12)
    plt.annotate(r'$\vec{b_2}$', (0.9,-2.2), size=12)
    fig1.savefig('grapheneBZ.png', transparent=True, dpi=300, bbox_inches='tight')
    


#plotBZ()
    
def plotflux():
    center = np.array([0,0,0])
    vec = np.array([[1,0,0],[0,1,0]])
    sqlat = np.empty((0,3))
    fig1 = plt.figure(1,figsize=(4,8))
    ax1 = plt.axes()
    for i in range(3):
        for j in range(5):
            newcents = center + i*vec[0] + j*vec[1]
            sqlat = np.vstack([sqlat,newcents])
    sqlat = np.vstack([sqlat,sqlat+np.array([5,0,0])])
    for i in sqlat:
#        Hex = RegularPolygon((i[0],i[1]), numVertices=6, radius=1, \
#                             facecolor=(0,0,0,0), edgecolor = (0,0,0,1), linewidth = 2)
        Hex = Rectangle([i[0],i[1]], 1, 1, angle=0.0,facecolor=(0,0,0,0),edgecolor = (0,0,0,1))
        ax1.add_patch(Hex)
    ax1.add_patch(Circle([1.5,1.5],0.4,facecolor="grey"))
#    plt.scatter(sqlat[:,0],sqlat[:,1])
#    plt.xlim((-1,9))
#    plt.ylim((-1,9))
    plt.xlim((4,9))
    plt.ylim((-1,9))
    plt.annotate(r'$\mathbf{m}$', (0.9,-0.35), size=15)
    plt.annotate(r'$\mathbf{m+1}$', (1.55,-0.35), size=15)
    plt.annotate(r'$\mathbf{n}$', (-0.5,0.92), size=15)
    plt.annotate(r'$\mathbf{n+1}$', (-0.9,1.92), size=15)
    plt.annotate(r'$\Phi_1$', (1.3,1.35), size=15)
#    plt.annotate(r'$\Phi_1$', (1.4,2.1), size=10)
#    plt.annotate(r'$\Phi_1$', (1.4,3.1), size=10)
#    plt.annotate(r'$\Phi_1$', (1.4,4.1), size=10)
#    plt.annotate(r'$\Phi_1$', (1.4,5.1), size=10)
    plt.annotate("", xy=(1.65,1.997), xytext =(1,1.997), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(1.65,2.997), xytext =(1,2.997), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(1.65,3.997), xytext =(1,3.997), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(1.65,4.997), xytext =(1,4.997), arrowprops = dict(arrowstyle="->"),size=15)
    
    ax1.add_patch(Circle([6.5,1.5],0.4,facecolor="grey"))
    ax1.add_patch(Circle([6.5,3.5],0.4,facecolor="grey"))
    plt.annotate(r'$\mathbf{m}$', (5.9,-0.35), size=15)
    plt.annotate(r'$\mathbf{m+1}$', (6.55,-0.35), size=15)
    plt.annotate(r'$\mathbf{n_1}$', (4.5,0.92), size=15)
#    plt.annotate(r'$n_1+1$', (4.3,1.92), size=12)
    plt.annotate(r'$\mathbf{n_2}$', (4.5,2.92), size=15)
    plt.annotate(r'$\Phi_1$', (6.3,1.35), size=15)
    plt.annotate(r'$\Phi_2$', (6.3,3.35), size=15)
#    plt.annotate(r'$\Phi_1$', (6.4,2.1), size=10)
#    plt.annotate(r'$\Phi_1$', (6.4,3.1), size=10)
#    plt.annotate(r'$\Phi_1$', (6.4,4.1), size=10)
#    plt.annotate(r'$\Phi_1$', (6.4,5.1), size=10)
    plt.annotate("", xy=(6.65,1.997), xytext =(6,1.997), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(6.65,2.997), xytext =(6,2.997), arrowprops = dict(arrowstyle="->"),size=15)
    
    plt.annotate("", xy=(6.7,3.997), xytext =(6,3.997), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(6.7,4.997), xytext =(6,4.997), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(6.6,3.997), xytext =(6,3.997), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(6.6,4.997), xytext =(6,4.997), arrowprops = dict(arrowstyle="->"),size=15)
    ax1.axis('off')
    plt.annotate(r'$(a)$', (-0.5,5.25), size=18)
#    plt.annotate(r'$(b)$', (4.5,5.25), size=18)
    
    fig1.savefig('optimalgaugelowres.png', transparent=True, bbox_inches = 'tight')
    return sqlat

#a = plotflux()
            
def plotsquare():
    centers = np.empty((0,3))
    fig = plt.figure(figsize=(16,6))
    ax = plt.axes()
    for i in range(4):
        for j in range(4):
            centers = np.vstack((centers,i*np.array([1,0,0])+j*np.array([0,1,0])))

    
    centers = np.vstack((centers, centers + np.array([5,0,0])))
    plt.scatter(centers[:,0],centers[:,1], s=150)
    ax.axhline(y = 0, xmin=0, xmax = 0.45, linestyle = "--")
    ax.axhline(y = 1, xmin=0, xmax = 0.45, linestyle = "--")
    ax.axhline(y = 2, xmin=0, xmax = 0.45, linestyle = "--")
    ax.axhline(y = 3, xmin=0, xmax = 0.45, linestyle = "--")
    ax.axvline(x = 0, ymin=0, ymax = 3, linestyle = "--")
    ax.axvline(x = 1, ymin=0, ymax = 3, linestyle = "--")
    ax.axvline(x = 2, ymin=0, ymax = 3, linestyle = "--")
    ax.axvline(x = 3, ymin=0, ymax = 3, linestyle = "--")
    
    ax.axhline(y = 0, xmin=0.55, xmax = 1, linestyle = "--")
    ax.axhline(y = 1, xmin=0.55, xmax = 1, linestyle = "--")
    ax.axhline(y = 2, xmin=0.55, xmax = 1, linestyle = "--")
    ax.axhline(y = 3, xmin=0.55, xmax = 1, linestyle = "--")
    ax.axvline(x = 5, ymin=0, ymax = 3, linestyle = "--")
    ax.axvline(x = 6, ymin=0, ymax = 3, linestyle = "--")
    ax.axvline(x = 7, ymin=0, ymax = 3, linestyle = "--")
    ax.axvline(x = 8, ymin=0, ymax = 3, linestyle = "--")
    Hex = Rectangle([1,1], 6, 1, angle=0.0,facecolor=(0,0.4,0.4,0.1),edgecolor = (0,0,0,1))
    ax.add_patch(Hex)
    plt.annotate("", xy=(1,2), xytext =(1,1), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(2,1), xytext =(1,1), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate(r'$\vec{a_2}$', (0.80,1.4), size=15)
    plt.annotate(r'$\vec{a_1}$', (1.45,0.80), size=15)
    plt.annotate(r'$q \times \vec{a_1}$', (3.85,0.56), size=20)
    plt.annotate("", xy=(7,0.5), xytext =(1,0.50), arrowprops = dict(arrowstyle="->"),size=15)
    ax.axis("off")
    plt.savefig("MUC.png", bbox_inches = 'tight', transparent = True)
    return centers
#a = plotsquare()

def brick():
    plt.figure(figsize=(6,6))
    ax = plt.axes()
    centers = np.array([[0,0,0],[0,2*np.sqrt(3)*latConst/2],[3*latConst/2,np.sqrt(3)*latConst/2,0],[3*latConst/2,3*np.sqrt(3)*latConst/2,0]])
    for item in centers:
        Hex = RegularPolygon((item[0],item[1]), numVertices=6, radius=latConst, orientation=np.radians(30), \
                             facecolor=(0,0,0,0), edgecolor = (0,0,0,1), linewidth = 2)
        ax.add_patch(Hex)
    plt.annotate("", xy=(5.7,2), xytext=(4.5,2), arrowprops = dict(arrowstyle="->"),size=15)
    centers = np.array([[7,-1.5,0], [9,0,0],[7,1.5,0],[9,3,0]])
    for item in centers:
        Hex = Rectangle([item[0],item[1]], 2, 3, angle=0.0,facecolor=(0,0.4,0.4,0),edgecolor = (0,0,0,1),linewidth = 2)
        ax.add_patch(Hex)
    plt.annotate(r'$m$', (8.9,-2.1), size=10)
    plt.annotate(r'$m-1$', (6.4,-2.1), size=10)
    plt.annotate(r'$m+1$', (10.4,-0.7), size=10)
    plt.annotate(r'$n$', (6.5,1.4), size=10)
    plt.annotate(r'$n-1$', (7.8,-0.1), size=10)
    plt.annotate(r'$n+1$', (7.8,2.8), size=10)
    plt.xlim((-2.5,12.5))
    ax.axis("off")
    plt.ylim((-5.5,9.5))
    plt.savefig("brick.png", transparent = True, bbox_inches = 'tight')
#brick()


def plotsquarenofield():
    centers = np.empty((0,3))
    fig = plt.figure(figsize=(8,8))
    ax = plt.axes()
    for i in range(4):
        for j in range(4):
            centers = np.vstack((centers,i*np.array([1,0,0])+j*np.array([0,1,0])))

    
#    centers = np.vstack((centers, centers + np.array([5,0,0])))
    plt.scatter(centers[:,0],centers[:,1], s=150)
    ax.axhline(y = 0, xmin=0, xmax = 1, linestyle = "--")
    ax.axhline(y = 1, xmin=0, xmax = 1, linestyle = "--")
    ax.axhline(y = 2, xmin=0, xmax = 1, linestyle = "--")
    ax.axhline(y = 3, xmin=0, xmax = 1, linestyle = "--")
    ax.axvline(x = 0, ymin=0, ymax = 3, linestyle = "--")
    ax.axvline(x = 1, ymin=0, ymax = 3, linestyle = "--")
    ax.axvline(x = 2, ymin=0, ymax = 3, linestyle = "--")
    ax.axvline(x = 3, ymin=0, ymax = 3, linestyle = "--")
#    ax.axvline(x = 4, ymin=0, ymax = 3, linestyle = "--")
#    ax.axhline(y = 0, xmin=0.55, xmax = 1, linestyle = "--")
#    ax.axhline(y = 1, xmin=0.55, xmax = 1, linestyle = "--")
#    ax.axhline(y = 2, xmin=0.55, xmax = 1, linestyle = "--")
#    ax.axhline(y = 3, xmin=0.55, xmax = 1, linestyle = "--")
#    ax.axvline(x = 5, ymin=0, ymax = 3, linestyle = "--")
#    ax.axvline(x = 6, ymin=0, ymax = 3, linestyle = "--")
#    ax.axvline(x = 7, ymin=0, ymax = 3, linestyle = "--")
#    ax.axvline(x = 8, ymin=0, ymax = 3, linestyle = "--")
#    Hex = Rectangle([1,1], 6, 1, angle=0.0,facecolor=(0,0.4,0.4,0.1),edgecolor = (0,0,0,1))
#    ax.add_patch(Hex)
    plt.annotate("", xy=(1,2), xytext =(1,1), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate("", xy=(2,1), xytext =(1,1), arrowprops = dict(arrowstyle="->"),size=15)
    plt.annotate(r'$\mathbf{a_2}$', (0.70,1.4), size=25)
    plt.annotate(r'$\mathbf{a_1}$', (1.40,0.75), size=25)
#    plt.annotate(r'$q \times \vec{a_1}$', (3.85,0.56), size=20)
#    plt.annotate("", xy=(7,0.5), xytext =(1,0.50), arrowprops = dict(arrowstyle="->"),size=15)
    ax.axis("off")
    plt.savefig("sqplot.png", bbox_inches = 'tight', transparent = True)
    return centers
#plotsquarenofield()


def plotleads():
    plt.figure()
    ax = plt.axes()
    ax.set_xlim((-2.5,2.5))
    ax.set_ylim((-4,4))
    plt.axis('off')
    Hex1 = Rectangle([-1,-1.5], 2, 3, angle=0.0,facecolor='cadetblue',edgecolor = (0,0,0,0.3))
    Hex2 = Rectangle([-10,-0.5], 9, 1, angle=0.0,facecolor=(0,0.4,0.4,0.1),edgecolor = (0,0,0,0.3))
    Hex3 = Rectangle([1,-0.5], 9, 1, angle=0.0,facecolor=(0,0.4,0.4,0.1),edgecolor = (0,0,0,0.3))
    plt.annotate('Lead 1', (-2.3,-0.3), size=25)
    plt.annotate('Lead 2', (1.25,-0.3), size=25)
    plt.annotate('Device', (-0.5,-0.2), size=25)
    ax.add_patch(Hex1)
    ax.add_patch(Hex2)
    ax.add_patch(Hex3)
    plt.savefig('devicelead.png',dpi=200, bbox_inches='tight', transparent = True)
#plotleads()
#plotsquarenofield()
    
#M1 = flake(5,"z")
#a = open("mag-z.txt")
#b = a.read()
#b = [item for item in b.split("\n")]
#del b[-1]
#b = [100*float(item) for item in b]
#plt.figure()
#plt.scatter(M1[:,0],M1[:,1], s = b, c = b, cmap = "coolwarm")
#plt.savefig("zigzag-flake-magnetization.png", dpi=200, transparent = True, bbox_inches='tight')
#
#
#M1 = flake(3,"a")
#a = open("mag-a.txt")
#b = a.read()
#b = [item for item in b.split("\n")]
#del b[-1]
#b = [400*float(item) for item in b]
#plt.figure()
#plt.scatter(M1[:,0],M1[:,1], s = b, c = b, cmap = "coolwarm")
#plt.savefig("armchair-flake-magnetization.png", dpi=200, transparent = True, bbox_inches='tight')

def nanostruct():
    fig = plt.figure(figsize=(13.63,6))
    cell = SimulationCell(1, 1, Lattice="custom", M = flake(7,"z"))
    cell.rotate(np.pi/6)
    for i in cell.atoms:
        plt.scatter(i.coord[0],i.coord[1], c="black", s=50)
        
#    cell = SimulationCell(1, 1, Lattice="custom", M = triangle(7,"z"))
#    cell.rotate(np.pi/6)
#    cell.translate(np.array([40,0,0]))
#    
#    for i in cell.atoms:
#        plt.scatter(i.coord[0],i.coord[1], c="grey", s=10)
    
    
    cell = SimulationCell(1, 1, Lattice="custom", M = junction(7,2,"z"))
    cell.rotate(np.pi/2)
    cell.translate(np.array([25,0,0]))
    
    for i in cell.atoms:
        plt.scatter(i.coord[0],i.coord[1], c="black", s=50)
        
    cell = SimulationCell(1, 1, Lattice="custom", M = crosssheet(6,2))
    cell.rotate(np.pi/2)
    cell.translate(np.array([57,0,0]))
    
    for i in cell.atoms:
        plt.scatter(i.coord[0],i.coord[1], c="black", s=50)
    
    ax = plt.axes()
    ax.axis("off")
    
    plt.savefig("nanostructures.png", bbox_inches='tight',transparent=True)

#nanostruct()









