#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 13 13:43:07 2019

@author: firaty
"""

from sympy import init_session
#from sympy.physics.vector import ReferenceFrame
#from sympy.physics.vector import curl

init_session()

A,B,sigma = symbols('A B sigma')
x_prime,y_prime,z_prime=symbols('x\',y\',z\'')
rho,phi,r_prime = symbols('rho,phi,r\'')
sigma = 1
MagFieldCart = Matrix([0,0,B*exp(-(x**2+y**2)/(2*sigma**2))])
MagFieldCyl = Matrix([0,0,B*exp((-rho**2)/(2*sigma**2))])

def curlcart(field):
    ax = diff(field[2],y)-diff(field[1],z)
    ay = diff(field[0],z)-diff(field[2],x)
    az = diff(field[1],x)-diff(field[0],y)
    return Matrix([ax,ay,az])
def curlcyl(field):
    arho = diff(field[2],phi)/rho - diff(field[1],z)
    aphi = diff(field[0],z) - diff(field[2],rho)
    az = (diff(rho*field[1],rho)-diff(field[0],phi))/rho
    return Matrix([arho,aphi,az])
    
VecPotCart = curlcart(MagFieldCart)
VecPotCyl = curlcyl(MagFieldCyl)
#a1 = integrate(VecPot[0],(x,-oo,oo))

#IntegrandCartx = VecPotCart[0]/sqrt((x-x_prime)**2+(y-y_prime)**2 + z**2)
IntegrandCyl = VecPotCyl[1]*rho/sqrt((rho-r_prime)**2+z**2)
import time
#ax = integrate(IntegrandCartx,(x,-oo,oo))

start=time.time()
aphi = integrate(IntegrandCyl,(phi,0,2*pi))
aphirho = integrate(aphi,(rho,0,oo))
#aphirhoz = integrate(aphirho,(z,-oo,oo))
end=time.time()
print(end-start)