#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 10:27:28 2018

@author: gerbotmd
"""

import numpy as np
import collections 
from hbtr_pressure import satpress, emax

# Geometrics: For the moment I am going to use "Aproximate" measures from myself,
# and some of the values from stolwijk and tanabe 
# as [core, head
#     upper arm, lower arm, hand, upper arm, lower arm, hand
#     upper leg, lower leg, foor, upper leg, lower leg, foot] 
orad = [0.100, 0.210, 
        0.045, 0.035, 0.035, 0.045, 0.035, 0.035,
        0.090, 0.060, 0.050, 0.090, 0.060, 0.050] # outer radius in meters
leng = [0.60, 0.10,
        0.32, 0.27, 0.10, 0.32, 0.27, 0.10,
        0.46, 0.50, 0.28, 0.46, 0.50, 0.28] # length in meters 
# skin fraction
#fskin = [0.0421, 0.04824] # fraction of skin using thickness from salloum

orad = np.asarray(orad)
leng = np.asarray(leng)
#fskin = np.asarray(fskin)

# skin thickness from salloum

tskin = np.array([0.01912, 0.0085, 
                  0.00451, 0.00451, 0.0074, 0.00451, 0.00451, 0.0074,
                  0.01064, 0.01064, 0.0117, 0.01064, 0.01064, 0.0117])

# compute surface area, volume, and inner radius 
sa = 2 * np.pi * orad * leng
irad = orad - tskin
sain = 2 * np.pi * irad * leng # Inner surface area

# base radiative and covective values
hrad = [5.2290, 6.3910,
        4.9966, 4.9966, 3.4860, 4.9966, 4.9966, 3.4860,
        4.6480, 4.6480, 4.6480, 4.6480, 4.6480, 4.6480] # From stolwijk (conv) W/(m^2 K)
hcon = [2.4983, 3.1955,
        3.4860, 3.4860, 3.8927, 3.4860, 3.4860, 3.8927,
        3.1955, 3.1955, 3.4860, 3.1955, 3.1955, 3.4860] # From stolwijk (conv) W/(m^2 K)

hrad = np.asarray(hrad)
hcon = np.asarray(hcon)

# base blood flow
mskn = [5.230, 2.240,
        0.860, 0.450, 0.910, 0.860, 0.450, 0.910,
        0.380, 0.110, 0.450, 0.380, 0.110, 0.450] # from tanabe l/h
mcor = [204.820, 46.210,
        1.760, 0.915, 0.211, 1.760, 0.915, 0.211,
        1.369, 0.160, 0.078, 1.369, 0.160, 0.078] # from tanabe l/h

mskn = np.asarray(mskn)
mcor = np.asarray(mcor)

# base metabolism
qskn = [0.591, 0.131,
        0.050, 0.026, 0.050, 0.050, 0.026, 0.050,
        0.122, 0.023, 0.100, 0.122, 0.023, 0.100] # from tanabe, but perhaps use values from anthropmetric measures
qcor = [58.945, 17.169,
        1.214, 0.345, 0.296, 1.214, 0.345, 0.296,
        1.318, 0.357, 0.213, 1.318, 0.357, 0.213] # from tanabe

qskn = np.asarray(qskn)
qcor = np.asarray(qcor)

# Heat capacities
corcap = np.array([126039.6, 11592.0,
                   6436.8, 4078.8, 615.6, 6436.8, 4078.8, 615.6,
                   20984.4, 9993.6, 910.8, 20984.4, 9993.6, 910.8]) # J/K converted from tanabe
skncap = np.array([5076.0, 1015.2,
                   543.6, 356.4, 356.4, 543.6, 356.4, 356.4,
                   1522.8, 734.4, 450.0, 1522.8, 734.4, 450.0]) # J/K converted from tanabe  
plcap  = 9396 # J/K, from tanabe 

# now compute the effective thermal conductivities from the above 
# effective thermal conductivities of the cores based on the thermal
# condudctivity and mass fraction of each tissue type in the segment
keff = [0.43, 0.66,
        0.4927, 0.4927, 0.4961, 0.4927, 0.4927, 0.4961,
        0.5035, 0.5035, 0.5195, 0.5035, 0.5035, 0.5195] # W/mK using mass fraction from stoljwick and conductivis from UCTI-Filia
keff = np.asarray(keff)

# core 2 skin resistance 
c2sr = 4 * np.pi * leng * (keff / (1 + orad / irad))

coret = np.array([36.0, 36.0,
                  36.0, 36.0, 36.0, 36.0, 36.0, 36.0,
                  36.0, 36.0, 36.0, 36.0, 36.0, 36.0]) 
skint = np.array([36.0, 36.0,
                  36.0, 36.0, 36.0, 36.0, 36.0, 36.0,
                  36.0, 36.0, 36.0, 36.0, 36.0, 36.0])

poolt = 36.0 

coreset = np.array([36.0] * 14)
skinset = np.array([36.0] * 14)

# control fractions
fsum = np.array([0.149 + 0.132 + 0.212, 0.070,
                 0.023, 0.012, 0.092, 0.023, 0.012, 0.092,
                 0.050, 0.025, 0.017, 0.050, 0.025, 0.017])
fswt  = np.array([0.481, 0.081, 
                  0.051, 0.026, 0.016, 0.051, 0.026, 0.016,
                  0.073, 0.036, 0.018, 0.073, 0.036, 0.018])
fshiv = np.array([0.850, 0.020,
                  0.004, 0.026, 0.000, 0.004, 0.026, 0.000,
                  0.023, 0.012, 0.000, 0.023, 0.012, 0.000])
fvcon = np.array([0.195, 0.022,
                  0.022, 0.022, 0.152, 0.022, 0.022, 0.152,
                  0.022, 0.022, 0.152, 0.022, 0.022, 0.152])
fvdil = np.array([0.322, 0.320,
                  0.031, 0.016, 0.061, 0.031, 0.016, 0.061,
                  0.092, 0.023, 0.050, 0.092, 0.023, 0.050])


a = 1.000
rhoc = 1.067 # Wh/(lK) blood heat capacity flow constant

# environment
envt = np.array([24.0] * 14)
airt = np.array([24.0] * 14)

rh = 0.30

# heat flow components:
# 1) Mbase
# 2) Mshiv
# 3) K == Conductivity
# 4) a. B core
#    b. B skin
# 5) Radiation
# 6) Conduction
# 7) Respiration for the core

# for the first part, I am going to use a globals
# but later I should put things into perhaps a parameters class to pass around

def errors():
    errc = coret - coreset
    errs = skint - skinset
    wrmc = np.clip(errc, 0, None)
    wrms = np.clip(errs, 0, None)
    cldc = -np.clip(errc, None, 0)
    clds = -np.clip(errs, None, 0)
    return errc, errs, wrmc, wrms, cldc, clds

def commands():
    # warning hard coded parameters 
    wrm = np.sum(wrms * fsum)
    cld = np.sum(clds * fsum)
    local = 2.0 ** (errs / 10)
    swt  = (371.2 * errc[1] + 33.6 * (wrm - cld)) * fswt * local
    shiv = (24.4 * cldc[1] * cldc) * fshiv
    vcon = (-11.5 * errc[1] - 11.5 * (wrm - cld))
    vdil = (117.0 * errc[1] + 7.5 * (wrm - cld))
    bfskn = (mskn + fvdil * vdil) / (1 + fvcon * vcon) * local
    return swt, shiv, bfskn

def bloodheat():
    """Compute the heat flows due to blood
    
    # todo, change mcore and mskin to use a controled value
    """
    bhfs_skin = a * rhoc * bfskn * (skint - poolt)
    bhfs_core = a * rhoc * (mcor + shiv / 1.16)  * (coret - poolt)
    return bhfs_skin, bhfs_core

def kondheat():
    """Conductive heat transfer between layers
    """
    return c2sr * (coret - skint)

def radconvheat():
    """Heat lost to the enviroment by radiation and conduction
    """
    radf = hrad * (skint - envt) * sa
    conf = hcon * (skint - airt) * sa
    return radf, conf

def respiration():
    """returns the repiration in the core
    from tanabe 
    """
    coreresp = (0.0014 * (34 - airt[1]) 
                + 0.017 * (5.867 - rh * satpress(airt[1])) 
                    * (qskn.sum() + qcor.sum() + shiv.sum()))
    resp = np.zeros(14)
    resp[0] = coreresp
    return resp

# logging 
def log_append(log, *args):
    """flatten and append all the args itno the log as a single line
    """
    logline = []
    for logp in args:
        if isinstance(logp, collections.Iterable):
            # we must flatten this into the log line
            for p in logp:
                logline.append(p)
        else:
            logline.append(logp)
    
    log.append(logline)

# logs
heatflowslog = []
templog = []

# the simulation loop 

dt = 1 # seconds

for i in range( 100000 ):
    # cycle through the timesteps computing the heat flows and then
    # updating the temperatures
    # errors
    errc, errs, wrmc, wrms, cldc, clds = errors()
    swt, shiv, bfskn = commands()
    
    bfs, bfc = bloodheat()
    kf = kondheat()
    rf, cf = radconvheat()
    resp = respiration()
    
    # new segment temperature     
    ncoret = coret + (qcor - resp + shiv - kf - bfc) * dt / corcap
    nskint = skint + (qskn + kf - bfs - cf - rf) * dt / skncap
    
    # new pool temp
    npoolt = poolt + (np.sum(bfs) + np.sum(bfc)) * dt / plcap 
    
    # Set temperatures     
    coret, skint, poolt = ncoret, nskint, npoolt
    
    # set logs
    log_append(heatflowslog, bfs, bfc, kf, rf, cf)
    log_append(templog, coret, skint, poolt)
    
# save out the matrix

np.savez("2layertest.npz", temp=np.asarray(templog), hfs=np.asarray(heatflowslog))
