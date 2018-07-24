#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 11:50:03 2018

@author: gerbotmd

Functions for partial pressure in the air
"""
import numpy as np

def satpress(temp):
    """Return the saturation pressure of water in air at the temperature
    temp. This function could be *anything* so a table lookup, a 
    function, etc. Here its going to be the Antione Equation for 
    water from 1 to 100 C
    
    Temp in C, returns pressure in mmHg
    """
    A, B, C = 8.07131, 1730.63, 233.426 # constants from ddbonline
    return 10 ** (A - B / (C + temp))

def emax(tskin, pair, lr, hc, area):
    """Returns the maximum sweating rate for each segment
    """
    pskin = satpress(tskin)
    return (pskin - pair) * lr * hc * area

