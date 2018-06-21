#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:31:11 2018

@author: gerbotmd
"""

import numpy as np
#import pdb

class Segment:
    def __init__(self, store, env, qmet, name="Segment", itemp=37.0):
        """builds a segment, that has children and parents
        """
        self.parent = None # should be a segment object or none
        self.childlinks = None # will end up being a list (segment, bloodcoff) tuples
        self.store = store
        self.env = env 
        self.qmet = qmet
        self.name = name
        self.temp = itemp
        
    def add_child(self, blcoeff, child):
        """Add a child to the segment tree
        """
        child.parent = self
        if self.childlinks is None:
            self.childlinks = [(child, blcoeff)]
        else:
            self.childlinks.append((child, blcoeff))
        return self # needed to allow construction of more compicated trees
    
    #################################################################
    # Coefficent Functions                                          #
    #################################################################
    def coeff_mat(self, dt):
        """Return the coefficent matrix with this segement as parent
        """
        size = self._get_size()
        cmat = np.zeros((size,size)) # here numpy should work
        cid = -1
        _, cmat = self._coeff_mat(cid, cmat, dt, 0.0)
        return cmat
        
    def _coeff_mat(self, cid, cmat, dt, pbcf):
        """Internal recursing function
        """
        cid += 1
        # self term (non, child depened)
        cmat[cid][cid] = self.store / dt + pbcf + self.env
        if self.childlinks is not None:
            for child, bcf in self.childlinks:
                sid = cid
                cmat[sid][sid]  += bcf
                cmat[sid][cid+1] = -bcf
                cmat[cid+1][sid] = -bcf
                cid, cmat = child._coeff_mat(cid, cmat, dt, bcf)
        return cid, cmat
    
    #################################################################
    # CONSTANTS VECTOR                                              #
    #################################################################
    def get_const(self, dt, tinf):
        """Get the constants terms for the system
        """
        cvec = np.zeros(self._get_size())
        cid = -1 #start before
        cvec, _ = self._get_const(cvec, cid, dt, tinf)
        return cvec
    
    def _get_const(self, cvec, cid, dt, tinf):
        """compute the constants
        """
        cid += 1
        cvec[cid] = self.qmet + (self.store * self.temp / dt) + (self.env * tinf)
        if self.childlinks is not None:
            for child, _ in self.childlinks:
                cvec, cid = child._get_const(cvec, cid, dt, tinf)
        return cvec, cid
    
    #################################################################
    # TEMPERATURE FUNCTIONS                                         #
    #################################################################
    def get_temps(self):
        """get temperatures
        """
        temps = np.zeros(self._get_size())
        cid = -1
        temps, _ = self._get_temps(temps, cid)
        return temps
    
    def _get_temps(self, temps, cid):
        """internal
        """
        cid += 1
        temps[cid] = self.temp
        if self.childlinks is not None:
            for child, _ in self.childlinks:
                temps, cid = child._get_temps(temps, cid)
        return temps, cid
    
    def set_temps(self, temps):
        """Set the temperatures, from a list
        """
        temps.reverse()
        self._set_temps(temps)
        
    def _set_temps(self, tstack):
        """Internal
        """
        self.temp = tstack.pop()
        if self.childlinks is not None:
            for child, _ in self.childlinks:
                tstack = child._set_temps(tstack)
        return tstack
    
    #################################################################
    # NAMES FUNCTIONS                                               #
    #################################################################
    def get_names(self):
        """get the names of the segments
        """
        names = []
        return self._get_names(names)
    
    def _get_names(self, names):
        """internal
        """
        names.append(self.name) 
        if self.childlinks is not None:
            for child, _ in self.childlinks:
                names = child._get_names(names)
        return names
            
    def set_names(self, names):
        """given a list of the proper length, set the names
        """
        # list size guard
        if len(names) != self._get_size():
            print("WRONG SIZE NAME LIST!")
            return 
        
        names.reverse()
        self._set_names(names)
        
    def _set_names(self, nstack):
        self.name = nstack.pop()
        if self.childlinks is not None:
            for child, _ in self.childlinks:
                nstack = child._set_names(nstack)
        return nstack
    
    def build_name_matrix(self):
        """Build a Matrix whoes elemnts are strings describing
        the interactions and testing the matrix construction
        """
        size = self._get_size()
        nmat = [["" for i in range(size)] for j in range(size)] #np.full((size,size), "", dtype=np.str)
        cid = -1 # start with -1 since each call sets itself
        _, nmat = self._build_name_matrix(cid, nmat)
        return nmat
    
    def _build_name_matrix(self, cid, nmat):
        """Internal function to build the name matrix as a
        template for building the matrix of the system
        """
        #pid = cid
        #pdb.set_trace()
        cid += 1
        # add self term to matrix, and return
        nmat[cid][cid] = "-".join((self.name, "self"))
        if self.childlinks is not None:
            # here we have to preform the recursions
            sid = cid # self ID
            for child, bcf in self.childlinks:
                nmat[sid][cid+1] = "-".join((self.name, child.name))
                nmat[cid+1][sid] = "-".join((self.name, child.name))
                cid, nmat = child._build_name_matrix(cid, nmat)
        return cid, nmat
                
        
    def _get_size(self):
        """Recursively walk the tree in depth first computing
        the size
        """
        if self.childlinks is None:
            return 1
        else:
            size = 1
            for child, bcf in self.childlinks:
                size += child._get_size()
            return size
        