#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 11:31:11 2018

@author: gerbotmd
"""

import numpy as np
import copy
#from collections import deque
#import pdb

class Segment:
    def __init__(self, store, env, qmet, name="Segment", itemp=37.0):
        """builds a segment, that has children and parents
        """
        self.parent = None # should be a segment object or none
        self.childlinks = None # will end up being a list (segment, bloodcoff) tuples
        # Base parameters
        self.name = name
        self.temp = itemp
        self.store = store
        self.env = env
        self.qmet = qmet
        self._idx = None # This index of the node
        self.matf = _matf # default matrix function
        # set the initail values 
        #self._eval_temperature()
        
    def add_child(self, blcoeff, child, reindex=True):
        """Add a child to the segment tree
        Reindex by default. Turn off to save walking the entier tree with each
        added child. This may be usefull if adding many children to the tree
        or when building a subtree that will be added to a main tree later
        """
        child.parent = self
        if self.childlinks is None:
            self.childlinks = [(child, blcoeff)]
        else:
            self.childlinks.append((child, blcoeff))
        # Reindex the tree as needed in perorder from the root
        if reindex:
            root = self._get_root()
            root.reindex()
        
        return self # needed to allow construction of more compicated trees
    
    #################################################################
    #                        Element Access                         #
    #################################################################
    
    def __getitem__(self, key):
        """ Allow attribute access by key (so that it is easy to build
        access to multiple different parmeters)
        """
        return getattr(self, key)
    
    def __setitem__(self, key, val):
        """Set the item by value by key
        """
        #TODO add checking (e.g. for spelling) so that only certaint
        # parmeters are setable
        setattr(self, key, val)
    
    #################################################################
    #                    Tree Walking functions                     #
    #################################################################
    
    def _walk_touch(self, toucher):
        """Walks the tree in Pre-order, touching each node with the toucher
        function, which generaly will call some method on each node
        """
        nodeq = [self,]
        while len(nodeq) != 0:
            cnode = nodeq.pop()
            toucher(cnode)
            if cnode.childlinks is not None:
                for child, _ in reversed(cnode.childlinks):
                    nodeq.append(child)
        
    def _walk_get(self, retriever):
        """Walks the tree in Pre-order, building a list to return of
        the values at each tree by running retriever, which takes a 
        node, and a collecter, appending the value to the collector
        """
        collector = []
        nodeq = [self,]
        while len(nodeq) != 0:   
            #print([x.name for x in nodeq])
            cnode = nodeq.pop()
            retriever(cnode, collector)
            if cnode.childlinks is not None:
                for child, _ in reversed(cnode.childlinks):
                    nodeq.append(child)
        return collector
    
    def _walk_set(self, setter, valstack):
        """Walks the tree in Per-order, with a setter poping from a reveresed
        valstack to set the tree. setter should be a function that takes a
        node and a valstack. In the setter, some value should be poped from
        valstack and then, the nodes member value should be set to that
        
        This method prints an error and
        returns early if valstack is the wrong size
        """
        if len(valstack) != self._get_size():
            print("ERROR: Wrong Sized Set!") # perhaps should throw an error
            return 
        # setup the stack of values in the right order
        valstack.reverse()
        nodeq = [self,]
        while len(nodeq) != 0:
            cnode = nodeq.pop()
            setter(cnode, valstack)
            if cnode.childlinks is not None:
                for child, _ in reversed(cnode.childlinks):
                    nodeq.append(child)
        
    def _get_root(self):
        """Returns the root of the tree, usefull for reseting the node indexes
        """
        if self.parent is None:
            return self
        else:
            return self.parent._get_root()
        
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
        
    def reindex(self):
        """Walk the tree setting the indexes on the nodes to be used for matrix
        generation. Indexes are set in preorder
        """
        idxer = Idxer(0)
        self._walk_touch(idxer)


    #################################################################
    #                       Coefficent Functions                    #
    #################################################################
    def coeff_mat(self, dt):
        """Return the coefficent matrix with this segement as parent
        
        DEPRICATED
        """
        size = self._get_size()
        cmat = np.zeros((size,size)) # here numpy should work
        cid = -1
        _, cmat = self._coeff_mat(cid, cmat, dt, 0.0)
        return cmat
        
    def _coeff_mat(self, cid, cmat, dt, pbcf):
        """Internal recursing function
        
        DEPRICATED
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
    
    def _walk_matrix(self, matf, matrix):
        """Walk the tree, and for each node, running a mat function to build
        a matrix. matf, should take a touple of the node and the parent
        link value, and the matrix that is going to be built (passed in in the
        correct form)
        """
        # unlike other walkers, each node is given as a touple of node
        # and parent link value (the blood flow coefficent)
        nodeqlnk = [(self, 0.0),]
        while len(nodeqlnk) != 0:
            nodel = nodeqlnk.pop()
            matf(nodel, matrix)
            if nodel[0].childlinks is not None:
                for clink in reversed(nodel[0].childlinks):
                    nodeqlnk.append(clink)
    
    def build_temp_matrix(self, dt):
        """New Function for building the matrix correctly using the walk
        Matrix function
        
        DEPRICATED
        """
        # build the per node function
        def matf(nlnk, matrix):
            node, pbcf = nlnk
            sid = node._idx 
            # self term
            matrix[sid][sid] = node.store / dt + node.env + pbcf
            if node.childlinks is not None:
                for child, lnk in node.childlinks:
                    cid = child._idx
                    # self term addition for children
                    matrix[sid][sid] += lnk
                    # cross terms
                    matrix[sid][cid] = -lnk
                    matrix[cid][sid] = -lnk
        # and then call the function correctly
        size = self._get_size()
        mat = np.zeros((size,size))
        self._walk_matrix(matf, mat)
        return mat 
    
    def build_func_matrix(self, dt):
        """In this case, each of the nodes is going to use its own 
        function, defined for/stored on the node as the .matf parameter
        which will point to a funciton. Otherwise this will call a defualt
        function that has the same properties as the build temp matrix 
        function above, matf must also take dt
        """
        # pass trhough for the matrix function
        def matf(nlnk, mat):
            node, _ = nlnk
            node.matf(nlnk, mat, dt)
        # call the function
        size = self._get_size()
        mat = np.zeros((size,size))
        self._walk_matrix(matf, mat)
        return mat
            
    #################################################################
    #  GET and SET each of the Parameters                           #
    #################################################################
    def get_param(self, pname, nptype=False):
        """Get the parameters in perorder as a list (perhaps add flag
        to return as Numpy vector)
        """
        params = self._walk_get(lambda n, c: c.append(n[pname]))
        if nptype:
            return np.asarray(params)
        else:
            return params
        
    def set_param(self, pname, params):
        """Set the parameters in perorder from the preorder stack
        """
        def pset(node, stack):
            node[pname] = stack.pop()
            
        params = copy.copy(params)
        self._walk_set(pset, params)
    
    #################################################################
    # CONSTANTS VECTOR                                              #
    #################################################################
    def get_const(self, dt, tinf):
        """Get the constants terms for the system
        """
        cvec = self._walk_get(
                lambda n, c: c.append(n.qmet + (n.store * n.temp / dt) + (n.env * tinf) 
                ))
        return np.asarray(cvec)
    
    #################################################################
    # TEMPERATURE FUNCTIONS                                         #
    #################################################################
    def get_temps(self):
        """get temperatures, returns numpy array
        """
        temps = self._walk_get(lambda n, c: c.append(n.temp))
        return np.asarray(temps)
    
    def set_temps(self, temps):
        """Set the temperatures, from a list
        """
        # setting func
        def tset(node, stack):
            node.temp = stack.pop()
        
        temps = copy.copy(temps)
        self._walk_set(tset, temps)
        #self._walk_touch(lambda n: n._eval_temperature())
        
    #################################################################
    # NAMES FUNCTIONS                                               #
    #################################################################
    def get_names(self):
        """get the names of the segments
        """
        return self._walk_get(lambda n, c: c.append(n.name))
    
    def set_names(self, names):
        """given a list of the proper length, set the names
        """
        def nset(node, stack):
            node.name = stack.pop()
        
        names = copy.copy(names)
        self._walk_set(nset, names)
        
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
                
    
###############################################################################
#                              BODY CLASS                                     #
###############################################################################
        
class Body:
    def __init__(self, body_tree=None):
        """Body is an object that will hold the "Whole Body" parameters for 
        the HBTR model (shuch as a chilled vs warm paramter). Update the
        tree (temperature and temperature depenant parameters). Build 
        a tree from the paramters 
        """
        self.body_tree = body_tree
        self.body_param_funcs = {}
        self.body_param_values = {}
        self.update_funcs = {}
        self.tree_log = []
        self.body_log = []
        
    #################################################################
    #                  BODY PARAMETERS                              #
    #################################################################
    
    def register_body_parameter(self, pname, pfunc, pstart):
        """Register a parameter to the body with name pname, and initial
        value pstart. pfunc should be a function that takes the body_tree, and
        old paramter as its argument and then returns the updated paramter. The
        function may not have ato return anything
        """
        self.body_param_funcs[pname] = pfunc
        self.body_param_values[pname] = pstart
        
    def register_body_constant(self, pname, pconst):
        """register a constant body parameter (i.e. a body paramter that
        remains constant for all time
        """
        self.body_param_funcs[pname] = lambda bt, pold: pconst
        self.body_param_values[pname] = pconst
        
    def update_body_parameters(self):
        """Update the body paramters, this should be called whenever the 
        paramters of the tree are updated (e.g. in the update phase of the 
        timestep). Cannonically this will be called after temperature updates,
        but before tree paramter updates (as the update functions for those
        might rely on control parmeters that are body paramters)
        """
        for pname in self.body_param_funcs:
            self.body_param_values[pname] = self.body_param_funcs[pname](
                    self.body_tree, self.body_param_values[pname])
    
    #################################################################
    #                      PARAMETER UPDATES                        #
    #################################################################
    
    def register_tree_update(self, tpname, tpfunc):
        """Resiter a function to update a tree parameter. This function should
        take the body_tree and a dictionary of body parameters and then
        calculate and return the updated body tree values as a preorder list
        """
        self.update_funcs[tpname] = tpfunc
        
    def update_tree_parameters(self):
        """Updatae the tree parameters by runing the functions registered with
        register_tree_update
        """
        for tparam in self.update_funcs:
            nparams = self.update_funcs[tparam](self.body_tree, 
                                       self.body_param_values)
            self.body_tree.set_param(tparam, nparams)
            
    #################################################################
    #                    DEF parameter loging                       #
    #################################################################
    
    def register_log_parameter(self, pname, ptype):
        """ Register a paramter to be loged, either a tree paramter
        or a body parameter. Calling log parameter 
        """
        if ptype == "tree":
            self.tree_log.append(pname)
        elif ptype == "body":
            self.body_log.append(pname)
            
    def _setup_log(self):
        """ Initialize a blank logging dictionary
        """
        self.log = {}
        for p in self.tree_log:
            self.log[p] = []
        for p in self.body_log:
            self.log[p] = []
        
    def _log_parameters(self):
        """ When called log out the parmeter into a dictionary. 
        """
        for p in self.tree_log:
            self.log[p].append(self.body_tree.get_param(p))
        for p in self.body_log:
            self.log[p].append(self.body_param_values[p])
    #################################################################
    #                   Setup and Run The Model                     #
    #################################################################
    
    def run_constant_temp(self, dt, env_temp, steps):
        """Run the simuation at a constant singular enviromental temperature
        for a number of steps
        """
        # To start, update the body and tree parameters
        self._setup_log()
        #self._log_parameters()
        #self.update_body_parameters()
        #self.update_tree_parameters()
        self._log_parameters()
        # create the storage for the temperatures
        output_temps = []
        # And then run the simulation loop 
        for i in range(steps):
            # compute the new temperatures 
            A = self.body_tree.build_func_matrix(dt) # temperature A in [[A]][t] = [c]
            c = self.body_tree.get_const(dt, env_temp) # constant vector
            ntemps = list(np.linalg.inv(A) @ c)
            # update the outputs, body temperature, body params, and matrix params
            output_temps.append(copy.copy(ntemps))
            self.body_tree.set_temps(ntemps)
            self.update_body_parameters()
            self.update_tree_parameters()
            self._log_parameters()
        #Return the temperatures
        return self.log
    
    
###############################################################################
#               Testing Functions and Utilities                               #
###############################################################################

def _matf(nlnk, mat, dt):
    """Default matrix function
    """
    node, pbcf = nlnk
    sid = node._idx 
    # self term
    mat[sid][sid] = node.store / dt + node.env + pbcf
    if node.childlinks is not None:
        for child, lnk in node.childlinks:
            cid = child._idx
            # self term addition for children
            mat[sid][sid] += lnk
            # cross terms
            mat[sid][cid] = -lnk
            mat[cid][sid] = -lnk
    
class Idxer:
    def __init__(self, sidx):
        """Utilitiy class to index nodes
        """
        self.cid = sidx
    def __call__(self, node):
        """Call for each node
        """
        node._idx = self.cid
        self.cid += 1


def _build_testbody():
    """return a test tree for checking the functionality of various refactors
    so that I can get the code to run nicely
    """
    tree = Segment(1, 2, 3,"A")
    tree.add_child(0.1, Segment(4, 5, 6, "B"))
    tree.add_child(0.2, 
                   Segment(7, 8, 9, "C").add_child(0.3, 
                   Segment(10, 11, 12, "D")))
    tree.add_child(0.5, Segment(13, 15, 16, "E"))
    return Body(tree)

