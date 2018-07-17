#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 28 15:37:30 2018

@author: Matthew D. Gerboth

New Segment model that operates with a stronger separation of concerns
and is agnostic to the parameters to be stored in it. Creating this
body this way will allow a greater diverstiy of segment types
"""

import numpy as np
import copy
from collections import OrderedDict

###############################################################################
#                    DEFAULT MATRIX CONSTRUCTOR                               #
###############################################################################

def _pass_func(*arg,**kwargs):
    """Function for default matrix constructor that just passes
    """
    pass

def _name_mat(node, matrix, *args, **kwargs):
    """Default matrix constructor, builds a matrix of the names of the nodes
    """
    sid = node.idx()
    matrix[sid][sid] = node.name
    for child in node.children:
        cid = child.idx()
        matrix[sid][cid] = "{}-{}".format(node.name, child.name)
        matrix[cid][sid] = "{}-{}".format(child.name, node.name)

def _name_vect(node, vector, *args, **kwargs):
    """Default vector constructor, builds a vector of the names of the nodes
    """
    sid = node.idx()
    vector[sid] = node.name
    
def sm_fcep_1_matrix(node, matrix, dt, *args, **kwargs):
    """Original matrix function. For each node the heat balance is of the form:
    store * (tnew - told) / dt = Qmet + MC_prox(T_n-1-T_n) - MC_dist(T_n-T_n+1)
    -env(T_n-T_env)
    
    This is the matrix function that was used in the previous segment model
    
    Assumes that each node has parameters:
        "store" and "env"
        and each link has a parameter "bfc"
    """
    # ID and values from parent
    sid = node.idx()
    if node.parent is not None:
        pbfc = node.parent.lnkprms[node]["bfc"]
    else:
        pbfc = 0.0
    # self term
    matrix[sid][sid] = node["store"] / dt + node["env"] + pbfc
    for child in node.children:
        cid = child.idx()
        cbfc = node.lnkprms[child]["bfc"]
        # additional self term
        matrix[sid][sid] += cbfc
        # cross terms
        matrix[sid][cid] = -cbfc
        matrix[cid][sid] = -cbfc
        
def sm_fcep_1_vect(node, vect, dt, tenv, *args, **kwargs):
    """Original vector function.
    
    Assumes that each node has parameters
        "qmet"
        "temp" # old / inital temperature
        "store"
    """
    nid = node.idx()
    nconst = (node["qmet"] 
              + (node["store"] * node["temp"] / dt )
              + (node["env"] * tenv))
    vect[nid] = nconst

###############################################################################
#                            BODY NODES                                       #
###############################################################################

class BodyNode:
    def __init__(self, name, params=None,  
                 matf=_pass_func, vectf=_pass_func, 
                 parent=None, plink={},
                 children=[], lnkprms=[], reindex=True):
        """Create a node in the body hierarchy that holds a number of
        parameters that can be used to compute the the matrix to solve the
        problem At = b, where t is the vector of the node temperatures 
        (and perhaps bloodflows), and b is some vector
        
        Each node will have children nodes, and these links may have parameters
        associated with them, such as the value of the blood flow etc.
        children are stored as a list: children
        and parameters are stored in a dictionary of dictionaries: lnkprms 
            that associates the children objects with link parameters 
        
        name::   Name of the body node (for debugging and printing)

        Optional Parameter:
        params::   Dictionary of parameters to be accessed by the node
        matf::     Function that takes a node and a matrix and fills in a number
                   of elements of the matrix
        vectf::    Returns the RHS of the matrix equation
        parent::   Parent Node of this node
        plink::    Dictionary of parmeters to add to parent node describing the
                   link to this node
        children:: List of children nodes for this node
        lnkprams:: List of dictionaries of parameters for each child node 
                   connection, this is converted to a dictionary between
                   child nodes and parent nodes
        """
        # node parameters 
        if params is None:
            self._params = {}
        else:
            self._params = params
        self.matf = matf
        self.vectf = vectf
        self.name = name
        self._idx = -1 # index of the node (for matrix construction)
        # Tree setup
        if parent is not None:
            # make sure parent is aware of its children at creation
            parent.add_child(self, plink, False) # don't reindex, we can't walk this node yet
        else:
            self.parent = None
        if len(lnkprms) != len(children):
            #if we have no link parmeters, then make an empty list of them
            if len(lnkprms) == 0:
                lnkprms = [{} for _ in range(len(children))]
            else:
                print("Opps! lnkprms and children not equal length")
        self.children = []
        self.lnkprms = {} # not a list this time! 
        for child, lnk in zip(children, lnkprms):
            self.add_child(child, lnk, False)
        if reindex:
            self.get_root().reindex()
        
    ## ATTRIBUTE SELECTION
    def __getitem__(self, key):
        """ Get a parameter from the node
        """
        if key == "name":
            return self.name    
        else:
            return self._params[key]
  
    def __setitem__(self, key, val):
        """Set a parameter of the node
        """
        if key == "name":
            self.name = val
        else:
            self._params[key] = val
    
    ## SAFE NODE INDEX
    def idx(self):
        """Returns the nodes index. This may check to make sure the tree is
        indexed correctly, that is erroring on a -1
        """
        if self._idx >= 0:
            return self._idx
        else:
            print("INDEX NOT SET! node: {}".format(self.name))
    
    ## TREE MANIPULATION
    def add_child(self, child, lnk=None, reindx=True):
        """Add a child to a node, and possibly adding link parameters
        """
        child.parent = self
        self.children.append(child)
        if lnk is None:
            lnk = {} # If this value is defualt, then all point to the same dictioanry
        self.lnkprms[child] = lnk
        if reindx:
            self.get_root().reindex()
        return self
        
    def get_root(self):
        """Return the root node of the tree
        """
        if self.parent is None:
            return self
        else:
            return self.parent.get_root()
        
    def reindex(self):
        """From self, reindex the tree
        """
        cid = [0] # fake pass by reference
        def sidx(node, cid):
            node._idx = cid[0]
            cid[0] += 1
        self._walk(sidx, cid)
        
    ## WALKING NODES
    def _walk(self, visitor, *visprms, **viskeys):
        """Walk the tree applying the visitor at each node, the visitor can
        then use the visprms and  viskeys to preform some action
        """
        nodeq = [self,]
        while len(nodeq) != 0:
            cnode = nodeq.pop()
            visitor(cnode, *visprms, **viskeys)
            # add its children to the tree
            for child in reversed(cnode.children):
                nodeq.append(child)
                
    def size(self):
        """Return the size of the tree (for alocating numpy arrays, etc.)
        """
        acc = [0] # fake pass by refrence with a list
        def vf(node, acc):
            acc[0] += 1
        self._walk(vf, acc)
        return acc[0]
        
    def get_param(self, pname):
        """Walk the tree in preorder returning the parameter names as a list
        """
        prms = []
        def pg(node, prms):
            prms.append(node[pname])
        self._walk(pg, prms)
        return prms
    
    def set_param(self, pname, pvals):
        """Walk the tree preorder setting the paramter values from 
        a copy of the list pvals, which should be in preorder
        """
        pvals = copy.copy(pvals)
        pvals.reverse() # treat as a stack poping from the end 
        def ps(node):
            node[pname] = pvals.pop()
        self._walk(ps)
        
    ## WALKING LINKS
    def get_lnk_param(self, pname):
        """Walk the tree in preorder, returning a list of thie link values.
        In general, this should be a list 1 shorter than the parent list
        """
        lvals = []
        def glval(node):
            if node.parent is None:
                return
            # from the node, get the link dict of its parent
            lvals.append(node.parent.lnkprms[node][pname])
        self._walk(glval)
        return lvals

    def set_lnk_param(self, pname, lvals):
        """Walk the tree in preorder, setting the link values. Remember link
        values will be a smaller list than the size of the tree
        """
        lvals = copy.copy(lvals)
        lvals.reverse()
        def slval(node):
            if node.parent is None:
                return
            nval = lvals.pop()
            node.parent.lnkprms[node][pname] = nval
        self._walk(slval)
    
    # MATRIX AND SIMULATION FUCTIONS
    def build_matrix(self, *runprms, **runkeys):
        """Build the matrix for the timestep using the per node matf. This
        this function should accept the runprms and the runkeys (arguments and
        keywords) that can describe information such as timestep, enviromental
        temperature, etc
        
        A data type can be passed in in order to build the matrix
        """
        # note: This may have to change if we are going to have nodes that have
        # more than one equation (say a mass flow balance and a temperature balance)
        size = self.size()
        dtype = runkeys.get('dtype', np.double)
        matrix = np.zeros((size,size), dtype)
            
        def nrun(node, *args, **kwargs):
            node.matf(node, matrix, *args, **kwargs)
        self._walk(nrun, *runprms, **runkeys)
        return matrix 
        
    def build_rhsvect(self, *args, **kwargs):
        """Build the vector for the right hand side of the matrix using the
        vectf fucntion on each node to set its value
        """
        size = self.size()
        dtype = kwargs.get('dtype', np.double)
        vect = np.zeros(size, dtype)
        def nrun(node, *args, **kwargs):
            node.vectf(node, vect, *args, **kwargs)
        self._walk(nrun, *args, **kwargs)
        return vect
  
## BODY CLASS
class Body:
    def __init__(self, bodytree=None):
        """Body is an object that will hold the "Whole Body" parameters for 
        the HBTR model (shuch as a chilled vs warm paramter). Update the
        tree (temperature and temperature depenant parameters). Build 
        a tree from the paramters 
        """
        self.bodytree = bodytree
        self.bodyfuncs = OrderedDict()
        self.bodyvalues = OrderedDict()
        self.treefuncs = OrderedDict()
        self.tree_log = []
        self.body_log = []
        
    #################################################################
    #                  BODY PARAMETERS                              #
    #################################################################
    
    def register_body_parameter(self, pname, pfunc, pstart):
        """Register a parameter to the body with name pname, and initial
        value pstart. pfunc should be a function that takes the bodytree, and
        old paramter as its argument and then returns the updated paramter. The
        function may not have ato return anything
        """
        self.bodyfuncs[pname] = pfunc
        self.bodyvalues[pname] = pstart
        
    def register_body_constant(self, pname, pconst):
        """register a constant body parameter (i.e. a body paramter that
        remains constant for all time
        """
        self.bodyfuncs[pname] = lambda bt, pold: pconst
        self.bodyvalues[pname] = pconst
        
    def update_body_parameters(self):
        """Update the body paramters, this should be called whenever the 
        paramters of the tree are updated (e.g. in the update phase of the 
        timestep). Cannonically this will be called after temperature updates,
        but before tree paramter updates (as the update functions for those
        might rely on control parmeters that are body paramters)
        """
        for pname in self.bodyfuncs:
            self.bodyvalues[pname] = self.bodyfuncs[pname](
                    self.bodytree, self.bodyvalues)
    
    #################################################################
    #                      PARAMETER UPDATES                        #
    #################################################################
    
    def register_tree_update(self, tpname, tpfunc):
        """Resiter a function to update a tree parameter. This function should
        take the bodytree and a dictionary of body parameters and then
        calculate and return the updated body tree values as a preorder list
        """
        self.treefuncs[tpname] = tpfunc
        
    def update_tree_parameters(self):
        """Updatae the tree parameters by runing the functions registered with
        register_tree_update
        """
        for tparam in self.treefuncs:
            nparams = self.treefuncs[tparam](self.bodytree, self.bodyvalues)
            self.bodytree.set_param(tparam, nparams)
            
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
            self.log[p].append(self.bodytree.get_param(p))
        for p in self.body_log:
            self.log[p].append(self.bodyvalues[p])
            
    #################################################################
    #                   Setup and Run The Model                     #
    #################################################################
    
    #TODO: rewrite this function <--> finish rewriting
    def run_constant_temp(self, dt, envlist, steps):
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
            # compute the inital parameters (error, warms, integratinos,
            # etc))
            self.update_tree_parameters()
            self.update_body_parameters()
            A = self.bodytree.build_matrix(dt, *envlist, self) # temperature A in [[A]][t] = [c]
            c = self.bodytree.build_rhsvect(dt, *envlist, self) # constant vector
            # TODO: Add some timestepping voodo? 
            ntemps = list(np.linalg.inv(A) @ c)
            # update the outputs, body temperature, body params, and matrix params
            output_temps.append(copy.copy(ntemps))
            self.bodytree.set_param("temp", ntemps)
            self._log_parameters()
        #Return the temperatures
        return self.log
        
        
###############################################################################
#                             GENERAL FUNCTIONS                               #
###############################################################################



###############################################################################
#                             TESTING FUNCTIONS                               #
###############################################################################

def _test_tree_of_bodynode():
    """Return a test tree's root node
    """
    a = BodyNode("A", {"foo": 1, "bar": 2})
    BodyNode("B", {"foo": 2, "bar": 3}, parent=a)
    a.add_child(BodyNode( "C", {"foo": 4, "bar": 5}, 
         children=[BodyNode("D", {"foo": 6, "bar": 7})]))
    a.add_child(BodyNode("E", {"foo": 8, "bar": 9}))
    return a 
        
def _test_tree_mimic_old():
    """Return a test tree's root node with the same parameterization and 
    functions as the previous segment model
    """
    a = BodyNode("A", {"store": 1.0, "env": 2.0, "qmet": 3.0},
                 sm_fcep_1_matrix, sm_fcep_1_vect)
    a.add_child(BodyNode("B",{"store": 4.0, "env": 5.0, "qmet": 6.0}, 
                         sm_fcep_1_matrix, sm_fcep_1_vect),
                {"bfc": 0.1})
    a.add_child(BodyNode("C", {"store": 7.0, "env": 8.0, "qmet": 9.0}, 
                         sm_fcep_1_matrix, sm_fcep_1_vect,
                         children=[BodyNode("D",
                                            {"store": 10.0, 
                                             "env": 11.0, 
                                             "qmet": 12.0}, 
                                            sm_fcep_1_matrix, 
                                            sm_fcep_1_vect)],
                         lnkprms=[{"bfc": 0.3}]),
                {"bfc": 0.2})
    a.add_child(BodyNode("E", {"store": 13.0, "env": 15.0, "qmet": 16.0}, 
                         sm_fcep_1_matrix, sm_fcep_1_vect),
                {"bfc": 0.5})
    
    a.set_param("temp", [37.0] * 5)
    
    return a
        
def _test_body():
    """Return the test tree above as a body matrix
    """
    tree = _test_tree_mimic_old()
    body = Body(tree)
    return body
