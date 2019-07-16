"""
run_init contains classes for initializing a run of a
double layered system

these classes are intended to be used by importing them into the simulation
system and using these to read in the run definition file

@author: Matthew D. Gerboth
@email: matthew.d.gerboth@vanderbilt.edu
@date: 2019
"""
import numpy as np

_N_SEGMENTS = 5
_N_NODES = 4

class Run:
    """Defines a run section of the simulation. i.e. a period of time in the
    simulation with a defined number of timesteps and a constant enviroment
    """
    def __init__(self, timesteps, condtemp, radtemp, humidity):
        self.timesteps = timesteps
        self.condtemp = condtemp
        self.radtemp = radtemp
        self.humidity = humidity


class SimulationParams:
    """SimulationParams holds the parameters for a simulation and constructs these 
    by reading in test files
    """
    def __init__(self, runfile):
        """Loads in a run file and build the simulation and build and
        stores the parameters

        """
        self.runs = []
        self.itemps = np.zeros(_N_SEGMENTS * _N_NODES) 
        with open(runfile) as f:
            f.readline()
            self.timestep = float(f.readline().strip())
            f.readline()
            for i in range(_N_SEGMENTS):
                segment_temps = [float(t) for t in f.readline().strip().split()]
                self.itemps[i*_N_NODES:i*_N_NODES + _N_NODES] = segment_temps
            f.readline()
            rline = f.readline()
            while rline:
                run_args = [float(a) for a in rline.strip().split()]
                self.runs.append(Run(*run_args))
                rline = f.readline()
