import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import re

class collison():
    def __init__(self, id, x, y, z, time, timestep):
        self.x = x
        self.y = y
        self.z = z
        self.time = time
        self.timestep = timestep
        self.id = id
        self.N = None

class Dump():
    def __init__(self, path):
        self.timestep, self.data = self.readfile(path)
        self.force_tolerance = 3.0
        self.collisions = []

    def readfile(self, path):
        file = open(path, 'r')
        timesteps = []
        cols = None
        ncols = None
        counter = 0
        nextTimestep = False
        nAtomsNext = False
        nAtoms = None
        propNext = False
        all_states = []
        id = None
        for line in file : 
            line = re.sub(r'\n', '', line)
            content = line.split(' ')
            if nextTimestep :
                timesteps += [int(line)]
                nextTimestep = False
            elif nAtomsNext :
                nAtoms = int(line)
                current_state = [None]*nAtoms
                nAtomsNext = False
            elif propNext :
                id = int(content[0]) - 1
                current_atom = [None]*ncols
                for i in range(ncols) :
                    current_atom[i] = float(content[i + 2])
                current_state[id] = current_atom
                counter += 1
                if counter >= nAtoms :
                    all_states.append(current_state)
                    propNext = False
            if content[0] == 'ITEM:' :
                if content[1] == 'TIMESTEP':
                    nextTimestep = True
                if content[1] == 'ATOMS' :
                    if not cols :
                        cols = content[4:]
                        ncols = len(cols)
                    counter = 0
                    propNext = True
                if content[1] == 'NUMBER' :
                    nAtomsNext = True 
        self.N = nAtoms
        return timesteps, all_states
          
    def pick_single_particle_prop(self, id, prop):
        prop_vec = []
        for i in range(len(self.timestep)):
            prop_vec += [self.data[i][id][prop]]
        return prop_vec
    
    def get_position(self, id, time_index):
        """
        Get position of atom 'id' at time = timestep[time_index].
        Returns: (x,y,z,t, time_index) of collision.
        """
        return (self.data[time_index][id][0], self.data[time_index][id][1],
                self.data[time_index][id][2], self.timestep[time_index], time_index)
    
    def get_single_profile(self, id, props):
        """
        prop: List of indices of x, y, and z component.
        id: Atom number

        Computes magnitute of vector as function of time, for
        a single atom.
        """
        res = 0
        for prop in props:
            res =+ np.array(self.pick_single_particle_prop(id, prop))**2
        return np.sqrt(res)

    def get_single_force_profile(self, id):
        """
        Single force profile
        """
        return self.get_single_profile(id, props = [6,7,8])
    
    def find_atom_collsions(self, id):
        """
        Returns time and locations of collisions of particles.
        """
        force = self.get_single_force_profile(id)
        colliding_force = force.copy()
        colliding_force[colliding_force < self.force_tolerance] = 0 
        peak_index, _ = find_peaks(colliding_force)
        positions = []
        for timestep in peak_index:
            positions += [self.get_position(id, timestep)]
        return force, colliding_force, positions
    
    def list_all_collisions(self):
        for id in range(self.N) :
            _, __, positions = self.find_atom_collsions(id)
            for position in positions :
                self.collisions += [collison(id, *position)]
        return
