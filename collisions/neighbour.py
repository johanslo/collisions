import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import copy

class Pair():
    def __init__(self, atom1, atom2, pairforce, isIncreasing):
        self.atoms = (atom1, atom2)
        self.pairforce = pairforce
        self.isIncreasing = isIncreasing

    def __eq__(self, otherPair):
        """
        Overloading == operator to check if two pairs cosists of the two same atoms.
        Inverting the order of the IDs will not affect the results.
        """
        if self.atoms[0] == otherPair.atoms[0] and self.atoms[1] == otherPair.atoms[1] :
            return True
        elif self.atoms[0] == otherPair.atoms[1] and self.atoms[1] == otherPair.atoms[0] :
            return True
        else:
            return False
        
    def __lt__(self, otherPair):
        """
        Checks if exactly one atom is common between two pairs
        """
        if self.atoms[0] == otherPair.atoms[1] or self.atoms[1] == otherPair.atoms[0] or \
           self.atoms[1] == otherPair.atoms[1] or self.atoms[0] == otherPair.atoms[0]:
            return True
        else:
            return False
        
    def checkForCollision(self, previousPair) :
        """
        Checks if pair force has reached a local maxima, i.e. pair-force has gone from increasing
        to decreasing when compared with the previous time-step. Additionally, the force must
        be positive (repulsive) to be considered 
        Updates the 'isIncreasing' member variable.
        """
        self.isIncreasing = self.pairforce > previousPair.pairforce
        if not self.isIncreasing and previousPair.isIncreasing and self.pairforce > 0 :
            return True
        else :
            return False
        
    def __str__(self):
        """
        Overloading str-opertator to make the Pair object compatible with print().
        Prints: [Atoms : #atoms, Pair-force : #pairforce, Increasing : #isIncreasing]
        """
        return f'[Atoms : ({self.atoms[0]}, {self.atoms[1]}), Pair-force : {self.pairforce}, Increasing : {self.isIncreasing}]'
        
class Pairlist():
    def __init__(self):
        self.pairlist = []
    
    def addToPairlist(self, pair):
        """
        Adds a Pair-type object into the pairlist member variable.
        """
        self.pairlist.append(pair)
        return
    
    def getFromPairlist(self, pair):
        """
        Checks pairlist for a specified Pair-type object. 
        Returns  the pair + True if it is in the list.
        Returns None + False if not.
        """
        inList = False
        res = None
        for p in self.pairlist :
            if p == pair :
                inList = True
                res = p
                break
        return res, inList
    
    def __str__(self):
        """
        Makes the pair-list compatible with the print()-function.
        Prints each Pair seperated by a newline.
        """
        res = ''
        for i in range(len(self.pairlist)):
            res += f'{i} : {self.pairlist[i]}'
            res += '\n'
        return res

class Position():
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def distanceTo(self, otherPosition, boxDimensions = None):
        """
        Distance between two positions. 
        TBA: taking account of periodic boundary conditions (dx=min(|x1-x2|,|x1+Lx-x2|, |x1-x2-Lx|))
        """
        if boxDimensions :
            raise NotImplementedError('Periodic boundaries not implemented in Position.distanceTo()')
        else:
            dx = self.x - otherPosition.x
            dy = self.y - otherPosition.y
            dz = self.z - otherPosition.z
            return np.sqrt(dx**2 + dy**2 + dz**2)
    def averagePosition(self, otherPosition, boxDimensions):
        """
        Computing average position between points. For equal one-component systems it
        is equivalent to the mass-center position.
        Returns: Position-object 
        TBA: taking into account 
        """
        if boxDimensions :
            raise NotImplementedError('Periodic boundaries not implemented in Position.distanceTo()')
        else:
            x_new = (self.x + otherPosition.x) / 2
            y_new = (self.y + otherPosition.y) / 2
            z_new = (self.z + otherPosition.z) / 2
            return Position(x_new, y_new, z_new)    

class Collision():
    def __init__(self, pair, timestep, order):
        self.pair = pair
        self.timestep = timestep
        self.maxforce = self.pair.pairforce
        self.order = order

    def __str__(self):
        return f'Pair : ({self.pair.atoms[0]}, {self.pair.atoms[1]}), \
Time : {self.timestep}, Force : {self.maxforce}, Order : {self.order}'

class NeighbourList():
    def __init__(self, path, f_index):
        """
        Reads the dump-file step-by-step, checking for collisions 

        Member variables:
        f_index (int): Index of force in the dump-file
        currentPairlist (Pairlist): List of interracting pairs in current timestep
        previousPairlist (Pairlist): Pairlist from previous timestep
        collisionList (list of Collision-object): All collisions that occur, 
                                                  stored in Collision-object
        timesteps (list of int): List of timestep-number 
        """
        self.f_index = f_index
        self.previousPairlist = None
        self.currentTimestep = 0
        self.currentPairlist = Pairlist()
        self.allPairlists = []
        self.collisionList = []
        self.timesteps = self.readfile(path)

    def readfile(self, path):
        """
        Reads dump-file step by step, storing the states/pairlists, and 
        checks for collisisions.
        """
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
                self.currentTimestep = int(line) 
                timesteps += [self.currentTimestep]
                nextTimestep = False
            elif nAtomsNext :
                nAtoms = int(line)
                current_state = [None]*nAtoms
                nAtomsNext = False
            elif propNext :
                # id1 = int(content[1])
                # id2 = int(content[2]) 
                # current_atom = [None]*ncols
                # current_atom[0] = id1
                # current_atom[1] = id2
                # for i in range(ncols-2) :
                #     current_atom[i+2] = float(content[i + 3])
                # current_state[counter] = current_atom
                # counter += 1
                # if counter >= nAtoms :
                #     all_states.append(current_state)
                #     propNext = False
                current_pair = Pair(int(content[1]), int(content[2]), float(content[self.f_index]), False)
                if self.previousPairlist :
                        self.checkIncreasingAndCollisions(current_pair)
                self.currentPairlist.addToPairlist(current_pair)
                counter += 1
                if counter >= nAtoms :
                    self.allPairlists.append(self.currentPairlist)
                    self.previousPairlist = copy.deepcopy(self.currentPairlist)
                    self.currentPairlist = Pairlist()
                    propNext = False
            if content[0] == 'ITEM:' :
                if content[1] == 'TIMESTEP':
                    nextTimestep = True
                if content[1] == 'ENTRIES' :
                    if not cols :
                        cols = content[3:]
                        ncols = len(cols)
                    counter = 0
                    propNext = True
                if content[1] == 'NUMBER' :
                    nAtomsNext = True 
        self.N = nAtoms
        return timesteps
    
    def checkIncreasingAndCollisions(self, pair):
        """
        Checks if a pair is in the previous pairlist. If yes, the isIncreasing parameter is updated,
        and it is checked whether a collision has occured.
        """
        if not self.previousPairlist :
            raise RuntimeError('checkIncreasingAndCollision() did not find a previous pair list.')
        prevPair, inPrevList = self.previousPairlist.getFromPairlist(pair)
        if inPrevList :
            hasCollided = pair.checkForCollision(prevPair)
            if hasCollided :
                order = self.checkCollisionOrder(pair)
                self.collisionList.append(Collision(pair, self.currentTimestep, order))
        return
    
    def checkCollisionOrder(self, pair):
        """
        Checks number of atoms within cutoff-range of either atom in the colliding pair.
        This includes both repulsive and attractive interactions.
        Returns: order (integer) >= 2.
        """
        order = 2
        for otherPair in self.currentPairlist.pairlist :
            if pair < otherPair :
                order += 1
        return order 
    
    def extractEvent(self, pair):
        """
        Looks through pairlists for a specified pair to retrive their collision event.
        Returns: t (timestep-array of collison), f (force-array of forces during collision)
        """
        t = []
        f = []
        for i in range(len(self.allPairlists)) :
            found_force = False
            current_list = self.allPairlists[i] 
            t += [self.timesteps[i]]
            for j in range(len(current_list.pairlist)) :
                if pair == current_list.pairlist[j]:
                    f += [current_list.pairlist[j].pairforce]
                    found_force = True
            if not found_force :
                f += [0.0]
        return np.array(t), np.array(f)


# Note: modify pair/local in LAMMPS to obtain the center-of-mass position of 
#       the two particles of a pair. 
#
#Note 2: It should be decided if Neighbourlist object stores a list of all states, or only two
#         objects at the same time (which uses less memory and is better for large files) 
#
#TO DO :
#       Gi Pair-objekt to Position memdlemsvariabler
#       Angi posisjon til Collision-object ved Position.averagePosition()
#       Histogram-funksjon for å se fordeling av kollisjonsposisjon
#       Mean-free-path utregner
#       Plottefunksjoner (i annen fil)
#       Egenskapfunksjoner flyttes til annen fil ?