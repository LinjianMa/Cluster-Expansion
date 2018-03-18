import numpy as np
import itertools
from copy import deepcopy

class atomset:

    """ The atomset class is designed to hold the position, neighboring information of a set of particles. """

    def __init__(self,ndim=3,name="A"):
        """ create a particle set with num atoms in ndim dimensions. """
        # allocate local variables and arrays
        # particles still need to be initialized before the particle set can be used
        if (name == "A"):
            num = 32
            numnei_A = 6
            numnei_B = 8
            numnei_O = 12
        if (name == "B"):
            num = 32
            numnei_A = 8
            numnei_B = 6
            numnei_O = 6
        if (name == "O"):
            num = 96
            numnei_A = 4
            numnei_B = 2
            numnei_O = 8
        self.name  = name           # name of particle set, only used in __str__ 
        self._num  = num
        self._ndim = ndim
        self._numnei_A  =  numnei_A       # number of A site neighbors
        self._numnei_B  =  numnei_B       # number of B site neighbors
        self._numnei_O  =   numnei_O      # number of O site neighbors
        self._alloc_pos_neigh()           # allocate memory for position and neighboring  
        self.initialized = 0   # must call a position initialization routine before this object can be used in a simulation
    # end def

    def __str__(self):
        # allow atomset object to be printed
        description = "atomset " + self.name + " has %d particles \n" % self._num
        description += str(self._pos) 
        if not self.initialized:
            description += "\n WARNING: not initialized "
        # end if
        return description
    # end def

    def _alloc_pos_neigh(self):
        # allocate arrays for position, velocity and acceleration
        self._pos   = np.zeros((self._num,self._ndim))
        self._species   = np.ones(self._num)      # 1 or 0
        """
        As to A site, 0=La, 1=Sr
        As to B site, 0=Co, 1=Fe
        As to O site, 0=O , 1=V
        """
        self._neigh_A   = np.zeros((self._num,self._numnei_A))
        self._neigh_B   = np.zeros((self._num,self._numnei_B))
        self._neigh_O   = np.zeros((self._num,self._numnei_O))
    # end def

    # ---------------------------
    # begin accessor methods
    #  these are trivial accessor methods
    # Note: These methods can be generalized better, but I chose to exapand them for readability.

    def size(self):
        return self._num
    # end def size

    def size_neiA(self):
        return self._numnei_A
    # end def size

    def size_neiB(self):
        return self._numnei_B
    # end def size

    def size_neiO(self):
        return self._numnei_O
    # end def size

    def name(self):
        return self.name
    # end def size

    def ndim(self):
        return self._ndim
    # end def

    def all_pos(self):
        # return the positions of all the atoms
        return deepcopy( self._pos )
    # end def 
    def all_species(self):
        # return the species of all the atoms
        return deepcopy( self._species )
    # end def 
    def all_neigh_A(self):
        # return the A-site neighbor of all the atoms
        return deepcopy( self._neigh_A )
    # end def 
    def all_neigh_B(self):
        # return the B-site neighbor of all the atoms
        return deepcopy( self._neigh_B )
    # end def 
    def all_neigh_O(self):
        # return the O-site neighbor of all the atoms
        return deepcopy( self._neigh_O )
    # end def 

    def pos(self,iat):
        # return position of the iat th atom
        return self._pos[iat,:].copy()
    # end def
    def species(self,iat):
        # return species of the iat th atom
        return self._species[iat].copy()
    # end def
    def neigh_A(self,iat):
        # return A-site neighbor of the iat th atom
        return self._neigh_A[iat,:].copy()
    # end def
    def neigh_B(self,iat):
        # return B-site neighbor of the iat th atom
        return self._neigh_B[iat,:].copy()
    # end def
    def neigh_O(self,iat):
        # return O-site neighbor of the iat th atom
        return self._neigh_O[iat,:].copy()
    # end def

    # end accessor methods
    # ---------------------------

    # ---------------------------
    # begin modifier methods
    #  these are trivial modifier methods

    # use numpy's built-in functions for fast array assignment
    def change_all_pos(self,new_pos):
        # guard against misuse
        assert new_pos.shape == self._pos.shape, "shape mismatch"
        # renew all entries of self._pos fast
        np.copyto(self._pos,new_pos)
    # end def

    def change_all_species(self,new_species):
        assert new_species.shape == self._species.shape, "shape mismatch"
        np.copyto(self._species,new_species)
    # end def

    def change_all_neigh_A(self,new_all_neigh_A):
        assert new_all_neigh_A.shape == self._neigh_A.shape, "shape mismatch"
        np.copyto(self._neigh_A,new_all_neigh_A)
    # end def

    def change_all_neigh_B(self,new_all_neigh_B):
        assert new_all_neigh_B.shape == self._neigh_B.shape, "shape mismatch"
        np.copyto(self._neigh_B,new_all_neigh_B)
    # end def

    def change_all_neigh_O(self,new_all_neigh_O):
        assert new_all_neigh_O.shape == self._neigh_O.shape, "shape mismatch"
        np.copyto(self._neigh_O,new_all_neigh_O)
    # end def

    def change_pos(self,iat,pos1):
        assert pos1.shape==self._pos[iat,:].shape, "shape mismatch"
        self._pos[iat,:] = pos1
    # end def

    def change_species(self,iat,new_species):
        #assert new_species.shape==self._species[iat].shape, "shape mismatch"
        self._species[iat] = new_species
    # end def

    def change_neigh_A(self,iat,neiA_input):
        assert neiA_input.shape==self._neigh_A[iat,:].shape, "shape mismatch"
        self._neigh_A[iat,:] = neiA_input
    # end def

    def change_neigh_B(self,iat,neiB_input):
        assert neiB_input.shape==self._neigh_B[iat,:].shape, "shape mismatch"
        self._neigh_B[iat,:] = neiB_input
    # end def

    def change_neigh_O(self,iat,neiO_input):
        assert neiO_input.shape==self._neigh_O[iat,:].shape, "shape mismatch"
        self._neigh_O[iat,:] = neiO_input
    # end def

    # end modifier methods
    # ---------------------------

    # ---------------------------
    # begin initializer methods
    #  these functions do NOT require modification

    def init_pos_3d(self):
        # initialize the particles in uniform 3d supercell

        # generate integer positions for the atoms inside the cube
        new_pos = np.zeros((self._num,self._ndim))
        if (self.name == "A"):
            #La, Sr, 32 atoms in all
            for k in range (4):
                new_pos[8*k]   = [0.00, 0.00, float(k/4.0)]
                new_pos[8*k+1] = [0.00, 0.50, float(k/4.0)]
                new_pos[8*k+2] = [0.25, 0.25, float(k/4.0)]
                new_pos[8*k+3] = [0.25, 0.75, float(k/4.0)]
                new_pos[8*k+4] = [0.50, 0.00, float(k/4.0)]
                new_pos[8*k+5] = [0.50, 0.50, float(k/4.0)]
                new_pos[8*k+6] = [0.75, 0.25, float(k/4.0)]
                new_pos[8*k+7] = [0.75, 0.75, float(k/4.0)]
        if (self.name == "B"):
            #Co, Fe, 32 atoms in all
            for k in range (4):
                new_pos[8*k]   = [0.00, 0.25, float(k/4.0+0.125)]
                new_pos[8*k+1] = [0.00, 0.75, float(k/4.0+0.125)]
                new_pos[8*k+2] = [0.25, 0.00, float(k/4.0+0.125)]
                new_pos[8*k+3] = [0.25, 0.50, float(k/4.0+0.125)]
                new_pos[8*k+4] = [0.50, 0.25, float(k/4.0+0.125)]
                new_pos[8*k+5] = [0.50, 0.75, float(k/4.0+0.125)]
                new_pos[8*k+6] = [0.75, 0.00, float(k/4.0+0.125)]
                new_pos[8*k+7] = [0.75, 0.50, float(k/4.0+0.125)]
        if (self.name == "O"):
            #O, 96 atoms in all
            for k in range (4):
                new_pos[8*k]   = [0.00, 0.25, float(k/4.0)]
                new_pos[8*k+1] = [0.00, 0.75, float(k/4.0)]
                new_pos[8*k+2] = [0.25, 0.00, float(k/4.0)]
                new_pos[8*k+3] = [0.25, 0.50, float(k/4.0)]
                new_pos[8*k+4] = [0.50, 0.25, float(k/4.0)]
                new_pos[8*k+5] = [0.50, 0.75, float(k/4.0)]
                new_pos[8*k+6] = [0.75, 0.00, float(k/4.0)]
                new_pos[8*k+7] = [0.75, 0.50, float(k/4.0)]                
            for k in range (4):
                for j in range (4):
                    for i in range (4):
                        new_pos[32+k*16+j*4+i]   = [float(i/4.0+0.125), float(j/4.0+0.125), float(k/4.0+0.125)]
        # renew all entries of self._pos fast
        self.change_all_pos(new_pos)
        self.initialized = True
    # end def init_pos_3d
    
    # end initializer methods
    # ---------------------------

    # ---------------------------
    # begin writing functions

    # writing the positions
    def posout(self):
        outFile=open("position_"+str(self.name)+".xyz","w")
        for i in range(self._num):
            outFile.write(str(self._pos[i][0])+"        "+str(self._pos[i][1])+"        "+str(self._pos[i][2])+"\n")
        outFile.close()
    # end def posout
    
    # writing the A site neighboring number
    def neigh_Aout(self):
        outFile=open("neigh_A_"+str(self.name)+".xyz","w")
        for i in range(self._num):
            for j in range(self._numnei_A):
                outFile.write(str(self._neigh_A[i][j])+"        ")
            outFile.write("\n")
        outFile.close()
    # end def neigh_Aout

    # writing the B site neighboring number
    def neigh_Bout(self):
        outFile=open("neigh_B_"+str(self.name)+".xyz","w")
        for i in range(self._num):
            for j in range(self._numnei_B):
                outFile.write(str(self._neigh_B[i][j])+"        ")
            outFile.write("\n")
        outFile.close()
    # end def neigh_Bout

    # writing the O site neighboring number
    def neigh_Oout(self):
        outFile=open("neigh_O_"+str(self.name)+".xyz","w")
        for i in range(self._num):
            for j in range(self._numnei_O):
                outFile.write(str(self._neigh_O[i][j])+"        ")
            outFile.write("\n")
        outFile.close()
    # end def neigh_Oout

    # end writing functions
    # ---------------------------
 
# end class ParticleSet

if __name__ == '__main__':

    pset = atomset(name="B")
    pset.init_pos_3d()
    print (pset)
    pset.posout()
    pset.neigh_Aout()

# end __main__
