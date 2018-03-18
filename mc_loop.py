#!/usr/bin/env python
import numpy as np
import os, itertools, random, math, sys, errno, time
from numpy.linalg import norm
from atomset import atomset
from neighbor import *
from bond import *
from read_bond import *
from struc_init import *
from least_square_mod import *
import pylab   

def compute_energy(mat_cluster,V):
    return np.dot(mat_cluster,V)
# end def

def update_pos(Aset, Bset, Oset, ran_A, ran_B, ran_O):
    # Aset
    A_species = Aset.species(ran_A[0])
    Aset.change_species(ran_A[0], Aset.species(ran_A[1]))
    Aset.change_species(ran_A[1], A_species)
    # Bset
    B_species = Bset.species(ran_B[0]) 
    Bset.change_species(ran_B[0], Bset.species(ran_B[1]))
    Bset.change_species(ran_B[1], B_species)
    # Oset
    O_species = Oset.species(ran_O[0]) 
    Oset.change_species(ran_O[0], Oset.species(ran_O[1]))
    Oset.change_species(ran_O[1], O_species)
# end def update_pos

def mc_loop( temperature, nsweeps, ndump, restart, seed, num_lattice, lattice):
    """ perform Monte Carlo simulation, file checks should be done outside of this function """

    beta = 1./(0.0257/298*temperature) # inverse temperature
    np.random.seed(seed)  # seed random number generator

    # initialize atomset, no need for velocities
    # ---------------------------
    Aset, Bset, Oset = init_neighbor()
    if restart:  
        """
        function of read_poscar:
        1. set the species of each position: 0 is La, Co, O, 1 is Sr, F
        2. determine how many atoms each species have
        """
        num_atom = read_poscar(Aset, Bset, Oset)  
    else: 
        write_coords(lattice[0], lattice[1], lattice[2], num_lattice[0], num_lattice[1], num_lattice[2], num_lattice[3], num_lattice[4])
        num_atom = read_poscar(Aset, Bset, Oset)
    # Load all the bond relations
    lines = load_bond()
    # Calculate the V_parameter and RMSE
    V_parameter, RMSE = CEV_calc(True)
    # end initialize atomset
   
    # setup outputs
    # ---------------------------
    print( "{0:10s}  {1:15s}  {2:15s}".format(
           "Step","potential","Acceptance") )
    print_fmt = "{isweep:4d}  {potential:15f}  {acceptance:15f}"
    # end setup outputs

    # MC loop
    # ---------------------------
    mat_cluster    = calc_cluster(Aset, Bset, Oset, num_atom, lines)
    num_cluster    = 24
    potential      = compute_energy(mat_cluster[:,:num_cluster],V_parameter) # Calculate energy of entire system 
    move           = 0
    num_accepted   = 0  
    energy         = []
    acc_rate       = 0

    for idump in range(nsweeps/ndump):
        for isweep in range(ndump):
            # choose one particle in Aset, Bset, Oset
            ran_A = np.random.choice(Aset.size(), 2,  replace=False)
            ran_B = np.random.choice(Bset.size(), 2,  replace=False)
            ran_O = np.random.choice(Oset.size(), 2,  replace=False)

            # change species
            update_pos(Aset, Bset, Oset, ran_A, ran_B, ran_O)

            # Calculate energy of entire system 
            potential_new = compute_energy(mat_cluster[:,:num_cluster],V_parameter)

            # calculate acceptance rate
            acc_rate = min(1,np.exp(-beta*(potential_new-potential)))
            # accept or decline
            if (np.random.rand() < acc_rate):       # accepted
                potential = potential_new
                num_accepted += 1
            else:                                   # declined
                update_pos(Aset, Bset, Oset, ran_A, ran_B, ran_O)

            energy.append(potential)

        # save snapshot and report
        #pset.append_frame()
        #print( print_fmt.format(**{'isweep':(isweep+idump*ndump),'potential':potential,'acceptance':acc_rate}))
        # end if

    # end isweep
    print("Acceptance Rate:", float(num_accepted)/nsweeps )

    sweep_plot = range(0,nsweeps)
    sweep_plot  = [i for i in sweep_plot]
    pylab.plot(sweep_plot,energy,'r',linewidth=1.5)
    pylab.xlabel('Monte Carlo steps')
    pylab.ylabel('energy')
    pylab.title('Potential Energy Trace')
    pylab.show()


# end def mc_loop

if __name__=="__main__":

    # !!!! change seed when collecting statistics (pass a different -s)
    seed        = 1
    temperature = 1000          # temperature in reduced units
    nsweeps     = 1000          # total number of MC steps
    restart     = 1             # restart particle positions from trajectory file
    ndump       = 10            # interval to dump particle positions and report
    trajfile    = 'POSCAR'      # trajectory file
    num_lattice = np.array([10,22,10,22,95])  # number of La, Sr, Co, Fe, O
    lattice     = np.array([10.86964544, 10.86964544, 7.686*2])

    # check input & output locations
    # ==============================
    if restart: # restart from trajfile
        # check trajectory file
        if not os.path.exists(trajfile):
            print( "WARNING: no trajectory file found. Starting from scratch" )
            restart = False
        # end if
    # end if restart

    # start MC
    # ==============================
    mc_loop(
        temperature = temperature,
        box_length  = box_length,
        nsweeps     = nsweeps,
        ndump       = ndump,
        restart     = restart,
        seed        = seed,
        num_lattice = num_lattice,
        lattice     = lattice,
    )

# end __main__