import numpy as np
from numpy.linalg import norm
from atomset import atomset

# a function to search the neighboring in a system
def neighbor (Aset, Bset, length, limit, num_neigh):
    """
    Aset, Bset are two atom sets
    length is the cell length
    limit is the three bond length
    """
    # neighboring
    neigh = np.zeros((Aset.size(),num_neigh))
    for i in range(Aset.size()):
        pos_i = Aset.pos(i)
        neigh_seq = 0
        for j in range(Bset.size()):
            pos_j = Bset.pos(j)
            for k in range(Aset.ndim()):      # considering the periodic boundary condition
                pos_j[k] = pos_j[k] - int(round((pos_j[k]-pos_i[k])))
            dist = norm(np.array([length[0]*(pos_j[0]-pos_i[0]),length[1]*(pos_j[1]-pos_i[1]),length[2]*(pos_j[2]-pos_i[2])]))
            if (dist/limit < 1.2 and dist/limit > 0.8):
                neigh[i][neigh_seq] = j
                neigh_seq += 1
    #print (neigh)
    return neigh
# end def 

# a function to search all kinds of neighborings
def neigh_sys (Aset, Bset, Oset, length, limit_AA, limit_BB, limit_OO, limit_AB, limit_AO, limit_BO):
    """
    Aset, Bset, Oset are two atom sets
    length is the cell length
    limit is the three bond length
    """
    # for A
    num_neigh = Aset.size_neiA()
    Aset.change_all_neigh_A( neighbor(Aset, Aset, length, limit_AA, num_neigh) )
    num_neigh = Aset.size_neiB()
    Aset.change_all_neigh_B( neighbor(Aset, Bset, length, limit_AB, num_neigh) )
    num_neigh = Aset.size_neiO()
    Aset.change_all_neigh_O( neighbor(Aset, Oset, length, limit_AO, num_neigh) )
    # for B
    num_neigh = Bset.size_neiA()
    Bset.change_all_neigh_A( neighbor(Bset, Aset, length, limit_AB, num_neigh) )
    num_neigh = Bset.size_neiB()
    Bset.change_all_neigh_B( neighbor(Bset, Bset, length, limit_BB, num_neigh) )
    num_neigh = Bset.size_neiO()
    Bset.change_all_neigh_O( neighbor(Bset, Oset, length, limit_BO, num_neigh) )
    # for O
    num_neigh = Oset.size_neiA()
    Oset.change_all_neigh_A( neighbor(Oset, Aset, length, limit_AO, num_neigh) )
    num_neigh = Oset.size_neiB()
    Oset.change_all_neigh_B( neighbor(Oset, Bset, length, limit_BO, num_neigh) )
    num_neigh = Oset.size_neiO()
    Oset.change_all_neigh_O( neighbor(Oset, Oset, length, limit_OO, num_neigh) )
# end def 

def init_neighbor():
    # initialize the three atom set
    Aset = atomset(name="A")
    Bset = atomset(name="B")
    Oset = atomset(name="O")
    Aset.init_pos_3d()
    Bset.init_pos_3d()
    Oset.init_pos_3d()
    # end initialization
    length = np.array([10.86964544, 10.86964544, 15.372]) # cell length
    delta_AA = np.array([0.25, 0.25, 0.00])  # distance between two neighboring A
    delta_BB = np.array([0.25, 0.25, 0.00])  # distance between two neighboring B
    delta_OO = np.array([0.00, 0.25, 0.00])  # distance between two neighboring oxygen
    delta_AB = np.array([0.00, 0.25, 0.125])  # distance between A-B
    delta_AO = np.array([0.25, 0.00, 0.00])  # distance between A-O
    delta_BO = np.array([0.125, 0.125, 0.00])  # distance between B-O
    limit_AA = norm(np.array([length[0]*delta_AA[0],length[1]*delta_AA[1],length[2]*delta_AA[2]]))    # the limit length when considering wheter oxygen has moved to another position
    limit_BB = norm(np.array([length[0]*delta_BB[0],length[1]*delta_BB[1],length[2]*delta_BB[2]]))    
    limit_OO = norm(np.array([length[0]*delta_OO[0],length[1]*delta_OO[1],length[2]*delta_OO[2]]))    
    limit_AB = norm(np.array([length[0]*delta_AB[0],length[1]*delta_AB[1],length[2]*delta_AB[2]]))    
    limit_AO = norm(np.array([length[0]*delta_AO[0],length[1]*delta_AO[1],length[2]*delta_AO[2]]))
    limit_BO = norm(np.array([length[0]*delta_BO[0],length[1]*delta_BO[1],length[2]*delta_BO[2]]))
    # update neighboring information
    neigh_sys (Aset, Bset, Oset, length, limit_AA, limit_BB, limit_OO, limit_AB, limit_AO, limit_BO)
    Aset.neigh_Aout()
    Aset.neigh_Bout()
    Aset.neigh_Oout()
    Bset.neigh_Aout()
    Bset.neigh_Bout()
    Bset.neigh_Oout()
    Oset.neigh_Aout()
    Oset.neigh_Bout()
    Oset.neigh_Oout()
    return Aset, Bset, Oset
# end def

if __name__ == '__main__':

    init_neighbor()

# end __main__
