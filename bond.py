import numpy as np
from numpy.linalg import norm
from atomset import atomset
from neighbor import *
import math

def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))
def length(v):
    return math.sqrt(dotproduct(v, v))
def angle(v1, v2):
    a = dotproduct(v1, v2) / (length(v1) * length(v2))
    if (a>1):
        a = 1
    if (a<-1):
        a = -1
    return math.acos(a)/math.pi*180

def angle_calc(pos_1,pos_2,pos_3,lattice,dim):
    for m in range(dim):      # considering the periodic boundary condition
        pos_1[m] = pos_1[m] - int(round((pos_1[m]-pos_3[m])))
        pos_2[m] = pos_2[m] - int(round((pos_2[m]-pos_3[m])))
    vector1 = (pos_1-pos_3)
    vector2 = (pos_2-pos_3)
    vector1 = [vector1[0]*lattice[0],vector1[1]*lattice[1],vector1[2]*lattice[2]]
    vector2 = [vector2[0]*lattice[0],vector2[1]*lattice[1],vector2[2]*lattice[2]]
    return vector1, vector2
#end def

""" bond definition for dimers"""
"""---------------------------"""
def bond_AO(Aset):
    bond = np.zeros((Aset.size()*Aset.size_neiO(),2))
    ite = 0
    for i in range(Aset.size()):
        neighbor = Aset.neigh_O(i)
        for j in range(len(neighbor)): 
            bond[ite] = [i,neighbor[j]]
            ite += 1
    # write out to the file
    outFile=open("bond_AO"+".xyz","w")
    for i in range(len(bond)):
        for j in range(len(bond[i])):
            outFile.write(str(bond[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def bond_BO(Bset):
    bond = np.zeros((Bset.size()*Bset.size_neiO(),2))
    ite = 0
    for i in range(Bset.size()):
        neighbor = Bset.neigh_O(i)
        for j in range(len(neighbor)): 
            bond[ite] = [i,neighbor[j]]
            ite += 1
    # write out to the file
    outFile=open("bond_BO"+".xyz","w")
    for i in range(len(bond)):
        for j in range(len(bond[i])):
            outFile.write(str(bond[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def bond_AB(Aset):
    bond = np.zeros((Aset.size()*Aset.size_neiB(),2))
    ite = 0
    for i in range(Aset.size()):
        neighbor = Aset.neigh_B(i)
        for j in range(len(neighbor)): 
            bond[ite] = [i,neighbor[j]]
            ite += 1
    # write out to the file
    outFile=open("bond_AB"+".xyz","w")
    for i in range(len(bond)):
        for j in range(len(bond[i])):
            outFile.write(str(bond[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def bond_OO(Oset):
    bond = np.zeros((Oset.size()*Oset.size_neiO()/2,2))
    ite = 0
    for i in range(Oset.size()):
        neighbor = Oset.neigh_O(i)
        for j in range(len(neighbor)):
            if (neighbor[j]>i):
                bond[ite] = [i,neighbor[j]]
                ite += 1
    # write out to the file
    outFile=open("bond_OO"+".xyz","w")
    for i in range(len(bond)):
        for j in range(len(bond[i])):
            outFile.write(str(bond[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def bond_AA(Aset):
    bond = np.zeros((Aset.size()*Aset.size_neiA()/2,2))
    ite = 0
    for i in range(Aset.size()):
        neighbor = Aset.neigh_A(i)
        for j in range(len(neighbor)):
            if (neighbor[j]>i):
                bond[ite] = [i,neighbor[j]]
                ite += 1
    # write out to the file
    outFile=open("bond_AA"+".xyz","w")
    for i in range(len(bond)):
        for j in range(len(bond[i])):
            outFile.write(str(bond[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def bond_BB(Bset):
    bond = np.zeros((Bset.size()*Bset.size_neiB()/2,2))
    ite = 0
    for i in range(Bset.size()):
        neighbor = Bset.neigh_B(i)
        for j in range(len(neighbor)):
            if (neighbor[j]>i):
                bond[ite] = [i,neighbor[j]]
                ite += 1
    # write out to the file
    outFile=open("bond_BB"+".xyz","w")
    for i in range(len(bond)):
        for j in range(len(bond[i])):
            outFile.write(str(bond[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def
"""---------------------------"""
"""end bond definition for dimers"""

""" bond definition for trimers"""
"""---------------------------"""
def trimer_AOA(Oset,Aset,lattice):
    trimer_90 = np.zeros((Oset.size()*4,3))
    trimer_00 = np.zeros((Oset.size()*2,3))
    ite_00 = 0
    ite_90 = 0
    for i in range(Oset.size()):
        neighbor = Oset.neigh_A(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                pos_j = Aset.pos(int(neighbor[j]))
                pos_k = Aset.pos(int(neighbor[k]))
                pos_O = Oset.pos(i)
                vector1, vector2 = angle_calc(pos_j,pos_k,pos_O,lattice,Aset.ndim())
                if (angle(vector1, vector2)>85 and angle(vector1, vector2)<95):
                    trimer_90[ite_90] = [neighbor[j],i,neighbor[k]]
                    ite_90 += 1
                else:
                    trimer_00[ite_00] = [neighbor[j],i,neighbor[k]]
                    ite_00 += 1
    # write out to the file for trimer_00
    outFile=open("trimer_AOA_00"+".xyz","w")
    for i in range(len(trimer_00)):
        for j in range(len(trimer_00[i])):
            outFile.write(str(trimer_00[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
    # write out to the file for trimer_90
    outFile=open("trimer_AOA_90"+".xyz","w")
    for i in range(len(trimer_90)):
        for j in range(len(trimer_90[i])):
            outFile.write(str(trimer_90[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def trimer_BOB(Oset,Bset,lattice):
    trimer = np.zeros((Oset.size(),3))
    ite = 0
    for i in range(Oset.size()):
        neighbor = Oset.neigh_B(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                trimer[ite] = [neighbor[j],i,neighbor[k]]
                ite += 1
    # write out to the file for trimer_00
    outFile=open("trimer_BOB_00"+".xyz","w")
    for i in range(len(trimer)):
        for j in range(len(trimer[i])):
            outFile.write(str(trimer[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def trimer_AOB(Oset,Aset,Bset,lattice):
    trimer = np.zeros((Oset.size()*8,3))
    ite = 0
    for i in range(Oset.size()):
        neighbor_A = Oset.neigh_A(i)
        neighbor_B = Oset.neigh_B(i)
        for j in range(len(neighbor_A)):
            for k in range(len(neighbor_B)):
                trimer[ite] = [neighbor_A[j],i,neighbor_B[k]]
                ite += 1
    # write out to the file for trimer_00
    outFile=open("trimer_AOB"+".xyz","w")
    for i in range(len(trimer)):
        for j in range(len(trimer[i])):
            outFile.write(str(trimer[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def trimer_OBO(Oset,Bset,lattice):
    trimer_00 = np.zeros((Bset.size()*3,3))
    trimer_90 = np.zeros((Bset.size()*12,3))
    ite_00 = 0
    ite_90 = 0
    for i in range(Bset.size()):
        neighbor = Bset.neigh_O(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                pos_j = Oset.pos(int(neighbor[j]))
                pos_k = Oset.pos(int(neighbor[k]))
                pos_B = Bset.pos(i)
                vector1, vector2 = angle_calc(pos_j,pos_k,pos_B,lattice,Bset.ndim())
                if (angle(vector1, vector2)>85 and angle(vector1, vector2)<95):
                    trimer_90[ite_90] = [neighbor[j],i,neighbor[k]]
                    ite_90 += 1
                else:
                    trimer_00[ite_00] = [neighbor[j],i,neighbor[k]]
                    ite_00 += 1
    # write out to the file for trimer_00
    outFile=open("trimer_OBO_00"+".xyz","w")
    for i in range(len(trimer_00)):
        for j in range(len(trimer_00[i])):
            outFile.write(str(trimer_00[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
    # write out to the file for trimer_90
    outFile=open("trimer_OBO_90"+".xyz","w")
    for i in range(len(trimer_90)):
        for j in range(len(trimer_90[i])):
            outFile.write(str(trimer_90[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def trimer_OOO_90(Oset,lattice):
    trimer_90 = np.zeros((Oset.size()*8,3))
    ite_90 = 0
    for i in range(Oset.size()):
        neighbor = Oset.neigh_O(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                pos_j = Oset.pos(int(neighbor[j]))
                pos_k = Oset.pos(int(neighbor[k]))
                pos_O = Oset.pos(i)
                vector1, vector2 = angle_calc(pos_j,pos_k,pos_O,lattice,Oset.ndim())
                if (angle(vector1, vector2)>85 and angle(vector1, vector2)<95):
                    trimer_90[ite_90] = [neighbor[j],i,neighbor[k]]
                    ite_90 += 1
    # write out to the file for trimer_00
    outFile=open("trimer_OOO_90"+".xyz","w")
    for i in range(len(trimer_90)):
        for j in range(len(trimer_90[i])):
            outFile.write(str(trimer_90[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def trimer_OOO_60(Oset,lattice):
    trimer_60 = np.zeros((512,3))
    ite_60 = 0
    for i in range(Oset.size()):
        for j in range(i+1,Oset.size()):
            for k in range(j+1,Oset.size()):
                pos_i = Oset.pos(i)
                pos_j = Oset.pos(j)
                pos_k = Oset.pos(k)
                vector1, vector2 = angle_calc(pos_i,pos_j,pos_k,lattice,Oset.ndim())
                vector3, vector4 = angle_calc(pos_i,pos_k,pos_j,lattice,Oset.ndim())
                vector5, vector6 = angle_calc(pos_j,pos_k,pos_i,lattice,Oset.ndim())
                if (angle(vector1, vector2)>55 and angle(vector1, vector2)<65 and angle(vector3, vector4)>55 and angle(vector3, vector4)<65 and angle(vector5, vector6)>55 and angle(vector5, vector6)<65):
                    trimer_60[ite_60] = [i,j,k]
                    ite_60 += 1
    # write out to the file for trimer_60
    outFile=open("trimer_OOO_60"+".xyz","w")
    for i in range(len(trimer_60)):
        for j in range(len(trimer_60[i])):
            outFile.write(str(trimer_60[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def trimer_BBB_90(Bset,lattice):
    trimer_90 = np.zeros((Bset.size()*12,3))
    ite_90 = 0
    for i in range(Bset.size()):
        neighbor = Bset.neigh_B(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                pos_j = Bset.pos(int(neighbor[j]))
                pos_k = Bset.pos(int(neighbor[k]))
                pos_B = Bset.pos(i)
                vector1, vector2 = angle_calc(pos_j,pos_k,pos_B,lattice,Bset.ndim())
                if (angle(vector1, vector2)>85 and angle(vector1, vector2)<95):
                    trimer_90[ite_90] = [neighbor[j],i,neighbor[k]]
                    ite_90 += 1
    # write out to the file for trimer_00
    outFile=open("trimer_BBB_90"+".xyz","w")
    for i in range(len(trimer_90)):
        for j in range(len(trimer_90[i])):
            outFile.write(str(trimer_90[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def trimer_BBB_180(Bset,lattice):
    trimer_180 = np.zeros((Bset.size()*3,3))
    ite_180 = 0
    for i in range(Bset.size()):
        neighbor = Bset.neigh_B(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                pos_j = Bset.pos(int(neighbor[j]))
                pos_k = Bset.pos(int(neighbor[k]))
                pos_B = Bset.pos(i)
                vector1, vector2 = angle_calc(pos_j,pos_k,pos_B,lattice,Bset.ndim())
                if (angle(vector1, vector2)>175 and angle(vector1, vector2)<185):
                    trimer_180[ite_180] = [neighbor[j],i,neighbor[k]]
                    ite_180 += 1
    # write out to the file for trimer_00
    outFile=open("trimer_BBB_180"+".xyz","w")
    for i in range(len(trimer_180)):
        for j in range(len(trimer_180[i])):
            outFile.write(str(trimer_180[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def trimer_AAA_90(Aset,lattice):
    trimer_90 = np.zeros((Aset.size()*12,3))
    ite_90 = 0
    for i in range(Aset.size()):
        neighbor = Aset.neigh_A(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                pos_j = Aset.pos(int(neighbor[j]))
                pos_k = Aset.pos(int(neighbor[k]))
                pos_A = Aset.pos(i)
                vector1, vector2 = angle_calc(pos_j,pos_k,pos_A,lattice,Aset.ndim())
                if (angle(vector1, vector2)>85 and angle(vector1, vector2)<95):
                    trimer_90[ite_90] = [neighbor[j],i,neighbor[k]]
                    ite_90 += 1
    # write out to the file for trimer_00
    outFile=open("trimer_AAA_90"+".xyz","w")
    for i in range(len(trimer_90)):
        for j in range(len(trimer_90[i])):
            outFile.write(str(trimer_90[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def trimer_ABA_90(Aset,Bset,lattice):
    trimer_90 = np.zeros((Bset.size()*12,3))
    ite_90 = 0
    for i in range(Bset.size()):
        neighbor = Bset.neigh_A(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                pos_j = Aset.pos(int(neighbor[j]))
                pos_k = Aset.pos(int(neighbor[k]))
                pos_B = Bset.pos(i)
                vector1, vector2 = angle_calc(pos_j,pos_k,pos_B,lattice,Aset.ndim())
                if (angle(vector1, vector2)<106):
                    trimer_90[ite_90] = [neighbor[j],i,neighbor[k]]
                    ite_90 += 1
    # write out to the file for trimer_00
    outFile=open("trimer_ABA_90"+".xyz","w")
    for i in range(len(trimer_90)):
        for j in range(len(trimer_90[i])):
            outFile.write(str(trimer_90[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def trimer_BAB_90(Aset,Bset,lattice):
    trimer_90 = np.zeros((Aset.size()*12,3))
    ite_90 = 0
    for i in range(Aset.size()):
        neighbor = Aset.neigh_B(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                pos_j = Bset.pos(int(neighbor[j]))
                pos_k = Bset.pos(int(neighbor[k]))
                pos_A = Aset.pos(i)
                vector1, vector2 = angle_calc(pos_j,pos_k,pos_A,lattice,Aset.ndim())
                if (angle(vector1, vector2)<90):
                    trimer_90[ite_90] = [neighbor[j],i,neighbor[k]]
                    ite_90 += 1
    # write out to the file for trimer_00
    outFile=open("trimer_BAB_90"+".xyz","w")
    for i in range(len(trimer_90)):
        for j in range(len(trimer_90[i])):
            outFile.write(str(trimer_90[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

"""---------------------------"""
"""end bond definition for trimers"""

"""bond definition for fourth cluster"""
"""---------------------------"""
def fourth_OOOO_90(Oset,lattice):
    fourth_90 = np.zeros((Oset.size()*8,4))
    ite_90 = 0
    for i in range(Oset.size()):
        neighbor = Oset.neigh_O(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                pos_j = Oset.pos(int(neighbor[j]))
                pos_k = Oset.pos(int(neighbor[k]))
                pos_O = Oset.pos(i)
                vector1, vector2 = angle_calc(pos_j,pos_k,pos_O,lattice,Oset.ndim())
                if (angle(vector1, vector2)>85 and angle(vector1, vector2)<95):
                    # search for the fourth O atom
                    neighbor_neighj = Oset.neigh_O(int(neighbor[j]))
                    neighbor_neighk = Oset.neigh_O(int(neighbor[k]))
                    for m in range(len(neighbor_neighj)):
                        if (int(neighbor_neighj[m])!=i and int(neighbor_neighj[m])!=int(neighbor[k]) and (neighbor_neighj[m] in neighbor_neighk) and (neighbor_neighj[m] not in neighbor)):
                            fourth_90[ite_90] = [neighbor[j],i,neighbor[k],neighbor_neighj[m]]
                            fourth_90[ite_90].sort()
                            ite_90 += 1
                            if (ite_90>1):
                                for ii in range(ite_90-1):
                                    if ((int(fourth_90[ite_90-1][0])==int(fourth_90[ii][0])) and (int(fourth_90[ite_90-1][1])==int(fourth_90[ii][1])) and (int(fourth_90[ite_90-1][2])==int(fourth_90[ii][2])) and (int(fourth_90[ite_90-1][3])==int(fourth_90[ii][3]))):
                                        ite_90 -= 1
                                        break
                                
    # write out to the file for trimer_00
    outFile=open("fourth_OOOO_90"+".xyz","w")
    for i in range(len(fourth_90)):
        for j in range(len(fourth_90[i])):
            outFile.write(str(fourth_90[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def fourth_OOOO_60(Oset,lattice):
    fourth_60 = np.zeros((Oset.size()*8,4))
    ite_60 = 0
    for i in range(Oset.size()):
        neighbor = Oset.neigh_O(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                pos_j = Oset.pos(int(neighbor[j]))
                pos_k = Oset.pos(int(neighbor[k]))
                pos_O = Oset.pos(i)
                vector1, vector2 = angle_calc(pos_j,pos_k,pos_O,lattice,Oset.ndim())
                if (angle(vector1, vector2)>85 and angle(vector1, vector2)<95):
                    # search for the fourth O atom
                    neighbor_neighj = Oset.neigh_O(int(neighbor[j]))
                    neighbor_neighk = Oset.neigh_O(int(neighbor[k]))
                    for m in range(len(neighbor_neighj)):
                        if (int(neighbor_neighj[m])!=i and int(neighbor_neighj[m])!=int(neighbor[k]) and (neighbor_neighj[m] in neighbor_neighk) and (neighbor_neighj[m] in neighbor)):
                            fourth_60[ite_60] = [min(neighbor[j],neighbor[k]),min(i,neighbor_neighj[m]),max(i,neighbor_neighj[m]),max(neighbor[j],neighbor[k])]
                            ite_60 += 1
                            if (ite_60>1):
                                for ii in range(ite_60-1):
                                    if ((int(fourth_60[ite_60-1][0])==int(fourth_60[ii][0])) and (int(fourth_60[ite_60-1][1])==int(fourth_60[ii][1])) and (int(fourth_60[ite_60-1][2])==int(fourth_60[ii][2])) and (int(fourth_60[ite_60-1][3])==int(fourth_60[ii][3]))):
                                        ite_60 -= 1
                                        break
                                
    # write out to the file for trimer_00
    outFile=open("fourth_OOOO_60"+".xyz","w")
    for i in range(len(fourth_60)):
        for j in range(len(fourth_60[i])):
            outFile.write(str(fourth_60[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def

def fourth_BOOO_60(Bset,Oset,lattice):
    fourth_60 = np.zeros((Bset.size()*8,4))
    ite_60 = 0
    for i in range(Bset.size()):
        neighbor = Bset.neigh_O(i)
        for j in range(len(neighbor)):
            for k in range(j+1,len(neighbor)):
                for m in range(k+1,len(neighbor)):
                    pos_j = Oset.pos(int(neighbor[j]))
                    pos_k = Oset.pos(int(neighbor[k]))
                    pos_m = Oset.pos(int(neighbor[m]))
                    pos_B = Bset.pos(i)
                    vector1, vector2 = angle_calc(pos_j,pos_k,pos_m,lattice,Oset.ndim())
                    if (angle(vector1, vector2)>55 and angle(vector1, vector2)<65):
                        fourth_60[ite_60] = [i, neighbor[j], neighbor[k], neighbor[m]]
                        ite_60 += 1                         
    # write out to the file for trimer_00
    outFile=open("fourth_BOOO_60"+".xyz","w")
    for i in range(len(fourth_60)):
        for j in range(len(fourth_60[i])):
            outFile.write(str(fourth_60[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end def
"""---------------------------"""
"""bond definition for fourth cluster"""

if __name__ == '__main__':

    Aset, Bset, Oset = init_neighbor()
    lattice = np.array([10.86964544, 10.86964544, 15.372]) # cell length
    #bond_AO(Aset)
    #bond_BO(Bset)
    #bond_OO(Oset)
    #bond_AB(Aset)
    #bond_AA(Aset)
    #bond_BB(Bset)
    #trimer_AOA(Oset,Aset,lattice)
    #trimer_BOB(Oset,Bset,lattice)
    #trimer_AOB(Oset,Aset,Bset,lattice)
    #trimer_OBO(Oset,Bset,lattice)
    #trimer_OOO_60(Oset,lattice)
    #trimer_OOO_90(Oset,lattice)
    #fourth_OOOO_90(Oset,lattice)
    #fourth_OOOO_60(Oset,lattice)
    #trimer_BBB_90(Bset,lattice)
    #trimer_AAA_90(Aset,lattice)
    #trimer_ABA_90(Aset,Bset,lattice)
    #trimer_BAB_90(Aset,Bset,lattice)
    #trimer_BBB_180(Bset,lattice)
    fourth_BOOO_60(Bset,Oset,lattice)
# end __main__
