import numpy as np
from numpy.linalg import norm
from atomset import atomset
from neighbor import *
from bond import *
import math, os, sys, random, errno, time

def silentchdir(filename):
    output = True     # Flag of chdir
    try:
        os.chdir(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        output = False
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred
    return output

def SubDirPath (d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])

def set_species(Aset, Bset, Oset, lines_poscar, num_atom):
    num_sum = np.zeros(5) # number of La, Sr, Co, Fe, O
    num_sum[0] = num_atom[0]
    for i in range(1, len(num_atom)):
        num_sum[i] = num_sum[i-1] + num_atom[i]
    # Set all as 1 very important !!!!!!!!!
    Aset.change_all_species(np.ones(Aset.size()))
    Bset.change_all_species(np.ones(Bset.size()))
    Oset.change_all_species(np.ones(Oset.size()))
    # Species of La
    for i in range(9,int(9+num_sum[0])):
        index = float(lines_poscar[i][2])/0.25*8+float(lines_poscar[i][0])/0.25*2
        if (float(lines_poscar[i][1]) == 0.75 or float(lines_poscar[i][1]) == 0.5):
            index += 1
        Aset.change_species(int(index),0)  
    # Species of Co
    for i in range(int(9+num_sum[1]),int(9+num_sum[2])):
        index = (float(lines_poscar[i][2])-0.125)/0.25*8+float(lines_poscar[i][0])/0.25*2
        if (float(lines_poscar[i][1]) == 0.75 or float(lines_poscar[i][1]) == 0.5):
            index += 1
        Bset.change_species(int(index),0)         
    # Species of O
    for i in range(int(9+num_sum[3]),int(9+num_sum[4])):
        if (float(lines_poscar[i][2]) % 0.25 == 0):
            index = float(lines_poscar[i][2])/0.25*8+float(lines_poscar[i][0])/0.25*2
            if (float(lines_poscar[i][1]) == 0.75 or float(lines_poscar[i][1]) == 0.5):
                index += 1     
        else:
            index = int((float(lines_poscar[i][2])-0.125)/0.25*16)+int((float(lines_poscar[i][1])-0.125)/0.25*4)+int((float(lines_poscar[i][0])-0.125)/0.25)+32
        Oset.change_species(int(index),0)  
#end def

def read_poscar (Aset, Bset, Oset):
    file = open('POSCAR','r')
    lines_poscar = [line.split() for line in file]
    num_atom = np.zeros(5) # number of La, Sr, Co, Fe, O
    for i in range(len(lines_poscar[6])):
        if (str(lines_poscar[5][i]) == "La"):
            num_atom[0] = int(lines_poscar[6][i])
        if (str(lines_poscar[5][i]) == "Sr"):
            num_atom[1] = int(lines_poscar[6][i])
        if (str(lines_poscar[5][i]) == "Co"):
            num_atom[2] = int(lines_poscar[6][i])
        if (str(lines_poscar[5][i]) == "Fe"):
            num_atom[3] = int(lines_poscar[6][i])
        if (str(lines_poscar[5][i]) == "O"):
            num_atom[4] = int(lines_poscar[6][i])
    set_species(Aset, Bset, Oset, lines_poscar, num_atom)
    return num_atom
#end def

def load_bond():
    lines = []
    # read bond_AA file
    file = open('bond_AA.xyz','r')     # 7. La-La    8. La-Sr    9. Sr-Sr
    lines_AA = [[int(float(x)) for x in line.split()] for line in file]     
    lines.append(lines_AA)
    # end read bond_AA
    # read bond_BB file
    file = open('bond_BB.xyz','r')     # 10. Co-Co    11. Co-Fe    12. Fe-Fe    
    lines_BB = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_BB)
    # end read bond_BB
    # read bond_OO file
    file = open('bond_OO.xyz','r')     # 13. O-O    
    lines_OO = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_OO)
    # end read bond_OO
    # read bond_AO file
    file = open('bond_AO.xyz','r')     # 14. La-O    15. Sr-O    
    lines_AO = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_AO)
    # end read bond_AO
    # read bond_BO file
    file = open('bond_BO.xyz','r')     # 16. Co-O    17. Fe-O
    lines_BO = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_BO)
    # end read bond_BO
    # read bond_AB file
    file = open('bond_AB.xyz','r')     # 18. La-Co    19. La-Fe    20. Sr-Co    21. Sr-Fe
    lines_AB = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_AB)
    # end read bond_AB
    # read trimer_AOA_90 file
    file = open('trimer_AOA_90.xyz','r')    # 22. La-O-La-90    23. La-O-Sr-90    24. Sr-O-Sr-90    
    lines_AOA90 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_AOA90)
    # end read trimer_AOA_90
    # read trimer_AOA_00 file
    file = open('trimer_AOA_00.xyz','r')    # 25. La-O-La-00    26. La-O-Sr-00    27. Sr-O-Sr-00    
    lines_AOA00 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_AOA00)
    # end read trimer_AOA_00 
    # read trimer_BOB file
    file = open('trimer_BOB_00.xyz','r')    # 28. Co-O-Co    29. Co-O-Fe    30. Fe-O-Fe
    lines_BOB = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_BOB)
    # end read trimer_BOB
    # read trimer_AOB file
    file = open('trimer_AOB.xyz','r')    # 31. La-O-Co    32. La-O-Fe    33. Sr-O-Co    34. Sr-O-Fe
    lines_AOB = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_AOB)
    # end read trimer_AOB
    # read trimer_OOO_60 file
    file = open('trimer_OOO_60.xyz','r')    # 35. O-O-O-60    
    lines_OOO60 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_OOO60)
    # end read trimer_OOO_60
    # read trimer_OOO_90 file
    file = open('trimer_OOO_90.xyz','r')    # 36. O-O-O-90
    lines_OOO90 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_OOO90)
    # end read trimer_OOO_90
    # read trimer_OBO_90 file
    file = open('trimer_OBO_90.xyz','r')    # 37. O-Co-O-90    38. O-Fe-O-90    
    lines_OBO90 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_OBO90)
    # end read trimer_OBO_90
    # read trimer_OBO_00 file
    file = open('trimer_OBO_00.xyz','r')   # 39. O-Co-O-00    40. O-Fe-O-00           
    lines_OBO00 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_OBO00)
    # end read trimer_OBO_00
    # read fourth_OOOO_90 file
    file = open('fourth_OOOO_90.xyz','r')            
    lines_OOOO_90 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_OOOO_90)
    # end read fourth_OOOO_90
    # read fourth_OOOO_60 file
    file = open('fourth_OOOO_60.xyz','r')              
    lines_OOOO_60 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_OOOO_60)
    # end read fourth_OOOO_60
    file = open('trimer_BBB_90.xyz','r')              
    lines_BBB_90 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_BBB_90)
    # end read trimer_BBB_90
    file = open('trimer_AAA_90.xyz','r')              
    lines_AAA_90 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_AAA_90)
    # end read lines_AAA_90
    file = open('trimer_ABA_90.xyz','r')              
    lines_ABA_90 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_ABA_90)
    # end read lines_ABA_90
    file = open('trimer_BAB_90.xyz','r')              
    lines_BAB_90 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_BAB_90)
    # end read lines_BAB_90
    file = open('trimer_BBB_180.xyz','r')              
    lines_BBB_180 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_BBB_180)
    # end read lines_BBB_180
    file = open('fourth_BOOO_60.xyz','r')              
    lines_BOOO_60 = [[int(float(x)) for x in line.split()] for line in file]
    lines.append(lines_BOOO_60)
    # end read fourth_BOOO_60
    return lines
#end def

def calc_cluster(Aset, Bset, Oset, num_atom, lines):
    array_cluster = np.zeros(27)
    array_cluster[0] = 1      # 1. V0
    # load lines
    lines_AA = lines[0]
    lines_BB = lines[1]
    lines_OO = lines[2]
    lines_AO = lines[3]
    lines_BO = lines[4]
    lines_AB = lines[5]
    lines_AOA90 = lines[6]
    lines_AOA00 = lines[7]
    lines_BOB = lines[8]
    lines_AOB = lines[9]
    lines_OOO60 = lines[10]
    lines_OOO90 = lines[11]
    lines_OBO90 = lines[12]
    lines_OBO00 = lines[13]
    lines_OOOO_90 = lines[14]
    lines_OOOO_60 = lines[15]
    lines_BBB_90 = lines[16]
    lines_AAA_90 = lines[17]
    lines_ABA_90 = lines[18]
    lines_BAB_90 = lines[19]
    lines_BBB_180 = lines[20]
    lines_BOOO_60 = lines[21]
    """
    number of fitting parameters:
    1. V0    2. La    3. Co    4. O
    5. La-La    6. Co-Co    
    7. O-O    8. La-O    9. Co-O
    10. La-Co
    11. La-O-La-90    12. La-O-La-00    
    13. Co-O-Co
    14. La-O-Co
    15. O-O-O-60    16. O-O-O-90
    17. O-Co-O-90
	18. O-Co-O-00    19. O-Fe-O-00        
    """
    array_cluster[1] = num_atom[0]  # 2. La    3. Co    4. O
    array_cluster[2] = num_atom[2]
    array_cluster[3] = num_atom[4]
    # read bond_AA file
    for i in range(len(lines_AA)):    # 5. La-La
        if ((Aset.species(lines_AA[i][0]) + Aset.species(lines_AA[i][1])) == 0):
            array_cluster[4] += 1     
    # end read bond_AA
    # read bond_BB file
    for i in range(len(lines_BB)):    # 6. Co-Co  
        if ((Bset.species(lines_BB[i][0]) + Bset.species(lines_BB[i][1])) == 0):
            array_cluster[5] += 1
    # end read bond_BB
    # read bond_OO file
    for i in range(len(lines_OO)):    # 7. O-O    
        if ((Oset.species(lines_OO[i][0]) + Oset.species(lines_OO[i][1])) == 0):
            array_cluster[6] += 1
    # end read bond_OO
    # read bond_AO file
    for i in range(len(lines_AO)):     # 8. La-O    
        if ((Aset.species(lines_AO[i][0]) == 0) and (Oset.species(lines_AO[i][1]) == 0)):
            array_cluster[7] += 1
    # end read bond_AO
    # read bond_BO file
    for i in range(len(lines_BO)):    # 9. Co-O   
        if ((Bset.species(lines_BO[i][0]) == 0) and (Oset.species(lines_BO[i][1]) == 0)):
            array_cluster[8] += 1
    # end read bond_BO
    # read bond_AB file
    for i in range(len(lines_AB)):    # 10. La-Co
        if ((Aset.species(lines_AB[i][0]) == 0) and (Bset.species(lines_AB[i][1]) == 0)):
            array_cluster[9] += 1
    # end read bond_AB
    # read trimer_AOA_90 file
    for i in range(len(lines_AOA90)):    # 11. La-O-La-90    
        if (Oset.species(lines_AOA90[i][1]) == 0):
            if ((Aset.species(lines_AOA90[i][0])+Aset.species(lines_AOA90[i][2])) == 0):
                array_cluster[10] += 1
    # end read trimer_AOA_90
    # read trimer_AOA_00 file
    for i in range(len(lines_AOA00)):    # 12. La-O-La-00   
        if (Oset.species(lines_AOA00[i][1]) == 0):
            if ((Aset.species(lines_AOA00[i][0])+Aset.species(lines_AOA00[i][2])) == 0):
                array_cluster[11] += 1
    # end read trimer_AOA_00 
    # read trimer_BOB file
    for i in range(len(lines_BOB)):    # 13. Co-O-Co
        if (Oset.species(lines_BOB[i][1]) == 0):
            if ((Bset.species(lines_BOB[i][0])+Bset.species(lines_BOB[i][2])) == 0):
                array_cluster[12] += 1
    # end read trimer_BOB
    # read trimer_AOB file
    for i in range(len(lines_AOB)):    # 14. La-O-Co
        if (Oset.species(lines_AOB[i][1]) == 0):
            if ((Aset.species(lines_AOB[i][0]) == 0) and (Bset.species(lines_AOB[i][2]) == 0)):
                array_cluster[13] += 1
    # end read trimer_AOB
    # read trimer_OOO_60 file
    for i in range(len(lines_OOO60)):    # 15. O-O-O-60    
        if ((Oset.species(lines_OOO60[i][0]) == 0) and (Oset.species(lines_OOO60[i][1]) == 0) and (Oset.species(lines_OOO60[i][2]) == 0)):
            array_cluster[14] += 1
    # end read trimer_OOO_60
    # read trimer_OOO_90 file
    for i in range(len(lines_OOO90)):    # 16. O-O-O-90
        if ((Oset.species(lines_OOO90[i][0]) == 0) and (Oset.species(lines_OOO90[i][1]) == 0) and (Oset.species(lines_OOO90[i][2]) == 0)):
            array_cluster[15] += 1
    # end read trimer_OOO_90
    # read trimer_OBO_90 file
    for i in range(len(lines_OBO90)):    # 17. O-Co-O-90    
        if ((Oset.species(lines_OBO90[i][0]) == 0) and (Oset.species(lines_OBO90[i][2]) == 0)):
            if (Bset.species(lines_OBO90[i][1]) == 0):
                array_cluster[16] += 1
    # end read trimer_OBO_90
    # read trimer_OBO_00 file
    for i in range(len(lines_OBO00)):    # 18. O-Co-O-00    19. O-Fe-O-00           
        if ((Oset.species(lines_OBO00[i][0]) == 0) and (Oset.species(lines_OBO00[i][2]) == 0)):
            if (Bset.species(lines_OBO00[i][1]) == 0):
                array_cluster[17] += 1
            else:
                array_cluster[18] += 1
    # end read trimer_OBO_00
    # read fourth_OOOO_90 file
    for i in range(len(lines_OOOO_90)):    # 20 O-O-O-O-90         
        if ((Oset.species(lines_OOOO_90[i][0]) == 0) and (Oset.species(lines_OOOO_90[i][1]) == 0) and (Oset.species(lines_OOOO_90[i][2]) == 0) and (Oset.species(lines_OOOO_90[i][3]) == 0)):
            array_cluster[19] += 1
    # end read fourth_OOOO_90
    for i in range(len(lines_AAA_90)):    # 23. A-A-A-90
        if ((Aset.species(lines_AAA_90[i][0]) == 0) and (Aset.species(lines_AAA_90[i][1]) == 0) and (Aset.species(lines_AAA_90[i][2]) == 0)):
            array_cluster[20] += 1
    # end read trimer_AAA_90
    for i in range(len(lines_ABA_90)):    # 23. A-A-A-90
        if ((Aset.species(lines_ABA_90[i][0]) == 0) and (Bset.species(lines_ABA_90[i][1]) == 0) and (Aset.species(lines_ABA_90[i][2]) == 0)):
            array_cluster[21] += 1
    # end read trimer_ABA_90
    for i in range(len(lines_BAB_90)):    # 23. A-A-A-90
        if ((Bset.species(lines_BAB_90[i][0]) == 0) and (Aset.species(lines_BAB_90[i][1]) == 0) and (Bset.species(lines_BAB_90[i][2]) == 0)):
            array_cluster[22] += 1
    # end read trimer_ABA_90
    for i in range(len(lines_BOOO_60)):    # 21 B-O-O-O-60         
        if ((Bset.species(lines_BOOO_60[i][0]) == 0) and (Oset.species(lines_BOOO_60[i][1]) == 0) and (Oset.species(lines_BOOO_60[i][2]) == 0) and (Oset.species(lines_BOOO_60[i][3]) == 0)):
            array_cluster[23] += 1
    # end read lines_BOOO_60
    # read fourth_OOOO_60 file
    for i in range(len(lines_OOOO_60)):    # 21 O-O-O-O-60         
        if ((Oset.species(lines_OOOO_60[i][0]) == 0) and (Oset.species(lines_OOOO_60[i][1]) == 0) and (Oset.species(lines_OOOO_60[i][2]) == 0) and (Oset.species(lines_OOOO_60[i][3]) == 0)):
            array_cluster[24] += 1
    # end read fourth_OOOO_60
    for i in range(len(lines_BBB_90)):    # 22. B-B-B-90
        if ((Bset.species(lines_BBB_90[i][0]) == 0) and (Bset.species(lines_BBB_90[i][1]) == 0) and (Bset.species(lines_BBB_90[i][2]) == 0)):
            array_cluster[25] += 1
    # end read trimer_BBB_90
    for i in range(len(lines_BBB_180)):    # 23. B-B-B-180
        if ((Bset.species(lines_BBB_180[i][0]) == 0) and (Bset.species(lines_BBB_180[i][1]) == 0) and (Bset.species(lines_BBB_180[i][2]) == 0)):
            array_cluster[26] += 1
    # end read trimer_BBB_90

    return array_cluster
#end def

def read_energy(d):
    with open('fe.dat') as myfile:     # write out information in fe.dat
        linesplit = list(myfile)[-1].split()
        energy1 = float(linesplit[2])
        energy2 = 0
    # Check wheter restart file is in the directory
    flag = silentchdir(d_2+"/restart")
    if (flag): 
        with open('fe.dat') as myfile:     # write out information in fe.dat
            linesplit = list(myfile)[-1].split()
            energy2 = float(linesplit[2])
    energy = min(energy1,energy2)
    return energy
#end def

if __name__ == '__main__':

    Aset, Bset, Oset = init_neighbor()
    path = "/mnt/a/u/sciteam/ma4/vasp_9.27/k_2/"
    path2 = "/mnt/a/u/sciteam/ma4/vasp_9.27/k_2/V7/La_15Sr_17Co_12Fe_20V_7"
    mat_cluster = np.zeros((296,27))
    mat_energy = np.zeros(296)
    filename = []
    lines = load_bond()
    """
    os.chdir(path2)
    num_atom = read_poscar(Aset, Bset, Oset)
    os.chdir(path)
    mat_cluster = calc_cluster(Aset, Bset, Oset, num_atom)
    print (mat_cluster)
    """
    """
    number of fitting parameters:
    1. V0    2. La    3. Co    4. O
    5. La-La    6. Co-Co    
    7. O-O    8. La-O    9. Co-O
    10. La-Co
    11. La-O-La-90    12. La-O-La-00    
    13. Co-O-Co
    14. La-O-Co
    15. O-O-O-60    16. O-O-O-90
    17. O-Co-O-90
	18. O-Co-O-00    19. O-Fe-O-00        
    """
    subpath = SubDirPath(path) # directories V0--V7
    ite = 0
    for d in subpath:
        os.chdir(d)
        subpath_2 = SubDirPath(d)  # each calculation in V0--V7
        for d_2 in subpath_2:
            os.chdir(d_2)
            #Aset, Bset, Oset = init_neighbor()
            num_atom = read_poscar(Aset, Bset, Oset)
            mat_energy[ite] = read_energy(d_2)
            filename.append(os.path.basename(d_2))
            os.chdir(path)
            time_initial = time.time()
            print (os.path.basename(d_2))
            mat_cluster[ite] = calc_cluster(Aset, Bset, Oset, num_atom, lines)
            time_end = time.time()
            print (time_end-time_initial)
            ite += 1
    outFile=open("filename","w")
    for i in range(len(filename)):
        outFile.write(str(filename[i])+"        ")
        outFile.write("\n")
    outFile.close()
    outFile=open("mat_energy","w")
    for i in range(len(mat_energy)):
        outFile.write(str(mat_energy[i])+"        ")
        outFile.write("\n")
    outFile.close()
    outFile=open("mat_cluster","w")
    for i in range(len(mat_cluster)):
        for j in range(len(mat_cluster[i])):
            outFile.write(str(mat_cluster[i][j])+"        ")
        outFile.write("\n")
    outFile.close()
# end __main__
