import numpy as np
import math, shutil
import os, sys, random
#import matplotlib

def write_coords(lattice_x, lattice_y, lattice_z, num_La, num_Sr, num_Co, num_Fe, num_O):
    La_Sr = []
    Co_Fe = []
    O = []
    #La, Sr, 32 atoms in all
    for k in range (4):
        La_Sr.append(str(0.00)+"	"+str(0.00)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")    
        La_Sr.append(str(0.00)+"	"+str(0.50)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        La_Sr.append(str(0.50)+"	"+str(0.00)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        La_Sr.append(str(0.50)+"	"+str(0.50)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        La_Sr.append(str(0.25)+"	"+str(0.25)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        La_Sr.append(str(0.25)+"	"+str(0.75)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        La_Sr.append(str(0.75)+"	"+str(0.25)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        La_Sr.append(str(0.75)+"	"+str(0.75)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
    #Co, Fe, 32 atoms in all
    for k in range (4):
        Co_Fe.append(str(0.00)+"	"+str(0.25)+"	"+str(float(k/4.0+0.125))+"	"+"T	T	T"+"\n")    
        Co_Fe.append(str(0.00)+"	"+str(0.75)+"	"+str(float(k/4.0+0.125))+"	"+"T	T	T"+"\n")   
        Co_Fe.append(str(0.25)+"	"+str(0.00)+"	"+str(float(k/4.0+0.125))+"	"+"T	T	T"+"\n")   
        Co_Fe.append(str(0.25)+"	"+str(0.50)+"	"+str(float(k/4.0+0.125))+"	"+"T	T	T"+"\n")   
        Co_Fe.append(str(0.50)+"	"+str(0.25)+"	"+str(float(k/4.0+0.125))+"	"+"T	T	T"+"\n")   
        Co_Fe.append(str(0.50)+"	"+str(0.75)+"	"+str(float(k/4.0+0.125))+"	"+"T	T	T"+"\n")   
        Co_Fe.append(str(0.75)+"	"+str(0.00)+"	"+str(float(k/4.0+0.125))+"	"+"T	T	T"+"\n")   
        Co_Fe.append(str(0.75)+"	"+str(0.50)+"	"+str(float(k/4.0+0.125))+"	"+"T	T	T"+"\n")   
    #O, 96 atoms in all
    for k in range (4):
        O.append(str(0.00)+"	"+str(0.25)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")    
        O.append(str(0.00)+"	"+str(0.75)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        O.append(str(0.25)+"	"+str(0.00)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        O.append(str(0.25)+"	"+str(0.50)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        O.append(str(0.50)+"	"+str(0.25)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        O.append(str(0.50)+"	"+str(0.75)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        O.append(str(0.75)+"	"+str(0.00)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")   
        O.append(str(0.75)+"	"+str(0.50)+"	"+str(float(k/4.0))+"	"+"T	T	T"+"\n")
    for k in range (4):
        for j in range (4):
            for i in range (4):
                O.append(str(float(i/4.0+0.125))+"	"+str(float(j/4.0+0.125))+"	"+str(float(k/4.0+0.125))+"	"+"T	T	T"+"\n")    
    # Write to file [POSCAR]
    outFile=open("POSCAR","w")
    #Writing headfile
    outFile.write("POSCAR_LSCF   "+"\n")
    outFile.write(" 1."+"\n")
    outFile.write("    "+str(lattice_x)+"    "+str(0)+"    "+str(0)+"    "+"\n")
    outFile.write("    "+str(0)+"    "+str(lattice_y)+"    "+str(0)+"    "+"\n")
    outFile.write("    "+str(0)+"    "+str(0)+"    "+str(lattice_z)+"    "+"\n")
    if (num_La==0 and num_Co==0):
        outFile.write("   "+"Sr"+"    "+"Fe"+"    "+"O"+"\n")
        outFile.write("    "+str(int(num_Sr))+"    "+str(int(num_Fe))+"    "+str(int(num_O))+"\n")
    elif (num_La==0 and num_Fe==0):
        outFile.write("   "+"Sr"+"    "+"Co"+"    "+"O"+"\n")
        outFile.write("    "+str(int(num_Sr))+"    "+str(int(num_Co))+"    "+str(int(num_O))+"\n")
    elif (num_Sr==0 and num_Co==0):
        outFile.write("   "+"La"+"    "+"Fe"+"    "+"O"+"\n")
        outFile.write("    "+str(int(num_La))+"    "+str(int(num_Fe))+"    "+str(int(num_O))+"\n")
    elif (num_Sr==0 and num_Fe==0):
        outFile.write("   "+"La"+"    "+"Co"+"    "+"O"+"\n")
        outFile.write("    "+str(int(num_La))+"    "+str(int(num_Co))+"    "+str(int(num_O))+"\n")
    elif (num_La==0):
        outFile.write("   "+"Sr"+"    "+"Co"+"    "+"Fe"+"    "+"O"+"\n")
        outFile.write("    "+str(int(num_Sr))+"    "+str(int(num_Co))+"    "+str(int(num_Fe))+"    "+str(int(num_O))+"\n")
    elif (num_Sr==0):
        outFile.write("   "+"La"+"    "+"Co"+"    "+"Fe"+"    "+"O"+"\n")
        outFile.write("    "+str(int(num_La))+"    "+str(int(num_Co))+"    "+str(int(num_Fe))+"    "+str(int(num_O))+"\n")
    elif (num_Co==0):
        outFile.write("   "+"La"+"    "+"Sr"+"    "+"Fe"+"    "+"O"+"\n")
        outFile.write("    "+str(int(num_La))+"    "+str(int(num_Sr))+"    "+str(int(num_Fe))+"    "+str(int(num_O))+"\n")
    elif (num_Fe==0):
        outFile.write("   "+"La"+"    "+"Sr"+"    "+"Co"+"    "+"O"+"\n")
        outFile.write("    "+str(int(num_La))+"    "+str(int(num_Sr))+"    "+str(int(num_Co))+"    "+str(int(num_O))+"\n")
    else:
        outFile.write("   "+"La"+"    "+"Sr"+"    "+"Co"+"    "+"Fe"+"    "+"O"+"\n")
        outFile.write("    "+str(int(num_La))+"    "+str(int(num_Sr))+"    "+str(int(num_Co))+"    "+str(int(num_Fe))+"    "+str(int(num_O))+"\n")
    outFile.write("Selective dynamics"+"\n")
    outFile.write("Direct"+"\n")
    #End headfile
    #Writing coordinates
    La_Sr = random.sample(La_Sr,32)
    Co_Fe = random.sample(Co_Fe,32)
    O = random.sample(O,96)
    for i in range(32):
        outFile.write(La_Sr[i])
    for i in range(32):
        outFile.write(Co_Fe[i])
    for i in range(96):
        outFile.write(O[i])
    outFile.close()
    # End writing [POSCAR]
    # Write to file [bwjobs.bat]
    outFile=open("bwjobs.bat","w")
    #Writing headfile
    outFile.write("#PBS -A baca"+"\n")
    outFile.write("#PBS -q low"+"\n")
    outFile.write("#PBS -l nodes=6:ppn=32:xe"+"\n")
    outFile.write("#PBS -l walltime=24:00:00"+"\n")
    outFile.write("#PBS -N "+"La_"+str(int(num_La))+"Sr_"+str(int(num_Sr))+"Co_"+str(int(num_Co))+"Fe_"+str(int(num_Fe))+"V_"+str(int(96-num_O))+"\n")
    outFile.write("#PBS -e $PBS_JOBID.err"+"\n")
    outFile.write("#PBS -o $PBS_JOBID.out"+"\n")
    outFile.write("#PBS -m bea"+"\n")
    outFile.write("#PBS -M lma16@illinois.edu"+"\n\n")
    outFile.write(". /opt/modules/default/init/bash"+"\n\n")
    outFile.write("cd $PBS_O_WORKDIR"+"\n\n")
    outFile.write("aprun -n 192 /u/sciteam/ma4/bin/vasp5.3.5_TS_neb/vasp_bw > CoT.out"+"\n\n")
    outFile.close()
    # End writing [bwjobs.bat]
    return

# start __main__
#construction of initial LSCF structure. Output is POSCAR
lattice_x = 10.86964544
lattice_y = 10.86964544
lattice_z = 7.686*2
num_La = 0
num_Sr = 32
num_SrLa = 32
num_Co = 0
num_Fe = 32
num_CoFe = 32
num_O = 96
num_V = 0

path = "/mnt/a/u/sciteam/ma4/vasp_9.27/k_2/neb/"

#When no vacancy is in the system
for i in (5,13,21): # Iteration for La and Sr
    for j in (25,17,9,5): # Iteration for Co and Fe
        for k in (0,1,2,3):
            num_La = i
            num_Sr = num_SrLa - num_La
            num_Fe = j
            num_Co = num_CoFe - num_Fe
            num_V = k
            num_O = 96 - num_V
            path2 = "La_"+str(int(num_La))+"Sr_"+str(int(num_Sr))+"Co_"+str(int(num_Co))+"Fe_"+str(int(num_Fe))+"V_"+str(int(96-num_O))+"_neb"
            os.mkdir(path+path2);
            os.chdir(path+path2);
            write_coords(lattice_x, lattice_y, lattice_z, num_La, num_Sr, num_Co, num_Fe, num_O);
            # Write to file [KPOINTS]
            shutil.copy(path+'KPOINTS',path+path2)
            # End writing [KPOINTS]
            if (num_La==0 and num_Co==0):
                shutil.copy(path+'INCAR_SF',path+path2)
                os.rename('INCAR_SF','INCAR')
                shutil.copy(path+'POTCAR_SF',path+path2)
                os.rename('POTCAR_SF','POTCAR')
            elif (num_La==0 and num_Fe==0):
                shutil.copy(path+'INCAR_SC',path+path2)
                os.rename('INCAR_SC','INCAR')
                shutil.copy(path+'POTCAR_SC',path+path2)
                os.rename('POTCAR_SC','POTCAR')
            elif (num_Sr==0 and num_Co==0):
                shutil.copy(path+'INCAR_LF',path+path2)
                os.rename('INCAR_LF','INCAR')
                shutil.copy(path+'POTCAR_LF',path+path2)
                os.rename('POTCAR_LF','POTCAR')
            elif (num_Sr==0 and num_Fe==0):
                shutil.copy(path+'INCAR_LC',path+path2)
                os.rename('INCAR_LC','INCAR')
                shutil.copy(path+'POTCAR_LC',path+path2)
                os.rename('POTCAR_LC','POTCAR')
            elif (num_La==0):
                shutil.copy(path+'INCAR_SCF',path+path2)
                os.rename('INCAR_SCF','INCAR')
                shutil.copy(path+'POTCAR_SCF',path+path2)
                os.rename('POTCAR_SCF','POTCAR')
            elif (num_Sr==0):
                shutil.copy(path+'INCAR_LCF',path+path2)
                os.rename('INCAR_LCF','INCAR')
                shutil.copy(path+'POTCAR_LCF',path+path2)
                os.rename('POTCAR_LCF','POTCAR')
            elif (num_Co==0):
                shutil.copy(path+'INCAR_LSF',path+path2)
                os.rename('INCAR_LSF','INCAR')
                shutil.copy(path+'POTCAR_LSF',path+path2)
                os.rename('POTCAR_LSF','POTCAR')
            elif (num_Fe==0):
                shutil.copy(path+'INCAR_LSC',path+path2)
                os.rename('INCAR_LSC','INCAR')
                shutil.copy(path+'POTCAR_LSC',path+path2)
                os.rename('POTCAR_LSC','POTCAR')
            else:
                shutil.copy(path+'INCAR_LSCF',path+path2)
                os.rename('INCAR_LSCF','INCAR')
                shutil.copy(path+'POTCAR_LSCF',path+path2)
                os.rename('POTCAR_LSCF','POTCAR')
            # Submit the job
            #os.system('qsub bwjobs.bat')

# end __main__
