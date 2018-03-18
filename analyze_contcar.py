import numpy as np
import math, shutil
import os, sys, random, errno
from numpy.linalg import norm
#import matplotlib

def SubDirPath (d):
    return filter(os.path.isdir, [os.path.join(d,f) for f in os.listdir(d)])

def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred
    return

def silentchdir(filename):
    output = True     # Flag of chdir
    try:
        os.chdir(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        output = False
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred
    return output

def delete_chargecar (path):
    subpath = SubDirPath(path) # directories V0--V7
    for d in subpath:
        os.chdir(d)
        subpath_2 = SubDirPath(d)  # each calculation in V0--V7
        for d_2 in subpath_2:
            os.chdir(d_2)
            silentremove("CHGCAR")
            silentremove("CHG")
            subpath_3 = SubDirPath(d_2)  # restart/restart1/...
            for d_3 in subpath_3:
                os.chdir(d_3)
                silentremove("CHGCAR")
                silentremove("CHG")
    return

def vefpl (path):              # run the ref.pl
    subpath = SubDirPath(path) # directories V0--V7
    for d in subpath:
        os.chdir(d)
        subpath_2 = SubDirPath(d)  # each calculation in V0--V7
        for d_2 in subpath_2:
            os.chdir(d_2)
            os.system("../.././vef.pl > /dev/null")
            subpath_3 = SubDirPath(d_2)  # restart/restart1/...
            for d_3 in subpath_3:
                os.chdir(d_3)
                os.system("../../.././vef.pl > /dev/null")
    return

def energy_output (path):
    filename = []
    energy = []
    error = []
    energy_restart = []
    error_restart = []    
    subpath = SubDirPath(path) # directories V0--V7
    for d in subpath:
        os.chdir(d)
        subpath_2 = SubDirPath(d)  # each calculation in V0--V7
        for d_2 in subpath_2:
            os.chdir(d_2)
            filename.append(os.path.basename(d_2))
            with open('fe.dat') as myfile:     # write out information in fe.dat
                linesplit = list(myfile)[-1].split()
                error.append(float(linesplit[1]))
                energy.append(float(linesplit[2]))
                error_restart.append(0.0)
                energy_restart.append(0.0)
            # Check wheter restart file is in the directory
            flag = silentchdir(d_2+"/restart")
            if (flag): 
                with open('fe.dat') as myfile:     # write out information in fe.dat
                    linesplit = list(myfile)[-1].split()
                    error_restart[-1] = float(linesplit[1])
                    energy_restart[-1] = float(linesplit[2])
    os.chdir(path)
    # Write to file [energy_out.txt]
    outFile=open("energy_out","w")
    outFile.write("Filename"+"        "+"error"+"        "+"energy"+"        "+"error_restart"+"        "+"energy_restart"+"\n")
    for i in range(len(filename)):
        outFile.write(str(filename[i])+"        "+str(error[i])+"        "+str(energy[i])+"        "+str(error_restart[i])+"        "+str(energy_restart[i])+"\n")
    outFile.close()
    # End writing [energy_out.txt]
    return

def distance_calc ():
    dist = 0.0         # distance between poscar and contcar
    flag = False       # wheter distance exceeds the limiting length
    length = [10.86964544, 10.86964544, 15.372] # cell length
    delta = [0, 0.25, 0]  # distance between two neighboring oxygen
    limit = 0.5*norm(np.dot(length,delta))    # the limit length when considering wheter oxygen has moved to another position
    file = open('POSCAR','r')
    lines_poscar = [line.split() for line in file]
    file = open('CONTCAR','r')
    lines_contcar = [line.split() for line in file]
    num_atom = 0       # number of atoms
    for i in range(len(lines_poscar[6])):
        num_atom = num_atom + int(lines_poscar[6][i])
    for i in range(9,9+num_atom):
        coord_pos = np.array([float(lines_poscar[i][0]),float(lines_poscar[i][1]),float(lines_poscar[i][2])])
        coord_cont = np.array([float(lines_contcar[i][0]),float(lines_contcar[i][1]),float(lines_contcar[i][2])])
        for j in range(3):      # considering the periodic boundary condition
            coord_cont[j] = coord_cont[j] - int(round((coord_cont[j]-coord_pos[j])))
	    vector = [(coord_cont[0]-coord_pos[0])*length[0],(coord_cont[1]-coord_pos[1])*length[1],(coord_cont[2]-coord_pos[2])*length[2]]
        dist = max(dist, norm(vector))
    if (dist > limit):
        flag = True
    print (flag)
    return (dist,flag)

def structure_analysis (path):    # Analyze whether the poscar and contcar files are similar
    filename = []
    distance = []
    flagdist = []
    distance_restart = []
    flagdist_restart = []
    subpath = SubDirPath(path) # directories V0--V7
    for d in subpath:
        os.chdir(d)
        subpath_2 = SubDirPath(d)  # each calculation in V0--V7
        for d_2 in subpath_2:
            os.chdir(d_2)
            filename.append(os.path.basename(d_2))
            dist,f = distance_calc()
            distance.append(dist)
            flagdist.append(f)
            distance_restart.append(0.0)
            flagdist_restart.append(False)
            # Check wheter restart file is in the directory
            flag = silentchdir(d_2+"/restart")
            if (flag): 
                dist,f = distance_calc()
                distance_restart[-1] = dist
                flagdist_restart[-1] = f
    os.chdir(path)
    # Write to file [distance_out.txt]
    outFile=open("distance_out","w")
    outFile.write("Filename"+"        "+"distance"+"        "+"flag"+"        "+"distance_restart"+"        "+"flag_restart"+"\n")
    for i in range(len(filename)):
        outFile.write(str(filename[i])+"        "+str(distance[i])+"        "+str(flagdist[i])+"        "+str(distance_restart[i])+"        "+str(flagdist_restart[i])+"\n")
    outFile.close()
    # End writing [distance_out.txt]
    return

if __name__ == '__main__':
    path = "/mnt/a/u/sciteam/ma4/vasp_9.27/k_2/"
    delete_chargecar(path)
    vefpl(path)
    energy_output(path)
    structure_analysis(path)
# end __main__
