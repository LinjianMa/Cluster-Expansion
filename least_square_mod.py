import numpy as np
from numpy.linalg import lstsq
import math, re, random
import pylab

def mix_energy(name, energy):
    x = int(name[0])
    y = int(name[2])
    delta = int(name[4])
    
    energy_LCO = -1138.954362/32
    energy_LFO = -1208.826435/32
    energy_SCO = -927.962282/32
    energy_SFO = -1029.69172/32
    energy_O = -9.86224/2
    energy = energy - energy_LCO*x*y/32.0 - energy_LFO*x*(32-y)/32.0 - energy_SCO*(32-x)*y/32.0 - energy_SFO*(32-x)*(32-y)/32.0 + energy_O*delta
    """
    energy_La2O3 = -671.12741355/16
    energy_SrO = -385.71932805/32
    energy_Fe2O3 = -541.511305319/16
    energy_Co3O4 = -336.275977085/8
    delta_La2O3 = x/2.0
    delta_SrO = 32-x
    delta_Fe2O3 = (32-y)/2.0
    delta_Co3O4 = y/3.0
    energy_O2 = -9.86224
    delta_O2 = 1.0/4*x+(2.0/3-3.0/4)*y-8+0.5*delta
    energy = energy - energy_La2O3*delta_La2O3 - energy_SrO*delta_SrO - energy_Fe2O3*delta_Fe2O3 - energy_Co3O4*delta_Co3O4 + energy_O2*delta_O2
    """
    return energy
#end def

def CEV_calc(flag):
    data_num = 296
    cluster_num = 27
    name_delete = []
    file = open('filename','r')
    lines_name = [line.split() for line in file]
    for i in range(len(lines_name)):
        name = re.findall(r"\d+", lines_name[i][0])
    num_delete = len(name_delete)
    print (num_delete)
    # initialize the cluster coefficient matrix and energy vector
    mat_cluster = np.zeros((data_num-num_delete,cluster_num))
    mat_energy = np.zeros(data_num-num_delete)
    mat_energy_orig = np.zeros(data_num-num_delete)
    # end initialization
    # load the energy vector
    file = open('mat_energy','r')
    lines_energy = [line.split() for line in file]
    delete = 0
    for i in range(data_num):
        if i not in name_delete:
            mat_energy[i-delete] = float(lines_energy[i][0])
            mat_energy_orig[i-delete] = mat_energy[i-delete]
            mat_energy[i-delete] = mix_energy(re.findall(r"\d+", lines_name[i][0]), mat_energy[i-delete])
        else: 
            delete += 1

    # load the cluster coefficient matrix
    file = open('mat_cluster','r')  
    lines_cluster = [line.split() for line in file]
    delete = 0
    for i in range(data_num):
        if i not in name_delete:
            for j in range(len(mat_cluster[i-delete])):
                mat_cluster[i-delete][j] = float(lines_cluster[i][j])
        else: 
            delete += 1
    """
    # random sample the overall space
    matorig = np.concatenate((mat_cluster, mat_energy.reshape((len(mat_energy),1))), axis=1)
    mat = matorig.copy()
    ran_num = np.random.choice(len(matorig), len(matorig),  replace=False)
    for i in range(len(matorig)):
        mat[i] = matorig[ran_num[i]]
        mat_cluster[i] = mat[i][:len(mat[0])-1]
        mat_energy[i] = mat[i][len(mat[0])-1]
    """

    # least square calculations
    V = lstsq(mat_cluster[:,:24], mat_energy[:])[0]
    print (V)
    error = (np.dot(mat_cluster[:,:24],V)-mat_energy[:])/160*1000
    np.set_printoptions(precision=8)
    np.set_printoptions(suppress=True)
    # calculating RMSE
    RMSE = 0
    for i in range(len(error)):
        RMSE += error[i]*error[i]
    RMSE = math.sqrt(RMSE/len(error))
    print (RMSE)
    if (flag == True):
        return V, RMSE
    else: 
        return V, RMSE, mat_cluster, mat_energy, mat_energy_orig, lines_name, error, name_delete 
#end def

if __name__ == '__main__':
    vacancy_0 = []
    vacancy_1 = []
    vacancy_2 = []
    vacancy_3 = []
    vacancy_4 = []
    vacancy_5 = []
    vacancy_6 = []
    vacancy_7 = []
    mat_cluster_vac_0 = []
    mat_cluster_vac_1 = []
    mat_cluster_vac_2 = []
    mat_cluster_vac_3 = []
    mat_cluster_vac_4 = []
    mat_cluster_vac_5 = []
    mat_cluster_vac_6 = []
    mat_cluster_vac_7 = []
    mat_energy_vac_0 = []
    mat_energy_vac_1 = []
    mat_energy_vac_2 = []
    mat_energy_vac_3 = []
    mat_energy_vac_4 = []
    mat_energy_vac_5 = []
    mat_energy_vac_6 = []
    mat_energy_vac_7 = []

    V, RMSE, mat_cluster, mat_energy, mat_energy_orig, lines_name, error, name_delete = CEV_calc(False)

    delete = 0
    for i in range(296):
        if i not in name_delete:
            # chech the number of vacancy
            name = re.findall(r"\d+", lines_name[i][0])
            if (int(name[4]) == 0):
                vacancy_0.append(i-delete)
            elif (int(name[4]) == 1):
                vacancy_1.append(i-delete)
            elif (int(name[4]) == 2):
                vacancy_2.append(i-delete)
            elif (int(name[4]) == 3):
                vacancy_3.append(i-delete)
            elif (int(name[4]) == 4):
                vacancy_4.append(i-delete)
            elif (int(name[4]) == 5):
                vacancy_5.append(i-delete)
            elif (int(name[4]) == 6):
                vacancy_6.append(i-delete)
            elif (int(name[4]) == 7):
                vacancy_7.append(i-delete)
            # end check the num of vacancy
        else: 
            delete += 1

    outFile=open("error.mod","w")
    for i in range(len(error)):
        outFile.write(str(i+1)+"        "+str(error[i])+"        ")
        outFile.write("\n")
    outFile.close()

    outFile=open("energy.out","w")
    for i in range(len(mat_energy_orig)):
        outFile.write(str(i+1)+"        "+str(mat_energy_orig[i])+"        ")
        outFile.write("\n")
    outFile.close()

   
    #Plot the relation between mat_cluster and mat_energy
    pylab.plot(np.dot(mat_cluster,V)/32,mat_energy/32,'o',linewidth=3.0,label='Experiment Curve')
    pylab.plot(mat_energy/32,mat_energy/32,'-',linewidth=3.0,label='Experiment Curve')
    #pylab.legend(bbox_to_anchor=(0, 0.7), loc=2, borderaxespad=0.)
    #pylab.title('Overpotential-Current density curve for SOEC/SOFC')
    pylab.ylabel('DFT formation energy (eV/unit cell)')
    pylab.xlabel('cluster expansion formation energy (eV/unit cell)')
    #pylab.axis([-1.5,1.5,-1,1])
    pylab.grid(True)
    pylab.show()

    # print the distribution of different vacancies
    for i in range(len(mat_energy)):
        if i in vacancy_0 : 
            mat_cluster_vac_0.append(mat_cluster[i,:])
            mat_energy_vac_0.append(mat_energy[i])
        if i in vacancy_1 : 
            mat_cluster_vac_1.append(mat_cluster[i,:])
            mat_energy_vac_1.append(mat_energy[i])
        if i in vacancy_2 : 
            mat_cluster_vac_2.append(mat_cluster[i,:])
            mat_energy_vac_2.append(mat_energy[i])
        if i in vacancy_3 : 
            mat_cluster_vac_3.append(mat_cluster[i,:])
            mat_energy_vac_3.append(mat_energy[i])
        if i in vacancy_4 : 
            mat_cluster_vac_4.append(mat_cluster[i,:])
            mat_energy_vac_4.append(mat_energy[i])
        if i in vacancy_5 : 
            mat_cluster_vac_5.append(mat_cluster[i,:])
            mat_energy_vac_5.append(mat_energy[i])
        if i in vacancy_6 : 
            mat_cluster_vac_6.append(mat_cluster[i,:])
            mat_energy_vac_6.append(mat_energy[i])
        if i in vacancy_7 : 
            mat_cluster_vac_7.append(mat_cluster[i,:])
            mat_energy_vac_7.append(mat_energy[i])
    mat_cluster_vac_0 = np.asarray(mat_cluster_vac_0)
    mat_energy_vac_0 = np.asarray(mat_energy_vac_0)
    mat_cluster_vac_1 = np.asarray(mat_cluster_vac_1)
    mat_energy_vac_1 = np.asarray(mat_energy_vac_1)
    mat_cluster_vac_2 = np.asarray(mat_cluster_vac_2)
    mat_energy_vac_2 = np.asarray(mat_energy_vac_2)
    mat_cluster_vac_3 = np.asarray(mat_cluster_vac_3)
    mat_energy_vac_3 = np.asarray(mat_energy_vac_3)
    mat_cluster_vac_4 = np.asarray(mat_cluster_vac_4)
    mat_energy_vac_4 = np.asarray(mat_energy_vac_4)
    mat_cluster_vac_5 = np.asarray(mat_cluster_vac_5)
    mat_energy_vac_5 = np.asarray(mat_energy_vac_5)
    mat_cluster_vac_6 = np.asarray(mat_cluster_vac_6)
    mat_energy_vac_6 = np.asarray(mat_energy_vac_6)
    mat_cluster_vac_7 = np.asarray(mat_cluster_vac_7)
    mat_energy_vac_7 = np.asarray(mat_energy_vac_7)

    #Plot the relation between mat_cluster and mat_energy
    pylab.plot(np.dot(mat_cluster_vac_0,V)/32,mat_energy_vac_0/32,'bo',linewidth=3.0,label='0 vacancy')
    pylab.plot(np.dot(mat_cluster_vac_1,V)/32,mat_energy_vac_1/32,'ro',linewidth=3.0,label='1 vacancy')
    pylab.plot(np.dot(mat_cluster_vac_2,V)/32,mat_energy_vac_2/32,'go',linewidth=3.0,label='2 vacancy')
    pylab.plot(np.dot(mat_cluster_vac_3,V)/32,mat_energy_vac_3/32,'co',linewidth=3.0,label='3 vacancy')
    pylab.plot(np.dot(mat_cluster_vac_4,V)/32,mat_energy_vac_4/32,'b^',linewidth=3.0,label='4 vacancy')
    pylab.plot(np.dot(mat_cluster_vac_5,V)/32,mat_energy_vac_5/32,'r^',linewidth=3.0,label='5 vacancy')
    pylab.plot(np.dot(mat_cluster_vac_6,V)/32,mat_energy_vac_6/32,'g^',linewidth=3.0,label='6 vacancy')
    pylab.plot(np.dot(mat_cluster_vac_7,V)/32,mat_energy_vac_7/32,'c^',linewidth=3.0,label='7 vacancy')
    pylab.plot(mat_energy/32,mat_energy/32,'-',linewidth=3.0)
    pylab.legend(bbox_to_anchor=(0, 0.9), loc=2, borderaxespad=0.)
    #pylab.title('Overpotential-Current density curve for SOEC/SOFC')
    pylab.ylabel('DFT formation energy (eV/unit cell)')
    pylab.xlabel('cluster expansion formation energy (eV/unit cell)')
    #pylab.axis([-1.5,1.5,-1,1])
    pylab.grid(True)
    pylab.show()
    
    """
    #Plot the relation between mat_cluster and mat_energy
    pylab.plot(np.dot(mat_cluster_vac_0,V)/32,mat_energy_vac_0/32,'bo',linewidth=3.0,label='0 vacancy')
    pylab.plot(mat_energy/32,mat_energy/32,'-',linewidth=3.0)
    pylab.legend(bbox_to_anchor=(0, 0.9), loc=2, borderaxespad=0.)
    #pylab.title('Overpotential-Current density curve for SOEC/SOFC')
    pylab.ylabel('DFT formation energy (eV/unit cell)')
    pylab.xlabel('cluster expansion formation energy (eV/unit cell)')
    #pylab.axis([-1.5,1.5,-1,1])
    pylab.grid(True)
    pylab.show()

    #Plot the relation between mat_cluster and mat_energy
    pylab.plot(np.dot(mat_cluster_vac_1,V)/32,mat_energy_vac_1/32,'ro',linewidth=3.0,label='1 vacancy')
    pylab.plot(mat_energy/32,mat_energy/32,'-',linewidth=3.0)
    pylab.legend(bbox_to_anchor=(0, 0.9), loc=2, borderaxespad=0.)
    #pylab.title('Overpotential-Current density curve for SOEC/SOFC')
    pylab.ylabel('DFT formation energy (eV/unit cell)')
    pylab.xlabel('cluster expansion formation energy (eV/unit cell)')
    #pylab.axis([-1.5,1.5,-1,1])
    pylab.grid(True)
    pylab.show()

    #Plot the relation between mat_cluster and mat_energy
    pylab.plot(np.dot(mat_cluster_vac_2,V)/32,mat_energy_vac_2/32,'go',linewidth=3.0,label='2 vacancy')
    pylab.plot(mat_energy/32,mat_energy/32,'-',linewidth=3.0)
    pylab.legend(bbox_to_anchor=(0, 0.9), loc=2, borderaxespad=0.)
    #pylab.title('Overpotential-Current density curve for SOEC/SOFC')
    pylab.ylabel('DFT formation energy (eV/unit cell)')
    pylab.xlabel('cluster expansion formation energy (eV/unit cell)')
    #pylab.axis([-1.5,1.5,-1,1])
    pylab.grid(True)
    pylab.show()

    #Plot the relation between mat_cluster and mat_energy
    pylab.plot(np.dot(mat_cluster_vac_3,V)/32,mat_energy_vac_3/32,'co',linewidth=3.0,label='3 vacancy')
    pylab.plot(mat_energy/32,mat_energy/32,'-',linewidth=3.0)
    pylab.legend(bbox_to_anchor=(0, 0.9), loc=2, borderaxespad=0.)
    #pylab.title('Overpotential-Current density curve for SOEC/SOFC')
    pylab.ylabel('DFT formation energy (eV/unit cell)')
    pylab.xlabel('cluster expansion formation energy (eV/unit cell)')
    #pylab.axis([-1.5,1.5,-1,1])
    pylab.grid(True)
    pylab.show()

    #Plot the relation between mat_cluster and mat_energy
    pylab.plot(np.dot(mat_cluster_vac_4,V)/32,mat_energy_vac_4/32,'b^',linewidth=3.0,label='4 vacancy')
    pylab.plot(mat_energy/32,mat_energy/32,'-',linewidth=3.0)
    pylab.legend(bbox_to_anchor=(0, 0.9), loc=2, borderaxespad=0.)
    #pylab.title('Overpotential-Current density curve for SOEC/SOFC')
    pylab.ylabel('DFT formation energy (eV/unit cell)')
    pylab.xlabel('cluster expansion formation energy (eV/unit cell)')
    #pylab.axis([-1.5,1.5,-1,1])
    pylab.grid(True)
    pylab.show()

    #Plot the relation between mat_cluster and mat_energy
    pylab.plot(np.dot(mat_cluster_vac_5,V)/32,mat_energy_vac_5/32,'r^',linewidth=3.0,label='5 vacancy')
    pylab.plot(mat_energy/32,mat_energy/32,'-',linewidth=3.0)
    pylab.legend(bbox_to_anchor=(0, 0.9), loc=2, borderaxespad=0.)
    #pylab.title('Overpotential-Current density curve for SOEC/SOFC')
    pylab.ylabel('DFT formation energy (eV/unit cell)')
    pylab.xlabel('cluster expansion formation energy (eV/unit cell)')
    #pylab.axis([-1.5,1.5,-1,1])
    pylab.grid(True)
    pylab.show()

    #Plot the relation between mat_cluster and mat_energy
    pylab.plot(np.dot(mat_cluster_vac_6,V)/32,mat_energy_vac_6/32,'g^',linewidth=3.0,label='6 vacancy')
    pylab.plot(mat_energy/32,mat_energy/32,'-',linewidth=3.0)
    pylab.legend(bbox_to_anchor=(0, 0.9), loc=2, borderaxespad=0.)
    #pylab.title('Overpotential-Current density curve for SOEC/SOFC')
    pylab.ylabel('DFT formation energy (eV/unit cell)')
    pylab.xlabel('cluster expansion formation energy (eV/unit cell)')
    #pylab.axis([-1.5,1.5,-1,1])
    pylab.grid(True)
    pylab.show()

    #Plot the relation between mat_cluster and mat_energy
    pylab.plot(np.dot(mat_cluster_vac_7,V)/32,mat_energy_vac_7/32,'c^',linewidth=3.0,label='7 vacancy')
    pylab.plot(mat_energy/32,mat_energy/32,'-',linewidth=3.0)
    pylab.legend(bbox_to_anchor=(0, 0.9), loc=2, borderaxespad=0.)
    #pylab.title('Overpotential-Current density curve for SOEC/SOFC')
    pylab.ylabel('DFT formation energy (eV/unit cell)')
    pylab.xlabel('cluster expansion formation energy (eV/unit cell)')
    #pylab.axis([-1.5,1.5,-1,1])
    pylab.grid(True)
    pylab.show()
    """
# end __main__