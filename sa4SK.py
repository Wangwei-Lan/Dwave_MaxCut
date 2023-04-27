def formulate_SK_qubo_sparse(filename):
    
    #! return -Q , where the sign is introduced to convert maximum problem to minimum problem
    #! The returned qubo matrix is symmetric 
    txt_file = open(filename)
    #num_nodes,num_edges = txt_file.readline().split()
    num_nodes = 100 
    #! initial qubo sparse dictionary 
    qubo_matrix = np.zeros((num_nodes,num_nodes)) 
    
    for line in txt_file:
        line_split = line.split(' ')
        node1,node2,weight = int(line_split[0])-1,int(line_split[1])-1,float(line_split[2])
        qubo_matrix[(node1,node2)] = -weight
        qubo_matrix[(node2,node1)] = -weight
    return qubo_matrix 


def formulate_SK_ising_sparse(filename):
    
    #! return h,J 
    #! where h is the linear biases of the Ising problems
    #! J does not have on-site terms for max cut problem
    #! 
    txt_file = open(filename)
    num_nodes = 100
    #! initial h
    h = np.zeros(int(num_nodes))
    J = np.zeros((num_nodes,num_nodes))
     
    for line in txt_file:
        line_split = line.split(' ')
        node1,node2,weight = int(line_split[0])-1,int(line_split[1])-1,float(line_split[2])
        J[(node1,node2)] = -weight
        J[(node2,node1)] = -weight

    return h,J



def get_cut_number(spin_vector,edges):
    
    cut_numer = 0
    for pair in edges.keys():
        if spin_vector[pair[0]] != spin_vector[pair[1]]:
            cut_numer += edges[pair]
    return cut_numer

def get_gs_energy(spin_vector,coupling):
    pass





import numpy as np
import networkx as nx 
import os
import random 
import dimod 
import re
import time 

# construct instance graph from txt file 
from dwave.samplers import SimulatedAnnealingSampler
sa_solver = SimulatedAnnealingSampler()

from dwave.samplers import TabuSampler
tabu_solver = TabuSampler()


from dwave.samplers import SteepestDescentSolver
sd_solver = SteepestDescentSolver()



# file contains edge information
file_path = os.getcwd()
gs_file_names = os.listdir(file_path+"/SK_gs/")
num_nodes = 100


for gs_file_name in gs_file_names:

    
    if gs_file_name == 'README.md':
        continue

    
    
    seed = re.split(r'\.|\_', gs_file_name)[1] 
    filename = file_path+"/SK_N100/100_SK_"+seed+".txt"
    gs_file  = file_path+"/SK_gs/"+gs_file_name

    spin_up_idx = np.array(np.loadtxt(gs_file),dtype=int)
    gs_vector = np.ones(100)*(-1)
    gs_vector[spin_up_idx - 1] = 1.0


    # (1) solve qubo model 
    #instance_dict = formulate_SK_qubo_sparse(filename)
    h,J = formulate_SK_ising_sparse(filename)
    exact_gs_energy = gs_vector@J@gs_vector
    numruns = 1000

    start_time = time.time()
    #sampleset = solver.sample_ising(np.zeros(100),J,num_reads=numruns,beta_schedule_type='geometric',num_sweeps=1000,beta_range=[0.0001,100000.0],num_sweeps_per_beta=1,)
    #sampleset = tabu_solver.sample_ising(np.zeros(100),J,num_reads=numruns,num_restarts=10,tenure=10)
    sampleset = sd_solver.sample_ising(np.zeros(100),J,num_reads=numruns)

    end_time = time.time()
    
    df = sampleset.to_pandas_dataframe()
    lowest_energy = np.min(df.iloc[:,num_nodes].values)
    average_energy = np.sum(df.iloc[:,num_nodes].values)/numruns
    print(seed,"  time {0:<8f}  lowest energy {1:<8f} average energy {2:<8f}   exact energy {3:<8f}".format((end_time-start_time)/numruns,lowest_energy,average_energy,exact_gs_energy))

    
    