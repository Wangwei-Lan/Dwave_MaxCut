
def read_graph(path):
    """
    Read the Max-Cut graph stored in a .txt/.SPARSE file and returns an adjacency matrix of the Max-Cut graph.
    path - filepath of the Max-Cut graph file
    adj_matrix - returns adjacency matrix of Max-Cut instance
    """
    with open(path) as f:
        line = f.readline().split()
        # get number of nodes and number of undirected edges
        N, entries = int(line[0]), int(line[1])
        # create a zeros QUBO matrix
        adj_matrix = np.zeros(shape=(N, N))
        for _ in range(entries):
            # extract the node indices and value; fill the adjacency matrix
            line = f.readline().split()
            node1, node2, value = int(line[0]), int(line[1]), int(line[2])
            # fill both Q[i,j] and Q[j,i] as the QUBO matrix is symmetric
            adj_matrix[node1 - 1, node2 - 1] = value
            adj_matrix[node2 - 1, node1 - 1] = value

        return adj_matrix, N


def formulate_qubo(adj_matrix):
    """
    Take an adjacency matrix of a Max-Cut instance as parameter and generate the QUBO matrix of that instance.
    adj_matrix - adjacency matrix of Max-Cut instance
    qubo - QUBO matrix of Max-Cut instance
    """
    qubo = adj_matrix.copy()
    # get the number of edges of each node by summing the respective row in G
    edge_count = np.sum(qubo, axis=1)
    # fill the diagonal positions with the corresponding edge_count * -1
    np.fill_diagonal(qubo, -edge_count)
    return qubo


def get_MaxCut_qubo_from_file(path):
    return formulate_qubo(read_graph(path)[0])




def formulate_maxcut_qubo_sparse(filename):
    
    #! return -Q , where the sign is introduced to convert maximum problem to minimum problem
    #! The returned qubo matrix is symmetric 
    qubo_dict = {}
    txt_file = open(filename)
    num_nodes,num_edges = txt_file.readline().split()
    
    #! initial qubo sparse dictionary 
    for i in range(int(num_nodes)):
        qubo_dict[(i,i)] = 0.0
     
    for line in txt_file:
        line_split = line.split(' ')
        node1,node2,weight = int(line_split[0])-1,int(line_split[1])-1,float(line_split[2])
        qubo_dict[(node1,node2)] = weight
        qubo_dict[(node2,node1)] = weight
        qubo_dict[(node1,node1)] -=weight
        qubo_dict[(node2,node2)] -=weight
    return int(num_nodes),int(num_edges),qubo_dict 


def formulate_maxcut_ising_sparse(filename):
    
    #! return h,J 
    #! where h is the linear biases of the Ising problems
    #! J does not have on-site terms for max cut problem
    #! 
    h = []
    J = {}
    txt_file = open(filename)
    num_nodes,num_edges = txt_file.readline().split()
    
    #! initial h
    h = np.zeros(int(num_nodes))
    
    for line in txt_file:
        line_split = line.split(' ')
        node1,node2,weight = int(line_split[0])-1,int(line_split[1])-1,float(line_split[2])
        J[(node1,node2)] = 0.25*weight
        J[(node2,node1)] = 0.25*weight
        h[node1] -=0.25*weight
        h[node2] -=0.25*weight

    return h,J



def get_cut_number(spin_vector,edges):
    
    cut_numer = 0
    for pair in edges.keys():
        if spin_vector[pair[0]] != spin_vector[pair[1]]:
            cut_numer += edges[pair]
    return cut_numer


import numpy as np
import networkx as nx 
import os
import random 
import dimod 
import time 

from dwave.samplers import SteepestDescentSolver
sd_solver = SteepestDescentSolver()
from dwave.samplers import SimulatedAnnealingSampler
sa_solver = SimulatedAnnealingSampler()
from dwave.samplers import TabuSampler 
tabu_solver = TabuSampler() 
from dwave.samplers import TreeDecompositionSolver 
tree_solver = TreeDecompositionSolver()

from dwave.system import LeapHybridSampler
hybrid_solver = LeapHybridSampler(solver={'category':'hybrid'})





all_instances = ['G0']#,'G12','G13','G14','G28','G35',
                      # 'G41','G50','G55','G56','G60',
                      # 'G61','G66','G67','G70','G72','G81']
# file contains edge information
file_path = os.getcwd()
#filename = file_path+"/myGenerate_100000_0.07"
#filename = file_path+"/myGenerate_10000_5.txt"
#filename = file_path+"/G81"

print("start calculations!")
#for instance_name in all_instances:

filename = file_path+"/"+instance_name

# construct instance graph from txt file 
# (1) solve qubo model 
num_nodes,num_edges,instance_dict = formulate_maxcut_qubo_sparse(filename)
running_time = []
#bqm = dimod.BQM.from_qubo(instance_dict)
start_time = time.time()
numreads = 1
#sampleset = sa_solver.sample_qubo(instance_dict,num_reads=20,beta_range=[0.01,10000],beta_schedule_type='geometric',num_sweeps=1000,num_sweeps_per_beta=1)
#sampleset = tabu_solver.sample_qubo(instance_dict,num_reads=20,tenure=100,num_restarts=100,coefficient_z_first=10000000,lower_bound_z=5000000)
#bqm = dimod.BQM.from_qubo(instance_dict)
#sampleset = sd_solver.sample_qubo(instance_dict,num_reads=numreads,large_sparse_opt=True)
bqm = dimod.BQM.from_qubo(instance_dict)
sampleset = hybrid_solver.sample(bqm)
end_time = time.time()
df = sampleset.to_pandas_dataframe()
lowest_energy = np.min(df.iloc[:,num_nodes].values)
average_energy = np.sum(df.iloc[:,num_nodes].values)/numreads
print(sampleset.info)
print("instance: ",instance_name,"  time {0:<8f}  lowest energy {1:<8} average energy {2:<8}".format((end_time-start_time)/numreads,lowest_energy,average_energy))


