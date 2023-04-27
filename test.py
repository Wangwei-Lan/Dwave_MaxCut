
import numpy as np
from dwave.samplers import SimulatedAnnealingSampler
solver = SimulatedAnnealingSampler()
import time 


H = np.array([[2,-1,-1,0,0],[-1,2,0,-1,-0],[-1,0,3,-1,-1],[0,-1,-1,3,-1],[0,0,-1,-1,2.0]])
sampleset = solver.sample_qubo(-H,num_reads=10,beta_schedule_type='geometric',num_sweeps=1000,beta_range=[0.001,100000000.0],num_sweeps_per_beta=1,)
print(sampleset.lowest())
#df = lowest.to_pandas_dataframe()