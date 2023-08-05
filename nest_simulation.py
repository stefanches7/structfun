import nest
import numpy as np
import matplotlib.pyplot as plt
import sys
from datetime import datetime
from pynestml.codegeneration.nest_code_generator_utils import NESTCodeGeneratorUtils

nest.local_num_threads = 2

nest.ResetKernel()
A = np.loadtxt("test_adjacency_mtx.txt")
N = np.shape(A)[0]
Pe = 0.8
Ne = int(np.ceil(N*Pe))
Ni = N - Ne
rand_Ne = np.random.rand(Ne)
rand_Ni = np.random.rand(Ni)

module_name = "nestml_d8a117e7d303458eb3abe87afe7f8135_module"
nest.Install(module_name)
neuron_model_name = 'izhikevich_customd8a117e7d303458eb3abe87afe7f8135_nestml__with_stdp_synapse_customd8a117e7d303458eb3abe87afe7f8135_nestml'
synapse_model_name = 'stdp_synapse_customd8a117e7d303458eb3abe87afe7f8135_nestml__with_izhikevich_customd8a117e7d303458eb3abe87afe7f8135_nestml'

nodes_exc = nest.Create(neuron_model_name, Ne, params=[{"a":0.02, "b":0.2, \
                                                        "c": -65+15*(rand_Ne[i])**2, "d": 8-6*(rand_Ne[i])**2} for i in range(Ne)])
nodes_inh = nest.Create(neuron_model_name, Ni, params=[{"a":0.02+0.08*(rand_Ni[i]), \
                                                         "b":0.25-0.05*(rand_Ni[i]), \
                                                        "c": -65, "d": 2} for i in range(Ni)])
MAX_EXC_WEIGHT=10
MAX_INH_WEIGHT=.5

nest.CopyModel(synapse_model_name, "excitatory_syn", {"Wmax" : 100, "Wmin":0})
nest.CopyModel(synapse_model_name, "inhibitory_syn", {"Wmax" : 0, "Wmin":-100})

all_ns = nodes_exc + nodes_inh
for i, pre in enumerate(nodes_exc):
    weights = A[:,i]
    
    nonzero_indices = np.where(weights != 0)[0]
    weights = weights[nonzero_indices]
    post = all_ns[nonzero_indices]
    pre_array = np.ones(len(nonzero_indices), dtype=int) * pre.get('global_id')
    nest.Connect(pre_array, post, conn_spec='one_to_one', syn_spec={"synapse_model": "excitatory_syn", "delay" : np.ones(len(post)), 'weight': weights})

for i, pre in enumerate(nodes_inh):
    weights = A[:,Ne + i]

    nonzero_indices = np.where(weights != 0)[0]
    weights = weights[nonzero_indices]
    post = all_ns[nonzero_indices]
    pre_array = np.ones(len(nonzero_indices), dtype=int) * pre.get('global_id')
    nest.Connect(pre_array, post, conn_spec='one_to_one', syn_spec={"synapse_model": "inhibitory_syn", "delay": np.ones(len(post)), 'weight': weights})

T = 2000*10**3 #ms
dt = 1 #ms

startsim = datetime.now()
nest.Simulate(T)
endsim = datetime.now()
