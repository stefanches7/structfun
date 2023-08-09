import nest
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
from pynestml.codegeneration.nest_code_generator_utils import NESTCodeGeneratorUtils

nest.local_num_threads = 14
dt = 1 #ms
nest.SetKernelStatus({"resolution":dt})

Path("spikes").mkdir(parents=True, exist_ok=True)
Path("weights").mkdir(parents=True, exist_ok=True)

spikes_exc_backend="spikes/spikes_exc_%s.tsv" % datetime.now()
spikes_inh_backend="spikes/spikes_inh_%s.tsv" % datetime.now()

# SIMULATION PARAMS
T = 50*10**3 #ms
synscaling_int = 1000 #ms
elapsed_time = 0 #ms
gamma = 0.2
time_washin = T/5
gamma_set = False
time_washoff = 3*T/5
gamma_unset = False

nest.ResetKernel()
A = np.loadtxt("test_adjacency_mtx.txt")
N = np.shape(A)[0]
Pe = 0.8
Ne = int(np.ceil(N*Pe))
Ni = N - Ne
rand_Ne = np.random.rand(Ne)
rand_Ni = np.random.rand(Ni)

module_name = "nestml_7849465c1f8d4551a665e78a5fcebef7_module"
nest.Install(module_name)
neuron_model_name = 'izhikevich_custom7849465c1f8d4551a665e78a5fcebef7_nestml__with_stdp_depression7849465c1f8d4551a665e78a5fcebef7_nestml'
synapse_model_name = 'stdp_depression7849465c1f8d4551a665e78a5fcebef7_nestml__with_izhikevich_custom7849465c1f8d4551a665e78a5fcebef7_nestml'

nodes_exc = nest.Create(neuron_model_name, Ne, params=[{"a":0.02, "b":0.2, \
                                                        "c": -65+15*(rand_Ne[i])**2, "d": 8-6*(rand_Ne[i])**2} for i in range(Ne)])
nodes_inh = nest.Create(neuron_model_name, Ni, params=[{"a":0.02+0.08*(rand_Ni[i]), \
                                                         "b":0.25-0.05*(rand_Ni[i]), \
                                                        "c": -65, "d": 2} for i in range(Ni)])
MAX_EXC_WEIGHT=10
MAX_INH_WEIGHT=.5

wr_exc = nest.Create('weight_recorder')
wr_inh = nest.Create('weight_recorder')
wr_exc.set(record_to='ascii', label='weights/excitatory_%s.tsv' % datetime.now())
wr_inh.set(record_to='ascii', label='weights/inhibitory_%s.tsv' % datetime.now())

nest.CopyModel(synapse_model_name, "excitatory_syn", {"Wmax" : 100, "Wmin":0, "weight_recorder": wr_exc})
nest.CopyModel(synapse_model_name, "inhibitory_syn", {"Wmax" : -100, "Wmin":0, "weight_recorder": wr_inh})

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

NOISE_MAX = 5
noise_exc = nest.Create("noise_generator",Ne, params={"std": NOISE_MAX})
noise_inh = nest.Create("noise_generator", Ni, params={"std": 2})
nest.Connect(noise_exc, nodes_exc, conn_spec = "one_to_one")
nest.Connect(noise_inh, nodes_inh, conn_spec = "one_to_one")
recorder_ex = nest.Create("spike_recorder")
recorder_in = nest.Create("spike_recorder")
recorder_ex.set(record_to="ascii", label=spikes_exc_backend)
recorder_in.set(record_to="ascii", label=spikes_inh_backend)

nest.Connect(nodes_exc, recorder_ex)
nest.Connect(nodes_inh, recorder_in)

w_targets = np.sum(A, axis = 0)
def normalize_weights(neurons_to_be_normalized, w_targets = None):
    if w_targets is None:
        w_targets = np.ones(len(neurons_to_be_normalized))

    for i, neuron in enumerate(neurons_to_be_normalized):
        conn = nest.GetConnections(target=neuron)
        w = np.array(conn.weight)
        w_normed = w / sum(abs(w))  # L1-norm
        conn.weight = w_targets[i] * w_normed

def set_gamma_excitatory(gamma):
    conn = nest.GetConnections(synapse_model = "excitatory_syn")
    conn.set({"gamma":gamma})

startsim = datetime.now()
with nest.RunManager():
    for _ in range(int(T/synscaling_int)):
        nest.Run(synscaling_int)
        elapsed_time += synscaling_int
        if (elapsed_time > time_washin and not gamma_set):
            set_gamma_excitatory(gamma)
            gamma_set = True
        elif (elapsed_time > time_washoff and not gamma_unset):
            set_gamma_excitatory(1)
            gamma_unset = True
        normalize_weights(all_ns, w_targets)

endsim = datetime.now()
