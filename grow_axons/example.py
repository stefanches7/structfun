import numpy as np
import datetime
import grow_axons

root_output_folder = "C:\\Users\\stefa\\Documents\\Masterarbeit\\grow_axons\\test_colonies"

# parameters:
r   = 1.5   # culture radius (mm)
rho = 400   # neuron density (neurons mm^-2)
Lmu = 1.0   # mean axon length

r_soma = 7.5e-3 # soma size (mm)

# place neurons:
L = Lmu / np.sqrt(np.pi/2.0)

M = int(np.pi*r**2*rho)
X,Y = np.zeros(M),np.zeros(M)

A = r*np.sqrt(np.random.rand())
theta = 2*np.pi*np.random.rand()
X[0],Y[0] = r+A*np.cos(theta),r+A*np.sin(theta)
for i in range(1,M):
    X[i],Y[i] = X[i-1],Y[i-1]
    while np.any(np.sqrt(np.power(X[:i]-X[i],2)+np.power(Y[:i]-Y[i],2)) < r_soma):  # any neuron overlaps?
        A = r*np.sqrt(np.random.rand())
        theta = 2*np.pi*np.random.rand()
        X[i],Y[i] = r+A*np.cos(theta),r+A*np.sin(theta)

W,_,Xi,Yi = grow_axons.grow_NC_grid(  X, Y,
                                    Pe = 0.8,                   # fraction of excitatory neurons
                                    alphaE = 0.5, alphaI = 0.5, # connectivity prob.
                                    L_mu_E = L, L_mu_I = L,     # mean exc/inh axon length
                                    Dl = 10e-3,                 # axon segment length
                                    phi_sd = 0.1,               # axon segment angle 
                                    r_d_mu_E =150e-3,           # excitatory dentritic tree radius (mm)
                                    r_d_mu_I =150e-3)           # inhibitory dentritic tree radius (mm)

timestamp = str(datetime.datetime.now())

def timestamped_txt(prefix):
    return "".join([prefix, timestamp, ".txt"]).replace(" ", "_").replace(":", "-")

np.savetxt("\\".join([root_output_folder, timestamped_txt("W_")]),W.astype(int))
np.savetxt("\\".join([root_output_folder, timestamped_txt("Xi_")]),Xi.astype(np.single))
np.savetxt("\\".join([root_output_folder, timestamped_txt("Yi_")]),Yi.astype(np.single))
