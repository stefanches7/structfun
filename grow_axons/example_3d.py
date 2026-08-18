import numpy as np
import datetime
import grow_axons_3d

root_output_folder = "C:\\Users\\stefa\\Documents\\Masterarbeit\\grow_axons\\test_colonies"

# parameters:
r   = 1   # culture radius (mm)
rho = 300   # neuron density (neurons mm^-2)
Lmu = 1.0   # mean axon length

r_soma = 7.5e-3 # soma size (mm)

# place neurons:
L = Lmu / np.sqrt(np.pi/2.0)

M = int(np.pi*r**2*rho)
X,Y, Z = np.zeros(M),np.zeros(M), np.zeros(M)

def sample_point_in_sphere(radius):
    A = radius * np.cbrt(np.random.rand())
    theta = 2 * np.pi * np.random.rand()      # azimuth
    cos_phi = 2 * np.random.rand() - 1        # cos(polar)
    sin_phi = np.sqrt(1.0 - cos_phi**2)
    return A, theta, cos_phi, sin_phi

A, theta, cos_phi, sin_phi = sample_point_in_sphere(r)
X[0] = r + A * sin_phi * np.cos(theta)
Y[0] = r + A * sin_phi * np.sin(theta)
Z[0] = r + A * cos_phi

for i in range(1,M):
    X[i],Y[i], Z[i] = X[i-1],Y[i-1], Z[i-1]
    while np.any(np.sqrt(np.power(X[:i]-X[i],2)+np.power(Y[:i]-Y[i],2) + np.power(Z[:i]-Z[i],2)) < r_soma):  # any neuron overlaps?
        A, theta, cos_phi, sin_phi = sample_point_in_sphere(r)
        X[i] = r + A * sin_phi * np.cos(theta)
        Y[i] = r + A * sin_phi * np.sin(theta)
        Z[i] = r + A * cos_phi
      

W,_,Xi,Yi,Zi = grow_axons_3d.grow_NC_grid_3d(  X, Y, Z,
                                    Pe = 0.8,                   # fraction of excitatory neurons
                                    alphaE = 0.4, alphaI = 0.2, # connectivity prob.
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
np.savetxt("\\".join([root_output_folder, timestamped_txt("Zi_")]),Zi.astype(np.single))
