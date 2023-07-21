import numpy as np
import matplotlib.pyplot as plt
import glob, os

# plot the profile for some variable at a selection of times
data = np.loadtxt("test2_HR/constraint_norms.dat")

# time data
timedata = data[:,0]
num_times = np.size(timedata)

# constraint data
hamdata = data[:,1]
momdata = data[:,2]

# Plot HR case
plt.semilogy(timedata, hamdata, label="Ham data HR")
plt.semilogy(timedata, momdata,'--', label="Mom data HR")
plt.xlabel('t')
plt.ylabel('value over time of constraints')

# plot the profile for some variable at a selection of times
data = np.loadtxt("test2_LR/constraint_norms.dat")

# time data
timedata = data[:,0]
num_times = np.size(timedata)

# constraint data
hamdata = data[:,1]
momdata = data[:,2]
#ymin = np.min(hamdata)
#ymax = np.max(hamdata)

# plot LR case
plt.semilogy(timedata, hamdata, label="Ham data LR")
plt.semilogy(timedata, momdata, '--', label="Mom data LR")
plt.xlabel('t')
#plt.ylim(1.1*ymin,1.1*ymax)

# tidy plot
plt.legend(loc=4)
plt.grid()
plt.tight_layout()
filename = "Convergence.png"
plt.savefig(filename)