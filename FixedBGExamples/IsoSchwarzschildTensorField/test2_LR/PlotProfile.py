import numpy as np
import matplotlib.pyplot as plt
import glob, os

# plot the profile for some variable at a selection of times
data = np.loadtxt("outputs.dat")

# radial points
rdataline = data[0,:]
num_points = np.size(rdataline) - 1
L = 128
dx = np.sqrt(3.0 * 128.0 * 128.0)/num_points
rdata = np.arange(0,num_points-1) * dx

# time data
timedata = data[:,0]
num_times = np.size(timedata)

outputfile_prefix = "plots/radial_profile"
ymax = np.max(data[:,1:num_points-1])
ymin = np.min(data[:,1:num_points-1])

for i, t_i in enumerate(timedata) :
    labelt = "t="+str(round(t_i,2))
    f_t = data[i, 1:num_points]
    plt.plot(rdata, f_t, label=labelt)
    plt.legend(loc=4)
    plt.xlabel('r')
    #plt.xlim(-0.2,35.0)
    plt.ylim(1.1*ymin,1.1*ymax)
    plt.ylabel('value over time of h_11')
    plt.grid()
    plt.tight_layout()
    filename = outputfile_prefix + ('%04d' % i) + ".png"
    plt.savefig(filename)
    plt.close()

# needs ffmpeg installed
os.system('ffmpeg -r 5 -f image2 -s 1920x1080 -i ' + outputfile_prefix + '%04d.png -crf 25 -pix_fmt yuv420p ' + outputfile_prefix + '.mp4')