import numpy as np

name = '/home/conrad/Dropbox/accretion/positions_density_cut_med_conrad'

data = np.loadtxt(name)

with open(name + '.ply', 'w') as f:
    f.write('ply\n')
    f.write('format ascii 1.0\n')
    f.write('element vertex {}\n'.format(data.shape[0]))
    f.write('property float x\n')
    f.write('property float y\n')
    f.write('property float z\n')
    f.write('end_header\n')

    for i in range(data.shape[0]):
        f.write('{:7.5f} {:7.5f} {:7.5f}\n'.format(data[i,0], data[i,1], data[i,2]))
