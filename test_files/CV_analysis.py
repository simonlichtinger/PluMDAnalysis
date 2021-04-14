""" This script loads a trajectory, and performs a PCA of which distances best describe the RMSD trend. """

import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.transformations import unwrap

import PluMDAnalysis.PluMDAnalysis

""" Calculate the angle between two lines, defined as tuples of MDAnalyses Atomgroup objects.
:params:
    a (tuple<MDAnalysis.AtomGroup>)  ---  Line 1
    b (tuple<MDAnalysis.AtomGroup>)  ---  Line 2
:returns:
    Angle in degrees
"""
def angle_between(a, b):
    vector1 =  a[1].center_of_geometry() - a[0].center_of_geometry()
    l1 = np.sqrt(vector1.dot(vector1))
    vector2 =  b[1].center_of_geometry() - b[0].center_of_geometry()
    l2 = np.sqrt(vector2.dot(vector2))
    return np.arcsin(vector1.dot(vector2) / (l1*l2)) * 360 / 2 / np.pi

""" Calculate the distance between two points, defined as MDAnalyses Atomgroup objects.
:params:
    a (MDAnalysis.AtomGroup)  ---  Selection 1
    b (MDAnalysis.AtomGroup)  ---  Selection 2
:returns:
    Distance in Angstrom
"""
def distance_between(a,b):
    return np.sqrt(np.sum((a.center_of_geometry() - b.center_of_geometry())**2))

# Define files of interest
coordinates = "equil.gro"
trajectory = "prod_align.xtc"
reference = "ref_wrong_res_numbers.pdb"



# How long is a frame in ns
frame_duration = 0.05

# What time periods are of interest?
times_of_interest = [(0,150)]

# Load trajectory and reference coordinates
trj = mda.Universe(coordinates, trajectory)
ref = mda.Universe(reference)
#trj.trajectory.add_transformations(unwrap(trj.select_atoms("protein")))

# Select atoms at end of helices

helices = [
    (trj.select_atoms('resid 68:73 and name CA'),
    trj.select_atoms('resid 46:51 and name CA')),
    (trj.select_atoms('resid 79:84 and name CA'),
    trj.select_atoms('resid 98:103 and name CA')),
    (trj.select_atoms('resid 125:130 and name CA'),
    trj.select_atoms('resid 110:115 and name CA')),
    (trj.select_atoms('resid 143:148 and name CA'),
    trj.select_atoms('resid 166:171 and name CA')),
    (trj.select_atoms('resid 199:204 and name CA'),
    trj.select_atoms('resid 177:182 and name CA')),
    (trj.select_atoms('resid 217:222 and name CA'),
    trj.select_atoms('resid 232:237 and name CA')),
    (trj.select_atoms('resid 325:330 and name CA'),
    trj.select_atoms('resid 290:295 and name CA')),
    (trj.select_atoms('resid 341:346 and name CA'),
    trj.select_atoms('resid 364:369 and name CA')),
    (trj.select_atoms('resid 397:402 and name CA'),
    trj.select_atoms('resid 376:381 and name CA')),
    (trj.select_atoms('resid 609:614 and name CA'),
    trj.select_atoms('resid 631:636 and name CA')),
    (trj.select_atoms('resid 660:665 and name CA'),
    trj.select_atoms('resid 642:647 and name CA')),
    (trj.select_atoms('resid 671:676 and name CA'),
    trj.select_atoms('resid 691:696 and name CA'))
]

rmsd_group = "(resid 68:73 or resid 46:51 or resid 79:84 or resid 98:103 or resid 125:130 or resid 110:115 or resid 143:148 or resid 166:171 or resid 199:204 or resid 177:182 or resid 217:222 or resid 232:237 or resid 325:330 or resid 290:295 \
or resid 341:346 or resid 364:369 or resid 397:402 or resid 376:381 or resid 609:614 or resid 631:636 or resid 660:665 or resid 642:647 or resid 671:676 or resid 691:696 or resid 417) and name CA"

test=trj.select_atoms(rmsd_group)

# Define interesting distances, between all tips and then all bases
distances_of_interest = []
for i in range(12):
    for j in range(i+1,12):
        distances_of_interest.append((helices[i][0], helices[j][0]))
for i in range(12):
    for j in range(i+1,12):
        distances_of_interest.append((helices[i][1], helices[j][1]))

# Gather distance and RMSD information
distances = np.zeros((len(distances_of_interest),len(trj.trajectory)))

# calculate the RMSD of the trajectory
#R = mda.analysis.rms.RMSD(trj, ref, select=rmsd_group)
#R.run()
#rmsds = np.array(R.rmsd.T[2])

# Calculate the interhelical tip/base distances over the trajectory
for time_index, ts in enumerate(trj.trajectory):
    for distance_index, dist in enumerate(distances_of_interest):
        distances[distance_index, time_index] = distance_between(dist[0], dist[1])

# Plot the RMSD progression for confirmation
#plt.plot(np.linspace(0, frame_duration*len(trj.trajectory), len(trj.trajectory)), rmsds)
#plt.xlabel("Time / ns")
#plt.ylabel("RMSD/ A")
#plt.show()

# Split trajectory into parts and set the mean of all distance measurements to zero separately for each block
frames_of_interest = [(x[0] / frame_duration, x[1] / frame_duration) for x in times_of_interest]

distances_parts = [distances[:,int(x[0]):int(x[1])] for x in frames_of_interest]
print(distances_parts, distances_parts[0].shape)
means = [ np.array([np.mean(dist[x,:]) for x in range(dist.shape[0])]).reshape(dist.shape[0],1) for dist in distances_parts]
for i in range(len(distances_parts)):
    distances_parts[i] = distances_parts[i] - means[i]

# compute the PCAs of the distance matrices
cov = [ np.cov(dist) for dist in distances_parts]
eigenvalues, eigenvectors = [], []
for i in range(len(distances_parts)):
    val, vec = np.linalg.eig(cov[i])
    eigenvalues.append(val)
    eigenvectors.append(vec)
print(f"The first principal component explains: {[eigenvalue[0]/sum(eigenvalue)*100 for eigenvalue in eigenvalues]} percent")

# Plot the first PCs
for i in range(len(times_of_interest)):
    plt.plot(range(1,133),eigenvectors[i][0], label = f"Time chunk {times_of_interest[i][0]} to {times_of_interest[i][1]} ns")

# Plot a grid guide as well
points = [1,12,22,31,39,46,52,57,61,64,66,67]
gridpoints = points.copy()
for p in points:
    gridpoints.append(p+66)

for x in gridpoints[1:-1]:
    plt.plot([x-0.5,x-0.5],[-0.5,0.5], color='black',linewidth=0.5)
plt.plot([0,134],[0,0], color='black',linewidth=0.5)

for x in range(1,12):
    pos = (gridpoints[x-1] + gridpoints[x])/2-0.1
    plt.text(pos-0.8,-0.3, str(x), fontweight='bold')
for x in range(1,12):
    pos = (gridpoints[x-1] + gridpoints[x])/2-0.1 + 66
    plt.text(pos-0.8,-0.3, str(x), fontweight='bold')

count = 2
for i in range(12):
    for j in range(i+2,13):
        plt.text(count-1.3,-0.28, str(j), fontsize='xx-small')
        count += 1


count = 2
for i in range(12):
    for j in range(i+2,13):
        plt.text(count-1.3+66,-0.28, str(j), fontsize='xx-small')
        count += 1

# And show it all
plt.ylim(-0.4,0.4)
plt.xlabel("TM helix pair")
plt.ylabel("PC 1 value")
plt.title("First principal component")
plt.legend()
fig = plt.gcf()
fig.set_size_inches(13,10)
plt.show()
