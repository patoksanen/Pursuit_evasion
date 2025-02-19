import matplotlib.pyplot as plt
import numpy as np
import scipy
import scipy.spatial
from scipy.spatial import voronoi_plot_2d
from ep_functions import *

points = np.array([
    [1,1],
    [1,-1],
    [-1,0],
    [0.5,2],
    [0,0]
])

points = list(points)

boundsx = [-2,2]
boundsy = [-2,3]

# Generate mirrored points
npoints = len(points)
for b in boundsx:
    for i in range(npoints):
        points.append(mirrorx(points[i],b))

for b in boundsy:
    for i in range(npoints):
        points.append(mirrory(points[i],b))

vor = scipy.spatial.Voronoi(points,qhull_options="Qx")
print("regions:",vor.regions)
print("region13",vor.regions[13])
print("vertices:",vor.vertices)
print("point_region:",vor.point_region)

voronoi_plot_2d(vor)

# Find and plot centroids
for i in range(npoints):
    vtcs = [ vor.vertices[j] for j in vor.regions[vor.point_region[i]] ]
    vtcs.append(vtcs[0])
    cent = polycenter(vtcs)
    plt.plot(cent[0],cent[1],"rx")

plt.show()