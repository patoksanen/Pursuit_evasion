import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi

def mirrorx(p, center = 0): # Mirror across the y-axis or a paralell line (center = X)
    pxnew = 2*center-p[0]
    return np.array([pxnew, p[1]])

def mirrory(p, center = 0):
    pynew = 2*center-p[1]
    return np.array([p[0], pynew])  


def polycenter(vs): # vs is list of vertices describing a polygon in 2D - space, returns numpy array corresponding to the polygon centroid
    cx = 0
    cy = 0
    A = 0
    for i in range(len(vs)-1):
        cx += (vs[i][0]+vs[i+1][0])*(vs[i][0]*vs[i+1][1] - vs[i+1][0]*vs[i][1])
        cy += (vs[i][1]+vs[i+1][1])*(vs[i][0]*vs[i+1][1] - vs[i+1][0]*vs[i][1])
        A += (1/2) * (vs[i][0]*vs[i+1][1] - vs[i+1][0]*vs[i][1])
    cx = cx/(6*A)
    cy = cy/(6*A)
    return np.array([cx, cy])

def voronoiCentroids_points(points, buffer = 1): # Takes a list of points and returns a list of corresponding centroids
    # Generate bounds for Voronoi cells
    xmax = 0
    xmin = 0
    ymax = 0
    ymin = 0
    for p in points:
        px = p[0]
        py = p[1]
        xmax = (px>xmax)*px + (px<=xmax)*xmax
        xmin = (px<xmin)*px + (px>=xmin)*xmin
        ymax = (py>ymax)*py + (py<=ymax)*ymax
        ymin = (py<ymin)*py + (py>=ymin)*ymin

    boundsx = [xmin-buffer,xmax+buffer]
    boundsy = [ymin-buffer,ymax+buffer]

    # Generate mirrored points
    npoints = len(points)
    for b in boundsx:
        for i in range(npoints):
            points.append(mirrorx(points[i],b))

    for b in boundsy:
        for i in range(npoints):
            points.append(mirrory(points[i],b))

    vor = Voronoi(points,qhull_options="Qx")

    # Find and plot centroids
    centroids = []
    for i in range(npoints):
        vtcs = [ vor.vertices[j] for j in vor.regions[vor.point_region[i]] ]
        vtcs.append(vtcs[0]) # Append first vertex to last, "closing" the polygon
        cent = polycenter(vtcs)
        centroids.append(cent)

    return centroids

def voronoiCentroids(agents,buffer = 1): # Takes a list of Agent-class objects and returns a list of corresponding centroids. The treated area is centered around the first object in the list.
    npoints = len(agents)

    points = []
    for a in agents:
        points.append([a.x,a.y])

    # Find bounds of x and y
    buffer = 0.1

    dxmax = 0
    dymax = 0

    centerx = points[0][0]
    centery = points[0][1]

    for p in points[1:]:
        dx = np.abs(p[0] - centerx)
        dy = np.abs(p[1] - centery)
        dxmax = (dx>dxmax)*dx + (dx<=dxmax)*dxmax
        dymax = (dy>dymax)*dy + (dy<=dymax)*dymax
        dall = max(dymax,dxmax)

    #boundsx = [centerx-dxmax-buffer,centerx+dxmax+buffer]
    #boundsy = [centery-dymax-buffer,centery+dymax+buffer]

    boundsx = [centerx-dall-buffer,centerx+dall+buffer]
    boundsy = [centery-dall-buffer,centery+dall+buffer]

    # Generate mirrored points
    for b in boundsx:
        for i in range(npoints):
            points.append(mirrorx(points[i],b))

    for b in boundsy:
        for i in range(npoints):
            points.append(mirrory(points[i],b))

    vor = Voronoi(points,qhull_options="Qx")

    # Find and list centroids
    centroids = []
    for i in range(npoints):
        vtcs = [ vor.vertices[j] for j in vor.regions[vor.point_region[i]] ]
        vtcs.append(vtcs[0]) # Append first vertex to last, "closing" the polygon.
        cent = polycenter(vtcs)
        centroids.append(cent)

    return centroids
        