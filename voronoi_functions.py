import numpy as np
from scipy.spatial import Voronoi

### Den här filen innehåller 2 nya funktioner:
# polyarea() - Ger arean av en polygon från dess hörn
# getVoronoiCells() - Ger ett scipy-Voronoi objekt med voronoi-celler för givna agenter inom ett begränsat 100x100-område

### Hjälpfunktioner ###
def mirrorx(p, center = 0): # Mirror across the y-axis or a paralell line (center = X)
    pxnew = 2*center-p[0]
    return np.array([pxnew, p[1]])

def mirrory(p, center = 0): # Same but with the x-axis
    pynew = 2*center-p[1]
    return np.array([p[0], pynew])  
### Dessa är kopierade direkt från ep_functions ###

### Nya funktioner ###
def polyarea(vs):
    '''vs is a list of vertices describing a simple polygon in 2D-space, returns area of the polygon'''
    # En modifierad verion av funktionen "polycenter()" i ep_functions
    vs.append(vs[0]) # Append first vertex to last, "closing" the polygon.
    A = 0
    for i in range(len(vs)-1):
        A += (1/2) * (vs[i][0]*vs[i+1][1] - vs[i+1][0]*vs[i][1])
    return abs(A)

def getVoronoiCells(agents,buffer = 1):
    '''Takes a list of Agent-class objects and returns the corresponding voronoi cells restricted to the area [0,100]x[0,100] as a scipy Voronoi-object. This function generates a bunch of additional voronoi-cells, but they do not need to be cared about.'''
    # En modifierad version av funktionen "voronoiCentroids()" i ep_functions
    # Notera att "buffer > 0" ger en något överskattad area för de resulterande voronoicellerna, men kan behövas för att hantera agenter vid kanterna. Prova annars "buffer = 0" för exakta areor.
    npoints = len(agents)

    points = []
    for a in agents:
        points.append([a.x,a.y])

    boundsx = [-buffer,100+buffer]
    boundsy = [-buffer,100+buffer]

    # Generate mirrored points
    for b in boundsx:
        for i in range(npoints):
            points.append(mirrorx(points[i],b))

    for b in boundsy:
        for i in range(npoints):
            points.append(mirrory(points[i],b))

    vor = Voronoi(points,qhull_options="Qx")

    return vor