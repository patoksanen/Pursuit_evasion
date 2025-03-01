import numpy as np

# Define a circle
c = ([0,0],1) # c = ( center point, radius )

def isInCircle(p,c):
    d = np.linalg.norm(np.array(p)-np.array(c[0]))
    return d<c[1]

# Define a rectangle
r = ([-10,-10],[10,10]) # r = (left bottom corner, right upper corner)

def isInRect(p,r):
    ret = np.all( (np.array(p) > np.array(r[0])) * (np.array(p) < np.array(r[1])) )
    return ret

# Define a triangle
t = ([-5,0],[5,0],[0,5]) # Vertices

def triangleArea(t):
    s1 = np.array(t[1])-np.array(t[0])
    s2 = np.array(t[2])-np.array(t[0])
    return np.linalg.norm(np.cross(s1,s2))/2

def isInTriangle(p,t):
    subsets = [(0,1),(1,2),(0,2)]
    A = triangleArea(t)
    A2 = 0
    for subset in subsets:
        p1 = t[subset[0]]
        p2 = t[subset[1]]
        t2 = (p, p1, p2)
        A2 += triangleArea(t2)
    return ( A == A2 )


print(isInRect([50,0.5],r))
print(isInTriangle([0,1],t))