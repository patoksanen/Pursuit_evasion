import numpy as np
import matplotlib.patches as patches

class Obstacle:

    def __init__(self):
        pass

    def isInside(self,point):
        pass

    def getPatch(self):
        pass

class Circle(Obstacle):

    def __init__(self,center,radius):
        self.c = np.array(center)
        self.r = radius

    def isInside(self,p):
        return (np.linalg.norm(self.c - np.array(p)) < self.r)
    
    def getPatch(self):
        return patches.Circle(self.c,self.r)
    
class Rectangle(Obstacle):

    def __init__(self,upperright,lowerleft):
        self.ur = np.array(upperright)
        self.ll = np.array(lowerleft)
        self.w = upperright[0]-lowerleft[0]
        self.h = upperright[1]-lowerleft[1]

    def isInside(self, p):
        ret = np.all( (np.array(p) > self.ll) * (np.array(p) < self.ur) )
        return ret
    
    def getPatch(self):
        return patches.Rectangle(self.ll,self.w,self.h)