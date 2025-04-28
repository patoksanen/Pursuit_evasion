import numpy as np
from scipy.spatial import Voronoi
from osbtacle_classes import *

class PursuitStrategy:
    def compute_target(self, pursuer, evader):
        raise NotImplementedError("Override this in subclasses")
    

class DirectPursuit(PursuitStrategy):
    def compute_target(self, pursuer, evader):
        return np.array([evader.x, evader.y]) - np.array([pursuer.x, pursuer.y]) 



    

