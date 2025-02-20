import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from ep_classes import *
from ep_functions import *

# Define functions
def isBlocked(x): # Checks if the given position x is inside an obstacle or outside the region.
    blocked = 0
    blocked += max( max(x > [100,100]) , max(x < [-100,-100]) ) # Check bounds of the area
    blocked += (np.linalg.norm(x) < 10) # Add a "forbidden zone" in the middle with radius 10
    return bool(blocked)

# Generate test agents
n = 5 # Number of agents (1 evader and n-1 pursuers)
speed = 1
evaders = []
pursuers = []

randpos = np.random.rand(2*n)*200-2*n*[100]
print(randpos)
evaders.append(Evader(randpos[0],randpos[1],speed))

for i in range(1,n):
    pursuers.append(Pursuer(randpos[2*i],randpos[2*i+1],speed,evaders[0]))

# Capture range
cr = 5

# Initialize plot
fig, ax = plt.subplots()
c1 = patches.Circle((0,0),10)
ax.add_patch(c1)

#plt.ylim(-100,100)
#plt.xlim(-100,100)

# Simulate
captured = False
#while not captured:
for i in range(300):

    # Move evader
    for e in evaders:
        e.move(pursuers)
        plt.plot(e.x,e.y,"go")
    
    # Move pursuers
    targetpoints = voronoiCentroids(evaders + pursuers) # Get target points
    print("tp:",targetpoints)
    for i in range (len(pursuers)):
        pur = pursuers[i]
        pur.move_towards_point(targetpoints[i+1])
        plt.plot(pur.x,pur.y,"r.")
    
    plt.pause(0.05)

plt.show()