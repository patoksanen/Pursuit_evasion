import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Define functions
def isBlocked(x): # Checks if the given position x is inside an obstacle or outside the region.
    blocked = 0
    blocked += max( max(x > [100,100]) , max(x < [-100,-100]) ) # Check bounds of the area
    # blocked += (np.linalg.norm(x) < 10) # Add a "forbidden zone" in the middle with radius 10
    return bool(blocked)

# Define starting values
E = []
n = 10 # Number of evaders
v = 2 # maximum speed of an evader
d = 20 # Desired distance
alpha = 0.01 # Strength of the "collective force"
plotpause = 0.05 # Time between plots, inverse simulation rate
time = 100 # Simulated time

# Generate random positions for evaders
for i in range(n):
    pos = np.random.rand(2)*200-np.array([100,100])
    E.append(pos)

# Initialize plot
fig, ax = plt.subplots()
# c1 = patches.Circle((0,0),10)
# ax.add_patch(c1)

plt.ylim(-100,100)
plt.xlim(-100,100)

# Plot initial positions
for e in E:
    plt.plot(e[0],e[1],"g.")
plt.pause(plotpause)

# Simulate
for i in range(time):

    # Move evaders
    newE = []
    for e in E:
        u = np.array([0,0])
        for e2 in E:
            if not np.all(e2 == e):
                u = u + alpha*(e2-e)*(1-d/np.linalg.norm(e2-e))

        if np.linalg.norm(u) > v: # Enforce max speed
            u = u*v/np.linalg.norm(u)
        
        newe = e + u

        if isBlocked(newe):
            newe = e
        
        newE.append(newe)
        
        plt.plot(newe[0],newe[1],"g.")

    # Update E
    E = newE
    plt.pause(plotpause)

for e in E:
    plt.plot(e[0],e[1],"ro")

plt.show()