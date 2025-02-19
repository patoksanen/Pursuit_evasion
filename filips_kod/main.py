import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# Define functions
def isBlocked(x): # Checks if the given position x is inside an obstacle or outside the region.
    blocked = 0
    blocked += max( max(x > [100,100]) , max(x < [-100,-100]) ) # Check bounds of the area
    blocked += (np.linalg.norm(x) < 10) # Add a "forbidden zone" in the middle with radius 10
    return bool(blocked)

# Define starting values
x = np.array([50,50])
y = np.array([-50,-50])
vx = 1      # Speed of x (pursuer)
vy = 5     # Speed of y (evader)

# Capture range
cr = 5

# Initialize plot
fig, ax = plt.subplots()
c1 = patches.Circle((0,0),10)
ax.add_patch(c1)

plt.ylim(-100,100)
plt.xlim(-100,100)

# Simulate
while np.linalg.norm(y-x) > cr:
#for i in range(300):

    # Move y (evader)
    uy = np.random.rand(2)*2-np.array([1,1])
    uy = uy*vy / np.linalg.norm(uy)
    newy = y + uy

    # Check if move is valid
    if not isBlocked(newy):
        y = newy

    # Move x (pursuer)
    ux = y-x # Move toward y
    ux = ux*vx / np.linalg.norm(ux)
    newx = x+ux

    # Check if move is valid
    if not isBlocked(newx):
        x = newx

    plt.plot(x[0],x[1],"r.")
    plt.plot(y[0],y[1],"g.")
    
    plt.pause(0.05)

plt.show()