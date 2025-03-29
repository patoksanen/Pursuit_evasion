from ep_classes2 import *
from obstacle_classes2 import *
import matplotlib.pyplot as plt
import numpy as np

agent = Evader(50,50,1)

obstacles = [Boundary(90,90),Circle((40,40),5)]

dirs = 4
angle = 2*np.pi/dirs
v = np.array([0,1])

rotmat = np.array([[np.cos(angle), -np.sin(angle)],
                   [np.sin(angle), np.cos(angle)]])

for i in range(dirs):
    R = np.linalg.matrix_power(rotmat,i)
    u = R @ v
    print(agent.getDistance(u,obstacles))

repvec = 300*agent.getRepulsiveForce(v,obstacles)
print(repvec)
plt.show()

fig, ax = plt.subplots()
plt.plot(agent.x,agent.y,"go")
plt.plot([agent.x,agent.x+repvec[0]],[agent.y,agent.y+repvec[1]])
for o in obstacles:
    ax.add_patch(o.getPatch())
plt.xlim(0,100)
plt.ylim(0,100)
plt.show()