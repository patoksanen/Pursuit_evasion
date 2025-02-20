import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.animation as animation
from ep_classes import *
from ep_functions import *
    
class Game:
    def __init__(self, evader, pursuers, max_steps=100):
        self.evader = evader
        self.pursuers = pursuers
        self.max_steps = max_steps
        self.history = []
    
    def step(self):
        self.evader.move(self.pursuers)
        targetpoints = voronoiCentroids([evader] + pursuers)
        for i in range(len(self.pursuers)):
            self.pursuers[i].move_towards_point(targetpoints[i+1])
    
    # Compute Voronoi diagram
        points = np.array([[p.x, p.y] for p in self.pursuers] + [[self.evader.x, self.evader.y]])
        if len(points) >= 4:  # Voronoi requires at least 4 points
            vor = Voronoi(points)
            vor_vertices = vor.vertices
            vor_regions = [vor.vertices[region] for region in vor.regions if region and -1 not in region]
        else:
            vor_vertices, vor_regions = [], []
    
        self.history.append((self.evader.x, self.evader.y, [(p.x, p.y) for p in self.pursuers], vor_vertices, vor_regions))
    
    def run(self):
        for _ in range(self.max_steps):
            self.step()
            if self.is_captured():
                break
    
    def is_captured(self):
        for pursuer in self.pursuers:
            if np.hypot(self.evader.x - pursuer.x, self.evader.y - pursuer.y) < 1.0:
                return True
        return False
    
    def plot(self):
        plt.figure(figsize=(8, 8))
    
        evader_path = np.array([(h[0], h[1]) for h in self.history])
        pursuer_paths = [np.array([(h[2][i][0], h[2][i][1]) for h in self.history]) for i in range(len(self.pursuers))]
    

        plt.plot(evader_path[:, 0], evader_path[:, 1], label='Evader', color='blue', linewidth=2)
        for i, path in enumerate(pursuer_paths):
            plt.plot(path[:, 0], path[:, 1], label=f'Pursuer {i+1}', linestyle='dashed', color='red')

    # Plot Voronoi cells at last step
        last_step = self.history[-1]
        vor_vertices, vor_regions = last_step[3], last_step[4]
    
        for region in vor_regions:
            polygon = np.array(region)
            plt.fill(polygon[:, 0], polygon[:, 1], edgecolor='black', fill=False, alpha=0.5)

    # Plot final positions
        plt.scatter(self.evader.x, self.evader.y, color='blue', marker='o', label='Evader (Final)')
        for pursuer in self.pursuers:
            plt.scatter(pursuer.x, pursuer.y, color='red', marker='x', label='Pursuer (Final)')

        plt.legend()
        plt.xlim(0, 100)  # Adjust bounds as needed
        plt.ylim(0, 100)
        plt.grid()
        plt.show()

# simple try
evader = Evader(70, 70, speed=0.2)
pursuers = [Pursuer(1, 10, speed=0.2, target=evader), Pursuer(1, 80, speed=0.2, target=evader), Pursuer(80, 80, speed=0.2, target=evader), Pursuer(80, 1, speed=0.2, target=evader), Pursuer(5, 5, speed=0.2, target=evader),  Pursuer(6, 5, speed=0.2, target=evader)  ]
game = Game(evader, pursuers, max_steps=20000)
game.run()
game.plot()

## CHAT GPT VIDEO ##

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_xlim(0, 100)
ax.set_ylim(0, 500)
ax.set_title("Pursuit-Evasion Simulation")

# Plot elements
evader_dot, = ax.plot([], [], 'bo', markersize=8, label="Evader")
pursuer_dots, = ax.plot([], [], 'rx', markersize=6, label="Pursuers")
evader_path, = ax.plot([], [], 'b-', linewidth=1)  # Evader trajectory
voronoi_patches = []

# Initialize function
def init():
    evader_dot.set_data([], [])
    pursuer_dots.set_data([], [])
    evader_path.set_data([], [])
    return evader_dot, pursuer_dots, evader_path

# Update function for animation
def update(frame):
    if frame >= len(game.history):
        return init()
    
    evader_x, evader_y, pursuer_positions, vor_vertices, vor_regions = game.history[frame]
    
    # Ensure evader position is a sequence
    evader_dot.set_data([evader_x], [evader_y])  # Wrap in a list
    
    # Ensure pursuer positions are sequences
    pursuer_xs, pursuer_ys = zip(*pursuer_positions) if pursuer_positions else ([], [])
    pursuer_dots.set_data(pursuer_xs, pursuer_ys)

    # Ensure evader path is a sequence
    evader_traj = np.array([(h[0], h[1]) for h in game.history[:frame+1]])
    if len(evader_traj) > 0:
        evader_path.set_data(evader_traj[:, 0], evader_traj[:, 1])
    
    # Clear old Voronoi regions
    for patch in voronoi_patches:
        patch.remove()
    voronoi_patches.clear()

    # Plot Voronoi regions
    for region in vor_regions:
        polygon = np.array(region)
        if len(polygon) > 0:
            patch = plt.Polygon(polygon, edgecolor='black', fill=False, alpha=0.5)
            ax.add_patch(patch)
            voronoi_patches.append(patch)
    
    return evader_dot, pursuer_dots, evader_path, *voronoi_patches

# Create animation
ani = animation.FuncAnimation(fig, update, frames=len(game.history), init_func=init, blit=False, interval=50)


plt.legend()
plt.show()
