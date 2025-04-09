import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.animation as animation
from ep_classes import *
from ep_functions import *
from osbtacle_classes import *
import random


class Game:
    def __init__(self, evader, pursuers, max_steps=100, field_size=(100, 100), obstacles=[]):
        self.evader = evader
        self.pursuers = pursuers
        self.max_steps = max_steps
        self.history = []
        self.field_size = field_size  # Define the playing field size
        self.obstacles = obstacles  # Store obstacles

    def enforce_boundaries(self, agent):
        """Clamp agent's position within the field boundaries."""
        agent.x = np.clip(agent.x, 0, self.field_size[0])
        agent.y = np.clip(agent.y, 0, self.field_size[1])

    def step(self):
        self.evader.move(self.pursuers, self.obstacles)  # Evader avoids obstacles
        self.enforce_boundaries(self.evader)

        targetpoints = voronoiCentroids([self.evader] + self.pursuers)

        for i, pursuer in enumerate(self.pursuers):
            pursuer.move(targetpoints[i+1], self.obstacles)  # Pass obstacles
            self.enforce_boundaries(pursuer)


        # Compute Voronoi diagram
        points = np.array([[p.x, p.y] for p in self.pursuers] + [[self.evader.x, self.evader.y]])
        if len(points) >= 4:
            vor = Voronoi(points)
            vor_vertices = vor.vertices
            vor_regions = [vor.vertices[region] for region in vor.regions if region and -1 not in region]
        else:
            vor_vertices, vor_regions = [], []

        self.history.append((self.evader.x, self.evader.y, [(p.x, p.y) for p in self.pursuers], vor_vertices, vor_regions))

        for pursuer in self.pursuers:
            old_x, old_y = pursuer.x, pursuer.y
            pursuer.move_towards_point(targetpoints[i+1], obstacles)
            if not is_valid_move(pursuer, self.obstacles):
                pursuer.x, pursuer.y = old_x, old_y  # Revert movement if inside an obstacle

    def run(self):
        for _ in range(self.max_steps):
            self.step()
            if self.is_captured():
                break
    def is_captured(self):
        for pursuer in self.pursuers:
            if np.hypot(self.evader.x - pursuer.x, self.evader.y - pursuer.y) < 3.0:
                return True
        return False

    def plot(self):
        plt.figure(figsize=(8, 8))
        ax = plt.gca()

    # Draw boundary
        boundary_patch = patches.Rectangle((0, 0), self.field_size[0], self.field_size[1], 
                                       fill=False, edgecolor='black', linewidth=2)
        ax.add_patch(boundary_patch)

    # Draw obstacles
        for obs in self.obstacles:
            ax.add_patch(obs.getPatch())

    # Plot evader path
        evader_path = np.array([(h[0], h[1]) for h in self.history])
        plt.plot(evader_path[:, 0], evader_path[:, 1], label='Evader Path', color='blue', linewidth=2)

    # Plot pursuer paths
        pursuer_paths = [np.array([(h[2][i][0], h[2][i][1]) for h in self.history]) 
                        for i in range(len(self.pursuers))]
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
        plt.xlim(0, self.field_size[0])
        plt.ylim(0, self.field_size[1])
        plt.grid()
        plt.show()

def generate_obstacles(num_circles=5, num_rectangles=5, field_size=(100, 100)):
    obstacles = []
    for _ in range(num_circles):
        center = (random.uniform(10, 90), random.uniform(10, 90))  # Avoid edges
        radius = random.uniform(5, 10)
        obstacles.append(Circle(center, radius))

    for _ in range(num_rectangles):
        lowerleft = (random.uniform(10, 80), random.uniform(10, 80))
        width = random.uniform(5, 15)
        height = random.uniform(5, 15)
        upperright = (lowerleft[0] + width, lowerleft[1] + height)
        obstacles.append(Rectangle(upperright, lowerleft))

    return obstacles

def is_valid_move(agent, obstacles):
    """Check if an agent's new position collides with any obstacle."""
    return all(not obs.isInside((agent.x, agent.y)) for obs in obstacles)


#simple try
evader = Evader(70, 70, speed=0.3)
pursuers = [Pursuer(1, 90, speed = 0.31, target = evader), Pursuer(1, 10, speed=0.31, target=evader), Pursuer(1, 20, speed=0.31, target=evader), Pursuer(80, 10, speed=0.21, target=evader), Pursuer(2, 10, speed=0.21, target=evader), Pursuer(90, 10, speed=0.21, target=evader)]
# Define non-overlapping obstacles
obstacle1 = Circle((30, 30), 10)  # Bottom left
obstacle2 = Circle((70, 30), 10)  # Bottom right
obstacle3 = Circle((40, 70), 10)  # Top left
obstacle4 = Circle((80, 70), 10)  # Top right
obstaclebound = Boundary(100,100)

# Store obstacles in a list
obstacles = [obstacle1, obstacle2, obstacle3, obstacle4,obstaclebound]

game = Game(evader, pursuers, max_steps=2000, obstacles=obstacles)
game.run()
game.plot()

#CHAT GPT VIDEO
fig, ax = plt.subplots(figsize=(8, 8))
# Draw obstacles
obstacle_patches = []
for obs in game.obstacles:
    if isinstance(obs, Rectangle):
        patch = plt.Rectangle(obs.ll, obs.ur[0] - obs.ll[0], obs.ur[1] - obs.ll[1], color='blue', alpha=0.5)
    elif isinstance(obs, Circle):
        patch = plt.Circle(obs.c, obs.r, color='blue', alpha=0.5)
    ax.add_patch(patch)
    obstacle_patches.append(patch)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_title("Pursuit-Evasion Simulation")

#Plot elements
evader_dot, = ax.plot([], [], 'bo', markersize=8, label="Evader")
pursuer_dots, = ax.plot([], [], 'rx', markersize=6, label="Pursuers")
evader_path, = ax.plot([], [], 'b-', linewidth=1)  # Evader trajectory
voronoi_patches = []

#Initialize function
def init():
    evader_dot.set_data([], [])
    pursuer_dots.set_data([], [])
    evader_path.set_data([], [])
    return evader_dot, pursuer_dots, evader_path

#Update function for animation
def update(frame):
    if frame >= len(game.history):
        return init()

    evader_x, evader_y, pursuer_positions, vor_vertices, vor_regions = game.history[frame]

#Ensure evader position is a sequence
    evader_dot.set_data([evader_x], [evader_y])  # Wrap in a list

#Ensure pursuer positions are sequences
    pursuer_xs, pursuer_ys = zip(*pursuer_positions) if pursuer_positions else ([], [])
    pursuer_dots.set_data(pursuer_xs, pursuer_ys)

    # Ensure evader path is a sequence
    evader_traj = np.array([(h[0], h[1]) for h in game.history[:frame+1]])
    if len(evader_traj) > 0:
        evader_path.set_data(evader_traj[:, 0], evader_traj[:, 1])

#Clear old Voronoi regions
    for patch in voronoi_patches:
        patch.remove()
    voronoi_patches.clear()
#Plot Voronoi regions
    for region in vor_regions:
        polygon = np.array(region)
        if len(polygon) > 0:
            patch = plt.Polygon(polygon, edgecolor='black', fill=False, alpha=0.5)
            ax.add_patch(patch)
            voronoi_patches.append(patch)

    return evader_dot, pursuer_dots, evader_path, *voronoi_patches

#Create animation
ani = animation.FuncAnimation(fig, update, frames=len(game.history), init_func=init, blit=False, interval=50)


plt.legend()
plt.show()
