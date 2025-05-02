import numpy as np
import matplotlib.pyplot as plt
from pursuit_evasion_encirclement import Game, Evader, Pursuer, generate_obstacles
import random
import matplotlib.animation as animation
from osbtacle_classes import *
from evasion_strategies import *
from pursuit_strategies import *
from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib.patches import Polygon

def check_distance(num_simulations = 1, field_size = (100, 100), max_steps = 1000):
    data = []
    evader_x, evader_y = 50, 50
    evader = Evader(evader_x, evader_x, speed = 0.5)
    pursuerpos = [
        (10,5),
        (40,5),
        (60,5),
        (90,5),
        (90,90)
    ]
    pursuers = []
    for pos in pursuerpos:
        pursuers.append(Pursuer(pos[0], pos[1], speed=0.5, target=evader))
    # pursuers = [Pursuer(5, 5, speed=0.5, target=evader), Pursuer(10, 5, speed=0.5, target=evader), Pursuer(15, 7, speed=0.5, target=evader), Pursuer(20, 9, speed=0.5, target=evader), Pursuer(25, 4, speed=0.5, target=evader), Pursuer(30, 1, speed = 0.5, target=evader)]
    obstacles = generate_obstacles(num_circles=1,num_rectangles=0,field_size=field_size)
    obstacle1 = Circle((30, 30), 5)  # Bottom left
    obstacle2 = Circle((70, 30), 5)  # Bottom right
    obstacle3 = Circle((30, 70), 5)  # Top left
    obstacle4 = Circle((70, 70), 5)  # Top right
    obstaclebound = Boundary(100,100)
    obstacles = [obstacle1, obstacle2, obstacle3, obstacle4, obstaclebound]
    #obstacles = [obstacle4, obstaclebound]
    #obstacles = [obstaclebound]
    game = Game(evader, pursuers, max_steps=max_steps, obstacles=obstacles)    
    game.run()
    for pursuer in pursuers:
        data.append(pursuer.history)

    plt.figure(figsize=(10, 5))
    for i, distances in enumerate(data):
        plt.plot(distances, label=f'Pursuer {i+1}')

    plt.xlabel("Timestep")
    plt.ylabel("Distance to Evader")
    plt.title("Pursuer Distance to Evader Over Time")
    plt.legend()
    plt.grid()
    plt.show()

    # ANIMATION PART
    fig, ax = plt.subplots(figsize=(8, 8))

    # Draw obstacles
    for obs in game.obstacles:
        if isinstance(obs, Rectangle):
            patch = plt.Rectangle(obs.ll, obs.ur[0] - obs.ll[0], obs.ur[1] - obs.ll[1], color='blue', alpha=0.5)
        elif isinstance(obs, Circle):
            patch = plt.Circle(obs.c, obs.r, color='blue', alpha=0.5)
        ax.add_patch(patch)

    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_title("Pursuit-Evasion Simulation")

    # Plot elements
    evader_dot, = ax.plot([], [], 'bo', markersize=8, label="Evader")
    pursuer_dots, = ax.plot([], [], 'rx', markersize=6, label="Pursuers")
    evader_path, = ax.plot([], [], 'b-', linewidth=1)  # Evader trajectory
    voronoi_patches = []

    # Initialize function
    def init():
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        evader_dot.set_data([], [])
        pursuer_dots.set_data([], [])
        evader_path.set_data([], [])
        return evader_dot, pursuer_dots, evader_path
#Update function for animation
    def update(frame):
        ax.clear()  # Clear the axes

    # Set limits and title again
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)
        ax.set_title("Pursuit-Evasion Simulation")

    # Redraw obstacles
        for obs in game.obstacles:
            if isinstance(obs, Rectangle):
                patch = plt.Rectangle(obs.ll, obs.ur[0] - obs.ll[0], obs.ur[1] - obs.ll[1], color='blue', alpha=0.5)
            elif isinstance(obs, Circle):
                patch = plt.Circle(obs.c, obs.r, color='blue', alpha=0.5)
            ax.add_patch(patch)

        if frame >= len(game.history):
            return []

    # Get evader and pursuer positions
        evader_x, evader_y, pursuer_positions, _, _ = game.history[frame]
    
    # Plot evader
        ax.plot(evader_x, evader_y, 'ro', markersize=8, label="Evader")  # red circle

    # Plot pursuers
        for x, y in pursuer_positions:
            ax.plot(x, y, 'bx', markersize=6, label="Pursuer")  # blue crosses

    # Plot evader trajectory
        evader_traj = np.array([(h[0], h[1]) for h in game.history[:frame+1]])
        ax.plot(evader_traj[:, 0], evader_traj[:, 1], 'r-', linewidth=1)

    # Avoid duplicate legend entries
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())

        return []


    # Create animation
    ani = animation.FuncAnimation(fig, update, frames=len(game.history), init_func=init, blit=False, interval=50)

    plt.legend()
    plt.show()
















    
check_distance()