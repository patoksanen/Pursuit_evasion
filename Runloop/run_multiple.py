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
import time
import json

# Run several games

n_runs = 100 # Amount of games to run
game_results = [] # List to store duration of all the simulated games
max_steps = 3000 # The maximum amount of steps in each simulated game
notes = "Evader uses naive strategy, capture range 2, obstacle configuration square formation circles r=5, distance-variant distance measurements, max steps 3000, v=0.5, n_p=6, wt=1, wa=10" # Optional notes about the particular simulation which will be saved to the results file

# Set time reference
start_time = time.time()

# Announce start
print(f"Running {n_runs} simulations...")

for i in range(n_runs):
    # Initialize a new game
    # evader_x, evader_y = 50, 50
    # evader = Evader(evader_x, evader_y, speed = 0.2)
    # pursuers = [Pursuer(5, 5, speed=0.2, target=evader), Pursuer(40, 5, speed=0.2, target=evader), Pursuer(5, 20, speed=0.2, target=evader), Pursuer(5, 95, speed=0.2, target=evader), Pursuer(95, 5, speed=0.2, target=evader), Pursuer(95, 95, speed = 0.3, target=evader)]
    obstacle1 = Circle((30, 30), 5)  # Bottom left
    obstacle2 = Circle((70, 30), 5)  # Bottom right
    obstacle3 = Circle((30, 70), 5)  # Top left
    obstacle4 = Circle((70, 70), 5)  # Top right
    obstaclebound = Boundary(100,100)
    obstacles = [obstacle1, obstacle2, obstacle3, obstacle4, obstaclebound]

    ## Generate random evader and pursuer positions
    n_pursuers = 6 # Amount of pursuers
    speed = 0.5 # Speed of pursuers and evaders
    pursuers = [] # Initialize empty list to be used later

    # Generate evader position
    valid_point_found = False
    while not valid_point_found:
        x_prop = random.random()*100
        y_prop = random.random()*100
        collisions = 0
        for o in obstacles:
            collisions += o.isInside( (x_prop,y_prop) )
        if collisions == 0:
            valid_point_found = True
    evader = Evader(x_prop, y_prop, speed = speed)

    # Generate pursuer positions
    for j in range(n_pursuers):
        valid_point_found = False
        while not valid_point_found:
            x_prop = random.random()*100
            y_prop = random.random()*100
            collisions = 0
            for o in obstacles:
                collisions += o.isInside( (x_prop,y_prop) )
            if collisions == 0:
                valid_point_found = True
        pursuers.append(Pursuer(x_prop, y_prop, speed = speed, target=evader))


    game = Game(evader, pursuers, max_steps=max_steps, obstacles=obstacles)
    game.run()

    # Evaluate the results of the game
    length = len(game.history)
    game_results.append(length)

    # Update average time
    avg_time = (time.time()-start_time)/(i+1)

    # Print status message
    print(f"Simulation {i+1}/{n_runs} completed")
    print(f"Estimated remaining time: {round( ((avg_time*(n_runs-i-1))/60), 1)} minutes")

print("Done!")
print(game_results)

## Save results to file (in json format)

# Package information with timestamp and notes
filesave = {"Time":time.ctime(),"Notes":notes,"Results":game_results}

# Write to file (If no file exists it will be created. I don't really know where though, but it should appear somewhere close to this file :) .)
with open("simulation_results.txt","a") as f:
    f.write(json.dumps(filesave)+"\n")
