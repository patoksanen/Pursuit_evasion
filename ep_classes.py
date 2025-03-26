import numpy as np
from scipy.spatial import Voronoi
from osbtacle_classes import *

class Agent:
    def __init__(self, x, y, speed):
        self.x = x
        self.y = y
        self.speed = speed

    def move(self):
        pass

    def avoid_collision(self, direction_x, direction_y, obstacles):
        """
    Avoid obstacles by applying a repulsive force if too close.
    """
        rectangle_repulsion_factor = 300000
        circle_repulsion_factor = 1000
        min_distance = 0.5
        avoid_x, avoid_y = 0, 0

        for obs in obstacles:
         
            if isinstance(obs, Rectangle):
                closest_x = np.clip(self.x, obs.ll[0], obs.ur[0])
                closest_y = np.clip(self.y, obs.ll[1], obs.ur[1])
                dx = self.x - closest_x
                dy = self.y - closest_y
                dist = np.hypot(dx, dy) + 1e-10
                
                if dist < min_distance:
                    
                    avoid_x += -(dx / dist) * rectangle_repulsion_factor / dist
                    avoid_y += -(dy / dist) * rectangle_repulsion_factor / dist
                   
                # If inside rectangle, push away from center
                if obs.ll[0] <= self.x <= obs.ur[0] and obs.ll[1] <= self.y <= obs.ur[1]:
                    dx = self.x - (obs.ll[0] + obs.ur[0]) / 2
                    dy = self.y - (obs.ll[1] + obs.ur[1]) / 2
                    dist = np.hypot(dx, dy) + 1e-5
                    print("h")

                    avoid_x += (dx / dist) * rectangle_repulsion_factor / dist
                    avoid_y += (dy / dist) * rectangle_repulsion_factor / dist

            elif isinstance(obs, Circle):
                dx = self.x - obs.c[0]
                dy = self.y - obs.c[1]
                dist = np.hypot(dx, dy) + 1e-5  

                if dist < obs.r + min_distance:
                    avoid_x += (dx / dist) * circle_repulsion_factor / dist
                    avoid_y += (dy / dist) * circle_repulsion_factor / dist

    # Modify direction with avoidance
        new_direction_x = direction_x + avoid_x
        new_direction_y = direction_y + avoid_y
        norm = np.hypot(new_direction_x, new_direction_y)
    
        if norm > 0:
            new_direction_x /= norm
            new_direction_y /= norm

        return new_direction_x, new_direction_y

class Evader(Agent):
    def move(self, pursuers, obstacles):
        if len(pursuers) <= 3:
            if pursuers:
                closest_pursuer = min(pursuers, key=lambda p: np.hypot(p.x - self.x, p.y - self.y))
                direction = np.array([self.x - closest_pursuer.x, self.y - closest_pursuer.y])
                norm = np.linalg.norm(direction)
                if norm > 0:
                    direction = direction / norm
                    direction[0], direction[1] = self.avoid_collision(direction[0], direction[1], obstacles)
                    self.move_towards(direction[0], direction[1])
            return

        points = np.array([[p.x, p.y] for p in pursuers] + [[self.x, self.y]])
        vor = Voronoi(points)

        max_area = 0
        best_move = np.array([self.x, self.y])
        for region in vor.regions:
            if not region or -1 in region:
                continue

            polygon = np.array([vor.vertices[i] for i in region])
            area = 0.5 * np.abs(np.dot(polygon[:, 0], np.roll(polygon[:, 1], 1)) - np.dot(polygon[:, 1], np.roll(polygon[:, 0], 1)))
            if area > max_area:
                max_area = area
                best_move = np.mean(polygon, axis=0)

        direction = best_move - np.array([self.x, self.y])
        norm = np.linalg.norm(direction)
        if norm > 0:
            direction = direction / norm
            direction[0], direction[1] = self.avoid_collision(direction[0], direction[1], obstacles)
            self.move_towards(direction[0], direction[1])

    def move_towards(self, direction_x, direction_y):
        x = self.x + self.speed * direction_x
        y = self.y + self.speed * direction_y
        if 0 <= x <= 100 and 0 <= y <= 100:
            self.x = x
            self.y = y



 

class Pursuer(Agent):
    def __init__(self, x, y, speed, target):
        super().__init__(x, y, speed)
        self.target = target

    def move(self, obstacles=[]):
        # Compute direction towards target
        direction = np.array([self.target.x - self.x, self.target.y - self.y])
        norm = np.linalg.norm(direction)
        if norm > 0:
            direction = direction / norm  # Normalize
            direction[0], direction[1] = self.avoid_collision(direction[0], direction[1], obstacles)
            self.x += self.speed * direction[0]
            self.y += self.speed * direction[1]

        # Keep within bounds
        self.x = np.clip(self.x, 0, 100)
        self.y = np.clip(self.y, 0, 100)

    def move_towards_point(self, point, obstacles):
        direction = np.array(point) - np.array([self.x, self.y])
        norm = np.linalg.norm(direction)
    
        if norm > 0:
            direction = direction / norm  # Normalize direction

        # Apply obstacle avoidance
        direction[0], direction[1] = self.avoid_collision(direction[0], direction[1], obstacles)

    # Move only if within bounds
        new_x = self.x + self.speed * direction[0]
        new_y = self.y + self.speed * direction[1]

        if 0 <= new_x <= 100 and 0 <= new_y <= 100:  # Adjust if needed
            self.x = new_x
            self.y = new_y