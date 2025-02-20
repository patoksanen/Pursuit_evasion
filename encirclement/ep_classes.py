import numpy as np
from scipy.spatial import Voronoi


class Agent:
    def __init__(self, x, y, speed):
        self.x = x
        self.y = y
        self.speed = speed
    
    def move(self):
        pass  # Implemented in Evader/Pursuer class



class Evader(Agent):
    
    def move(self, pursuers):
        if len(pursuers) <= 3:
            # If unable to create voronoi diagram
            if pursuers:
                closest_pursuer = min(pursuers, key=lambda p: np.hypot(p.x - self.x, p.y - self.y)) # Find closes pursuer
                direction = np.array([self.x - closest_pursuer.x, self.y - closest_pursuer.y]) # Direction away
                norm = np.linalg.norm(direction)
                if norm > 0:
                    direction = direction / norm
                    #self.x += self.speed * direction[0]
                    #self.y += self.speed * direction[1]
                    self.move_towards(direction[0], direction[1])
            return  # Skip Voronoi if not enough points
        
        # Create Voronoi diagram
        points = np.array([[p.x, p.y] for p in pursuers] + [[self.x, self.y]])
        vor = Voronoi(points)

        # Find the largest Voronoi cell
        max_area = 0
        best_move = np.array([self.x, self.y])
        for region in vor.regions:
            if not region or -1 in region:
            #if not region:    
                self.adjust_movement(pursuers)
                continue 
                
            polygon = np.array([vor.vertices[i] for i in region])
            area = 0.5 * np.abs(np.dot(polygon[:, 0], np.roll(polygon[:, 1], 1)) - np.dot(polygon[:, 1], np.roll(polygon[:, 0], 1)))
            if area > max_area:
                max_area = area
                best_move = np.mean(polygon, axis=0)
               
        # Move towards the best region
        direction = best_move - np.array([self.x, self.y])
        norm = np.linalg.norm(direction)
        if norm > 0:
            direction = direction / norm
            #self.x += self.speed * direction[0]
            #self.y += self.speed * direction[1]
            self.move_towards(direction[0], direction[1])
    def adjust_movement(self,pursuers):
           # If unable to create voronoi diagram
            if pursuers:
                closest_pursuer = min(pursuers, key=lambda p: np.hypot(p.x - self.x, p.y - self.y)) # Find closes pursuer
                direction = np.array([self.x - closest_pursuer.x, self.y - closest_pursuer.y]) # Direction away
                norm = np.linalg.norm(direction)
                if norm > 0:
                    direction = direction / norm
                    #self.x += self.speed * direction[0]
                    #self.y += self.speed * direction[1]
                    self.move_towards(direction[0], direction[1])
            return  # Skip Voronoi if not enough points
    
    def move_towards(self, direction_x, direction_y):
        x = self.x + self.speed * direction_x
        y = self.y + self.speed * direction_y
        if 0 <= x <= 100 and 0 <= y <= 500:  # Ensure it stays within bounds
            self.x = x
            self.y = y                

class Pursuer(Agent):                                    
    def __init__(self, x, y, speed, target):
        super().__init__(x, y, speed)
        self.target = target
    
    def move(self):
        # pursuit strategy (Simple direction)
        direction = np.array([self.target.x - self.x, self.target.y - self.y])
        norm = np.linalg.norm(direction)
        if norm > 0:
            direction = direction / norm  # Normalize tehe direction
            self.x += self.speed * direction[0]
            self.y += self.speed * direction[1]

        self.y += self.speed * direction[1]

    def move_towards_point(self, point): # Point is array of length 2
        u = np.array(point)-np.array([self.x,self.y])
        u = u*self.speed/np.linalg.norm(u)
        self.x = self.x + u[0]
        self.y = self.y + u[1]
