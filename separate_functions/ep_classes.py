import numpy as np
from scipy.spatial import Voronoi
from osbtacle_classes import *

class Agent:
    def __init__(self, x, y, speed):
        self.x = x
        self.y = y
        self.speed = speed

    def move(self, data, obstacles):
        target = self.target_vector(data)
        if np.linalg.norm(target) > 0:
            tarv = target / np.linalg.norm(target)
        else:
            tarv = np.array([0,1])
        avoidance = self.getRepulsiveForce(tarv,obstacles)

        # Define weights for both vectors
        wt = 1
        wa = 20
        control = wt*target + wa*avoidance

        # Enforce max speed
        proposedspeed = np.linalg.norm(control)
        if proposedspeed > self.speed:
            control = control*self.speed / proposedspeed

        dx = control[0]
        dy = control[1]
        self.x += dx
        self.y += dy

    def target_vector(self,data):
        pass

    def avoid_collision(self, obstacles):
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

    # Return avoidance vector

        return avoid_x, avoid_y
    
    def getDistance(self, v, obstacles, step = 1):
        v = v/np.linalg.norm(v) # Make sure v is a unit vector
        p = np.array([self.x,self.y])
        pc = p
        d = step
        blocked = False
        while not blocked:
            d += step
            pc = p+d*v
            for o in obstacles:
                blocked += o.isInside(pc) # Check if we have reached an obstacle.
        return d
    
    def getClosestDistance(self,agents):
        d = float('inf')
        for a in agents:
            ad = np.hypot(self.x-a.x, self.y-a.y)
            ismin = (ad < d)
            d = ismin*ad + (1-ismin)*d
        return d
    
    def getRepulsiveForce(self,u,obstacles):
        """Get the repulsive force acting on the agent. u should be a unit vector."""
        g = lambda x:1/((x+1e-1)**2) # The distance to force function
        alpha = 1 # Scale factor

        dirs = 8
        repforce = 0

        angle = 2*np.pi/dirs
        rotmat = np.array([[np.cos(angle), -np.sin(angle)],
                        [np.sin(angle), np.cos(angle)]])

        for i in range(dirs):
            R = np.linalg.matrix_power(rotmat,i)
            v = R @ u
            dist = self.getDistance(v,obstacles)
            repforce += -v*g(dist)
        
        return alpha*repforce

class Evader(Agent):

    def target_vector(self,pursuers):
        if len(pursuers) <= 3:
            if pursuers:
                closest_pursuer = min(pursuers, key=lambda p: np.hypot(p.x - self.x, p.y - self.y))
                direction = np.array([self.x - closest_pursuer.x, self.y - closest_pursuer.y])
            return direction

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
        return direction

    def move_old(self, pursuers, obstacles):
        target = self.target_vector(pursuers)
        tarv = target / np.linalg.norm(target)
        avoidance = self.getRepulsiveForce(tarv,obstacles)

        # Define weights for both vectors
        wt = 1
        wa = 1
        control = wt*target + wa*avoidance

        # Enforce max speed
        proposedspeed = np.linalg.norm(control)
        if proposedspeed > self.speed:
            control = control*self.speed / proposedspeed

        dx = control[0]
        dy = control[1]
        self.x += dx
        self.y += dy

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

    def target_vector(self, targetpoint):
        direction = np.array(targetpoint) - np.array([self.x, self.y])
        return direction

    def move_old(self, pursuers, obstacles):
        target = self.target_vector(pursuers)
        tarv = target / np.linalg.norm(target)
        avoidance = self.getRepulsiveForce(tarv,obstacles)

        # Define weights for both vectors
        wt = 1
        wa = 1
        control = wt*target + wa*avoidance

        # Enforce max speed
        proposedspeed = np.linalg.norm(control)
        if proposedspeed > self.speed:
            control = control*self.speed / proposedspeed

        dx = control[0]
        dy = control[1]
        self.x += dx
        self.y += dy

    def move_towards_point(self, point, obstacles):
        direction = np.array(point) - np.array([self.x, self.y])
        norm = np.linalg.norm(direction)
    
        if norm > 0:
            direction = direction / norm  # Normalize direction

    # Move only if within bounds
        new_x = self.x + self.speed * direction[0]
        new_y = self.y + self.speed * direction[1]

        if 0 <= new_x <= 100 and 0 <= new_y <= 100:  # Adjust if needed
            self.x = new_x
            self.y = new_y