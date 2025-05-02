import numpy as np
from scipy.spatial import Voronoi
from osbtacle_classes import *
from pursuit_strategies import *
from shapely.geometry import Polygon, box
from shapely.geometry import Point
import voronoi_functions
import ep_functions

class Agent:
    def __init__(self, x, y, speed):
        self.x = x
        self.y = y
        self.speed = speed
        self.history = []
        

    def move_test(self, data, obstacles):
        control = self.target_vector(data)

        control = control/np.linalg.norm(control)
        
        dx = control[0]
        dy = control[1]
        self.x += dx
        self.y += dy

    def move(self, data, obstacles, pursuers):
        pursuerrep = np.array([0,0])
        if isinstance(self, Evader):
            wt = 1
            wa = 10
            wd = 1
        if isinstance(self, Pursuer):
            distance = np.linalg.norm([self.target.x - self.x, self.target.y - self.y])
            self.history.append(distance)
            wt = 1
            wa = 10
            wd = 1
            # for p in pursuers:
            #     if p != self:
            #         distvec = np.array([p.x-self.x,p.y-self.y])
            #         if np.linalg.norm(distvec) < 5:
            #             pursuerrep = pursuerrep - distvec
        # target = self.target_vector(data)
        # target = self.naive_strategy(data) # Use instead for naive strategy, also check that the correct data is provided in the Game-class
        if isinstance(self, Pursuer):
            target = self.naive_strategy(data)
        else:
            target = self.target_vector(data)
        
        if np.linalg.norm(target) > 0:
            tarv = target / np.linalg.norm(target)
        else:
            tarv = np.array([0,1])
       
        avoidance = self.getRepulsiveForce(tarv,obstacles)
      
        # Define weights for both vectors
  
        control = wt*target + wa*avoidance + wd*pursuerrep

        # Enforce max speed
        proposedspeed = np.linalg.norm(control)
        if proposedspeed > self.speed:
            control = control*self.speed / proposedspeed

        dx = control[0]
        dy = control[1]
        #self.x += dx
        #self.y += dy
        new_x = self.x + dx
        new_y = self.y + dy
        collision = any(obs.isInside([new_x, new_y]) for obs in obstacles)

        if not collision:
            self.x = new_x
            self.y = new_y
        else:
            # Optionally: try partial move (e.g., reduce step)
            self.x = self.x  # No-op for clarity
            self.y = self.y
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
    
    def getDistance(self, v, obstacles, step = 0.25):
        v = v/np.linalg.norm(v) # Make sure v is a unit vector
        p = np.array([self.x,self.y])
        pc = p
        d = step
        step_number = 0
        blocked = False
        while not blocked:
        
            d += step*(2**step_number)
            # d += step
            pc = p+d*v
            if len(obstacles) == 0:
                break
            for o in obstacles:
                blocked += o.isInside(pc) # Check if we have reached an obstacle.
            step_number += 1
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
    def __init__(self, x, y, speed):
        super().__init__(x, y, speed)
        self.step_counter = 0
        self.cached_target = np.array([0.0, 1.0])  # Default forward
        self.storedcell = 0

    def target_vector(self, pursuers):
        # Get all Voronoi cells
        vor = voronoi_functions.getVoronoiCells([self] + pursuers)

        # Extract vertices of interesting voronoi cells
        cellvertices = []
        for i in range(len(pursuers)+1):
            cellvertices.append([ vor.vertices[j] for j in vor.regions[vor.point_region[i]] ])
        
        # Get areas of all Voronoi cells
        vorareas = []
        for cell in cellvertices:
            vorareas.append(voronoi_functions.polyarea(cell))

        # Check if evader is in large enough cell
        # if vorareas[0] > 500: # Evader cell has index 0 in the list because we passsed the evader position first when creating the cells
        if vorareas[0] == max(vorareas):
            self.storedcell = 0
            print("Evader free!",vorareas[0])
            # return self.naive_strategy(pursuers)
            return ep_functions.polycenter(cellvertices[0]) - np.array([self.x,self.y])
        
        ## If not in largest cell
        # Find largest cell

        if self.storedcell == 0:
            maxarea = max(vorareas[0:])
            largeindex = vorareas.index(maxarea)
            print("Largest cell: ",largeindex)
            self.storedcell = largeindex

        # Find closest vertex of largest cell
        min_dist = 1e4 # Placeholder
        closestvertex = np.array([0,0]) # Placeholder
        for vertex in cellvertices[self.storedcell]:
            dist = np.linalg.norm(vertex-np.array([self.x,self.y]))
            if dist < min_dist:
                min_dist = dist
                closestvertex = vertex
        
        return closestvertex - np.array([self.x,self.y])
    
    def target_vector2(self, pursuers):
        # Get all Voronoi cells
        vor = voronoi_functions.getVoronoiCells([self] + pursuers)

        # Extract vertices of interesting voronoi cells
        cellvertices = []
        for i in range(len(pursuers)+1):
            cellvertices.append([ vor.vertices[j] for j in vor.regions[vor.point_region[i]] ])
        
        # Get areas of all Voronoi cells
        vorareas = []
        for cell in cellvertices:
            vorareas.append(voronoi_functions.polyarea(cell))

        # Check if evader is in large enough cell
        # if vorareas[0] > 500: # Evader cell has index 0 in the list because we passsed the evader position first when creating the cells
        if self.storedcell != 0:
            if np.linalg.norm(self.storedcell - np.array([self.x,self.y])) < 2:
                self.storedcell = 0
        
        if vorareas[0] == max(vorareas) and self.storedcell == 0:
            print("Evader free!",vorareas[0])
            # return self.naive_strategy(pursuers)
            return ep_functions.polycenter(cellvertices[0]) - np.array([self.x,self.y])
        
        ## If not in largest cell
        # Find largest cell

        if self.storedcell == 0:
            maxarea = max(vorareas[0:])
            largeindex = vorareas.index(maxarea)
            print("Largest cell: ",largeindex)
            self.storedcell = largeindex

            # Find closest vertex of largest cell
            min_dist = 1e4 # Placeholder
            closestvertex = np.array([0,0]) # Placeholder
            for vertex in cellvertices[self.storedcell]:
                dist = np.linalg.norm(vertex-np.array([self.x,self.y]))
                if dist < min_dist:
                    min_dist = dist
                    self.storedcell = vertex # Store the target point
        targetp = self.storedcell
        
        return targetp - np.array([self.x,self.y])

    # Get as far away from pursuers as possible - maximize distance
    def naive_strategy(self,pursuers):
        direction = np.array([0,0])
        for p in pursuers:
            direction = direction - np.array([p.x-self.x,p.y-self.y])
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

    def plot_voronoi(self, pursuers):
        import matplotlib.pyplot as plt
        from scipy.spatial import Voronoi, voronoi_plot_2d

        points = np.array([[p.x, p.y] for p in pursuers] + [[self.x, self.y]])
        vor = Voronoi(points)

        fig = voronoi_plot_2d(vor)
        plt.plot(self.x, self.y, 'ro')  # Evader in red
        for p in pursuers:
            plt.plot(p.x, p.y, 'bo')  # Pursuers in blue
        plt.show()

    def get_clipped_voronoi_region(self, pursuers):
        from scipy.spatial import Voronoi
        from shapely.geometry import Polygon, box

    # Build point set
        points = np.array([[p.x, p.y] for p in pursuers] + [[self.x, self.y]])
        vor = Voronoi(points)
        evader_region_index = vor.point_region[-1]
        region = vor.regions[evader_region_index]

        if -1 in region or len(region) == 0:
            print("Unbounded region!")
            return None  # Infinite or empty region

        polygon = Polygon([vor.vertices[i] for i in region])
        clipped = polygon.intersection(box(0, 0, 100, 100))

        if clipped.is_empty:
            return None

        return clipped
 

class Pursuer(Agent):
    def __init__(self, x, y, speed, target):
        super().__init__(x, y, speed)
        self.target = target

    def target_vector(self, targetpoint):
        direction = np.array(targetpoint) - np.array([self.x, self.y])
        return direction
    
    # Travel directly toward evader
    def naive_strategy(self, evader):
        direction = np.array([evader.x-self.x,evader.y-self.y])
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