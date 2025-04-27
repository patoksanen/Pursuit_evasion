import numpy as np
from scipy.spatial import Voronoi
from osbtacle_classes import *
from pursuit_strategies import *
from shapely.geometry import Polygon, box
from shapely.geometry import Point

class Agent:
    def __init__(self, x, y, speed):
        self.x = x
        self.y = y
        self.speed = speed
        self.history = []
        

    def move(self, data, obstacles):
        if isinstance(self, Evader):
            wt = 1
            wa = 20
        if isinstance(self, Pursuer):
            distance = np.linalg.norm([self.target.x - self.x, self.target.y - self.y])
            self.history.append(distance)
            wt = 1
            wa = 50
        target = self.target_vector(data)
        if np.linalg.norm(target) > 0:
            tarv = target / np.linalg.norm(target)
        else:
            tarv = np.array([0,1])
       
        avoidance = self.getRepulsiveForce(tarv,obstacles)
      
        # Define weights for both vectors
  
        control = wt*target + wa*avoidance

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
    
    def getDistance(self, v, obstacles, step = 1):
        v = v/np.linalg.norm(v) # Make sure v is a unit vector
        p = np.array([self.x,self.y])
        pc = p
        d = step
        blocked = False
        while not blocked:
        
            d += step
            pc = p+d*v
            if len(obstacles) == 0:
                break
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
    def __init__(self, x, y, speed):
        super().__init__(x, y, speed)
        self.step_counter = 0
        self.cached_target = np.array([0.0, 1.0])  # Default forward

    def target_vector(self, pursuers):
        self.step_counter += 1

        if self.step_counter % 10 != 1:
        # Return cached target vector on non-update steps
            return self.cached_target




        if len(pursuers) <= 3:
            if pursuers:
                closest_pursuer = min(pursuers, key=lambda p: np.hypot(p.x - self.x, p.y - self.y))
                direction = np.array([self.x - closest_pursuer.x, self.y - closest_pursuer.y])
                return direction
            else:
                return np.array([0.0, 0.0])  # No pursuers case
    
    # Generate Voronoi diagram
        points = np.array([[p.x, p.y] for p in pursuers] + [[self.x, self.y]])
        vor = Voronoi(points)
        evader_point = np.array([self.x, self.y])
    
        # Get evader's region
        from matplotlib.path import Path

        clipped_region = self.get_clipped_voronoi_region(pursuers)
        evader_point = np.array([self.x, self.y])

        if clipped_region is None:
    # fallback: flee from closest pursuer
            closest_pursuer = min(pursuers, key=lambda p: np.hypot(p.x - self.x, p.y - self.y))
            return np.array([self.x - closest_pursuer.x, self.y - closest_pursuer.y])

# Determine where to move
        in_largest = clipped_region.contains(Point(evader_point))

        if in_largest:
            target = np.array(clipped_region.centroid.coords[0])
        else:
    # Find closest point on boundary
            distances = [np.linalg.norm(np.array(pt) - evader_point) for pt in clipped_region.exterior.coords]
            closest_idx = np.argmin(distances)
            target = np.array(clipped_region.exterior.coords[closest_idx])

        max_area = 0
        largest_region = None
        largest_polygon = None
        areas = 0
        for region in vor.regions:
            if not region or -1 in region:
                continue
            polygon = np.array([vor.vertices[i] for i in region])
            area = 0.5 * np.abs(np.dot(polygon[:, 0], np.roll(polygon[:, 1], 1)) -
                            np.dot(polygon[:, 1], np.roll(polygon[:, 0], 1)))
            areas += area
            if area > max_area:
                max_area = area
                largest_region = region
                largest_polygon = polygon
        #print(areas)
        #print(self.step_counter)
        #if self.step_counter == 421:
            #self.plot_voronoi(pursuers)
        if largest_polygon is None:
        # fallback: flee from closest pursuer
            closest_pursuer = min(pursuers, key=lambda p: np.hypot(p.x - self.x, p.y - self.y))
            return np.array([self.x - closest_pursuer.x, self.y - closest_pursuer.y])
    # Check if evader is inside the largest-area Voronoi cell
        from matplotlib.path import Path
        in_largest = False
        if largest_polygon is not None:
            path = Path(largest_polygon)
            in_largest = path.contains_point(evader_point)

        if in_largest:
        # Move towards the centroid of the largest region
            target = np.mean(largest_polygon, axis=0)
        else:
        # Move towards the closest vertex of the largest region
            dists = np.linalg.norm(largest_polygon - evader_point, axis=1)
            target = largest_polygon[np.argmin(dists)]

        direction = target - evader_point
        #self.plot_voronoi(pursuers)
        alpha = 0.7  # Adjust to control smoothness
        smoothed_target = alpha * direction + (1 - alpha) * self.cached_target
        print(self.step_counter)
        if self.step_counter == 450:
            self.plot_voronoi()
        self.cached_target = smoothed_target
        return smoothed_target

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