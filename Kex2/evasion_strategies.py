import numpy as np
from scipy.spatial import Voronoi
from osbtacle_classes import *
from abc import ABC, abstractmethod


class EvaderStrategy(ABC):
    @abstractmethod
    def compute_target(self, evader, pursuers):
        pass


class VoronoiEvade(EvaderStrategy):
    def compute_target(self, evader, pursuers):
        points = np.array([[p.x, p.y] for p in pursuers] + [[evader.x, evader.y]])
        vor = Voronoi(points)
        region_index = vor.point_region[-1]
        region = vor.regions[region_index]

        if not region or -1 in region:
            return np.random.randn(2)

        polygon = np.array([vor.vertices[i] for i in region])
        best_score = -np.inf
        best_point = np.array([evader.x, evader.y])

        for point in polygon:
            score = min(np.linalg.norm(point - np.array([p.x, p.y])) for p in pursuers)
            if score > best_score:
                best_score = score
                best_point = point

        return best_point - np.array([evader.x, evader.y])

class FleeClosestEvade(EvaderStrategy):
    def compute_target(self, evader, pursuers):
        if not pursuers:
            return np.zeros(2)

        closest = min(pursuers, key=lambda p: np.hypot(evader.x - p.x, evader.y - p.y))
        return np.array([evader.x - closest.x, evader.y - closest.y])
    

class VoronoiEvasion(EvaderStrategy):
    def compute_target(self, evader, pursuers):
        if len(pursuers) <= 3:
            if pursuers:
                closest_pursuer = min(pursuers, key=lambda p: np.hypot(p.x - self.x, p.y - self.y))
                direction = np.array([self.x - closest_pursuer.x, self.y - closest_pursuer.y])
            return direction
        
        points = np.array([[p.x, p.y] for p in pursuers] + [[self.x, self.y]])
        vor = Voronoi(points)
        region_index = vor.point_region[-1]  # Voronoi region index for evader
        region = vor.regions[region_index]
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

class VoronoiEvasion2(EvaderStrategy):
    pass


