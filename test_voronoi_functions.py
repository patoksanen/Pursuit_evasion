from ep_classes import *
import voronoi_functions
import scipy.spatial
import matplotlib.pyplot as plt

# Generera punkter att testa med
e = Evader(10,10,1)
p1 = Pursuer(50,50,1,e)
p2 = Pursuer(60,70,1,e)

v = voronoi_functions.getVoronoiCells([e,p1,p2])
print("Punkter: ",v.points)
print("Antal punkter: ",len(v.points))
# Märk att vi har tre punkter men får 15 punkter tillbaka (3 + 3*4 = 15). Det är för att vi speglar punkterna i gränslinjerna.
# Vi behöver bara bry oss om de första 3 punkterna, som i ordning motsvarar e, p1, p2
print("Viktiga punkter:", v.points[:3])
# Motsvarar e, p1, p2

# Plot
scipy.spatial.voronoi_plot_2d(v)
plt.show()

# Area för p1:s voronoi-cell
vertices = [ v.vertices[i] for i in v.regions[v.point_region[0]] ]
print("Vertices: ", vertices)
print("Area p1: ", voronoi_functions.polyarea(vertices))

# Vilken area har alla voronoi-celler sammanlagt?
A_sum = 0
for a in range(3):
    vtcs = [ v.vertices[i] for i in v.regions[v.point_region[a]] ]
    A_sum += voronoi_functions.polyarea(vtcs)
print("Total area: ", A_sum)
# Den totala arean är lite större än 10000 eftersom vi har en buffer i funktionen getVoronoiCells()
# Jämför med koden nedan där "buffer=0" för getVoronoiCells()

v = voronoi_functions.getVoronoiCells([e,p1,p2], buffer=0)
A_sum = 0
for a in range(3):
    vtcs = [ v.vertices[i] for i in v.regions[v.point_region[a]] ]
    A_sum += voronoi_functions.polyarea(vtcs)
print("Total area (utan buffer): ", A_sum)
