import json
import matplotlib.pyplot as plt

# This code reads the results data from a file created through the code in "run_multiple.py" and presents it as a python list

# Specify the name (or path) of the file to retrieve data from
filename = "simulation_results.txt"
# Which line to retrieve the data from, 1 means first line, and so forth
active_line = 5

with open(filename,"r") as f:
    lines = f.readlines()
data = json.loads(lines[active_line-1])["Results"]
# The data should norw be in the variable "data" as a list, and whatever you want can be done to it, for example:

print(data)
plt.hist(data)
plt.show()