import json
import matplotlib.pyplot as plt
import numpy as np

# This code reads the results data from a file created through the code in "run_multiple.py" and presents it as a python list

# Specify the name (or path) of the file to retrieve data from
filename = "simulation_results.txt"
# Which line to retrieve the data from, 1 means first line, and so forth
active_line = 11 # Set to 0 for last line

with open(filename,"r") as f:
    lines = f.readlines()
data = json.loads(lines[active_line-1])["Results"]
# The data should norw be in the variable "data" as a list, and whatever you want can be done to it, for example:

cleandata = [d for d in data if d < 3000]

print(len(cleandata))

print(data)
print("Mean: ", np.mean(cleandata))
plt.hist(cleandata)
plt.show()