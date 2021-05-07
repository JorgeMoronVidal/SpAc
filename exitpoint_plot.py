import pandas as pd
import matplotlib.pyplot as plt 
import os
file = "Debug/ExitPoints.csv"
print(file)
fig, ax = plt.subplots()
table = pd.read_csv(file)
x = table['x']
y = table['y']
ax.scatter(table['x'], table['y'])
ax.set_aspect('equal', 'box')
plt.show()