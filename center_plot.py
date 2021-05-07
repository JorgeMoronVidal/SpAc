import pandas as pd
import matplotlib.pyplot as plt 
import os
indir = 'Python'
fig, ax = plt.subplots()
for file in os.listdir(indir):
    print(file)
    table = pd.read_csv(indir + '/' + file,',')
    x = table['x']
    y = table['y']
    label = table['i']
    ax.plot(table['x'], table['y'],'-o')
    for i, txt in enumerate(label):
        ax.annotate(txt, (x[i], y[i]))
ax.set_aspect('equal', 'box')
plt.show()