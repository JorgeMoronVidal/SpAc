import pandas as pd
import matplotlib.pyplot as plt 
import numpy as np
import os
indir = 'Python'
fig, ax = plt.subplots()
for file in os.listdir(indir):
    print(file)
    table = pd.read_csv(indir + '/' + file,',')
    x = table['x']
    y = table['y']
    label = table['i']
    if '12' in file:
        ax.plot(table['x'], table['y'],'-+',color = 'red')
    elif '_0_' in file:
        ax.plot(table['x'], table['y'],'-+',color = 'red')
    elif '_0.' in file:
        ax.plot(table['x'], table['y'],'-+',color = 'red')
    else:
        ax.plot(table['x'], table['y'],'-+', color = 'blue')
    #for i, txt in enumerate(label):
        #ax.annotate(txt, (x[i], y[i]))
ax.set_aspect('equal', 'box')
x_dom = np.array([-1, -1, 1, 1, -1]) * 6
y_dom = np.array([-1, 1, 1, -1, -1]) * 6
ax.plot(x_dom, y_dom, color ='black',linewidth=4)
plt.title('[6x6] SpAc Domain Decomposition')
plt.xlabel('x')
plt.ylabel('y')
plt.show()