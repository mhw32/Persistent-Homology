# -- Figure 3 : Silhouette Model -- 

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style('whitegrid')

def draw_triangle(p, m=1, c1='k', c2='k'):
    x, y = p[0], p[1]
    b = y - m*x
    x2 = 0
    y2 = b
    plt.plot([x, x2], [y, y2], '-', c=c1)
    plt.plot([x, 2*x], [y, y2], '-', c=c2)

fig = plt.figure(figsize=(7, 2.5))

draw_triangle([3.3, 1.0], m=1, c1='#1e90ff')
draw_triangle([4.1, 0.4], m=1)
draw_triangle([4.6, 0.5], m=1)
draw_triangle([4.5, 1.4], m=1, c2='#1e90ff')
draw_triangle([4.2, 2.9], m=1, c1='orange', c2='orange')    

plt.scatter(
    [4.1, 4.6], 
    [0.4, 0.5], 
    marker='o', 
    edgecolor='black', 
    linewidth='1', 
    s=40,
    facecolor='grey',
    zorder=10
)

plt.scatter(
    [4.2], 
    [2.9], 
    marker='o', 
    edgecolor='black', 
    linewidth='1', 
    s=40,
    facecolor='orange',
    zorder=10
)

plt.scatter(
    [3.3, 4.5], 
    [1.0, 1.4], 
    marker='o', 
    edgecolor='black', 
    linewidth='1', 
    s=40,
    facecolor='#1e90ff',
    zorder=10
)

plt.plot([3.3, 3.7], [1.0, 0.6], '#1e90ff')
plt.plot([3.7, 4.5], [0.6, 1.4], '#1e90ff')
plt.plot([0, 10], [0, 0], 'k')
plt.plot([7, 10], [0, 0], 'orange')
plt.plot([0, 1.3], [0, 0], 'orange')
plt.plot([1.3, 2.3], [0, 0], '#1e90ff')
plt.plot([5.9, 7], [0, 0], '#1e90ff')
plt.xlabel('(Birth+Death)/2', fontsize=18)
plt.ylabel('(Death-Birth)/2', fontsize=18)
plt.xlim([0, 10])
plt.ylim([0, 3.5])
plt.tick_params(labelsize=16)
plt.tight_layout()
plt.savefig('figure_3_silh.pdf')
# plt.show()
