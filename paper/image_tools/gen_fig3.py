# -- Figure 3 : Silhouette Model -- 

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style('whitegrid')

points = np.array([
    [3.3, 1.0],
    [4.1, 0.4],
    [4.6, 0.5],
    [4.2, 2.9],
    [4.5, 1.4]
])

def draw_triangle(p, m=1, c='k'):
    x, y = p[0], p[1]
    b = y - m*x
    x2 = 0
    y2 = b
    plt.plot([x, x2], [y, y2], '-', c=c)
    plt.plot([x, 2*x], [y, y2], '-', c=c)

fig = plt.figure(figsize=(7, 2.5))

for point in points:
    if point[1] == 2.9:
        draw_triangle([4.2, 2.9], m=1, c='orange')    
    else:
        draw_triangle(point, m=1, c='k')

plt.scatter(
    points[:, 0], 
    points[:, 1], 
    marker='o', 
    edgecolor='black', 
    linewidth='1', 
    s=40,
    facecolor='orange',
    zorder=10
)

plt.plot([0, 10], [0, 0], 'k')
plt.plot([7, 10], [0, 0], 'orange')
plt.plot([0, 1.3], [0, 0], 'orange')
plt.xlabel('(Birth+Death)/2', fontsize=18)
plt.ylabel('(Death-Birth)/2', fontsize=18)
plt.xlim([0, 10])
plt.ylim([0, 3.5])
plt.tick_params(labelsize=18)
plt.tight_layout()
plt.savefig('figure_3_silh.pdf')
# plt.show()
