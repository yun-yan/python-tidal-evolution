import matplotlib.pyplot as plt
from matplotlib import patches


plt.rcParams['figure.figsize'] = (4.0, 4.0)

circle1 = plt.Circle((0, 0), 0.2, color='r')
circle2 = plt.Circle((0.5, 0.5), 0.1, color='blue')
circle3 = plt.Circle((1, 1), 0.2, color='g', clip_on=False)
ellipse=patches.Ellipse((0.5, 0.5), 0.2, 0.2,
                     angle=-90, linewidth=2, fill=False, zorder=2)

fig, ax = plt.subplots() # note we must use plt.subplots, not plt.subplot
# (or if you have an existing figure)
# fig = plt.gcf()
# ax = fig.gca()

ax.add_artist(circle1)
ax.add_artist(circle2)
ax.add_artist(circle3)
ax.add_artist(ellipse)

fig.savefig('ellipse_compare')
